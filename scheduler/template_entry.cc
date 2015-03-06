/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 /*
  * This is TemplateEntry module to hold and instantiate the templates.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/template_entry.h"
#include "scheduler/job_manager.h"
#include "shared/helpers.h"

using namespace nimbus; // NOLINT

#define MAX_DEPTH 100
#define ACTIVE_INSTANTIATION true

TemplateEntry::TemplateEntry(std::string template_name) {
  finalized_ = false;
  template_name_ = template_name;
  // TODO(omidm): currently we do not support future job id in templates!
  future_job_id_ptr_ = boost::shared_ptr<job_id_t>(new job_id_t(0));
}

TemplateEntry::~TemplateEntry() {
  // Used shared ptr for allocated pointers. -omidm

  // Clean the access pattern meta data
  AccessIndex::iterator iter = access_pattern_.begin();
  for (; iter != access_pattern_.end(); ++iter) {
    VersionIndex::iterator it = iter->second->begin();
    for (; it != iter->second->end(); ++it) {
      delete it->second;
    }
    delete iter->second;
  }
}

bool TemplateEntry::finalized() {
  return finalized_;
}

size_t TemplateEntry::compute_jobs_num() {
  return compute_jobs_.size();
}

std::string TemplateEntry::template_name() {
  return template_name_;
}

boost::shared_ptr<VersionMap> TemplateEntry::vmap_base() const {
  return vmap_base_;
}

void TemplateEntry::set_vmap_base(boost::shared_ptr<VersionMap> vmap_base) {
  vmap_base_ = vmap_base;
}

bool TemplateEntry::Finalize() {
  if (finalized_) {
    dbg(DBG_WARN, "WARNING: template has been already finalized!\n");
    return true;
  }

  if (job_id_ptrs_map_.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: referenced jobs are not equal to defined jobs!\n");
    return false;
  }

  assert(entry_list_.size() == job_id_ptrs_.size());

  CompleteParentJobIndices();

  CompleteBreadthFirstSearch();

  finalized_ = true;
  return true;
}

void TemplateEntry::CompleteParentJobIndices() {
  size_t index = 0;
  parent_job_indices_.clear();
  TemplateJobEntryVector::iterator iter = compute_jobs_.begin();
  for (; iter != compute_jobs_.end(); ++iter) {
    if (!(*iter)->sterile()) {
      parent_job_indices_.push_back(index);
    }
    ++index;
  }
}

void TemplateEntry::CompleteBreadthFirstSearch() {
  std::cout << "OMID STARTED: " << template_name_ << std::endl;
  size_t depth = 0;
  while (assign_ordered_indices_.size() != compute_jobs_.size()) {
    std::cout << "OMID depth: " << depth << std::endl;
    std::cout << "OMID compute size: " << compute_jobs_.size() << std::endl;
    std::cout << "OMID assigned ordered size: " << assign_ordered_indices_.size() << std::endl;
    std::cout << "OMID traverese queue size: " << traverse_queue_.size() << std::endl;
    assert(depth++ < MAX_DEPTH);
    assert(traverse_queue_.size() != 0);
    std::list<Vertex<TemplateJobEntry, job_id_t>*> temp_traverse_queue;
    assign_batch_mark_indices_.insert(assign_ordered_indices_.back());

    std::list<Vertex<TemplateJobEntry, job_id_t>*>::iterator iter;
    for (iter = traverse_queue_.begin(); iter != traverse_queue_.end(); ++iter) {
      typename Edge<TemplateJobEntry, job_id_t>::Iter it;
      for (it = (*iter)->outgoing_edges()->begin(); it != (*iter)->outgoing_edges()->end(); ++it) {
        Vertex<TemplateJobEntry, job_id_t>* start = it->second->start_vertex();
        Vertex<TemplateJobEntry, job_id_t>* end = it->second->end_vertex();
        bool removed_edge = job_graph_.RemoveEdge(start, end);
        assert(removed_edge);

        if (end->incoming_edges()->size() == 0) {
          temp_traverse_queue.push_back(end);
          assign_ordered_indices_.push_back(end->entry()->index());
	  std::cout << "OMID for: " << template_name_ << " NEXT:  " << end->entry()->job_name() << " " << end->entry()->job_id() << std::endl; // NOLINT
        }
      }
    }

    traverse_queue_.swap(temp_traverse_queue);
  }

  assign_batch_mark_indices_.insert(assign_ordered_indices_.back());
  last_assign_index_ = assign_ordered_indices_.back();
  assert(assign_ordered_indices_.size() == compute_jobs_.size());
}

bool TemplateEntry::CleanPartiallyFilledTemplate() {
  if (finalized_) {
    dbg(DBG_ERROR, "ERROR: template has been finalized and cannot get cleaned!\n");
    return false;
  }

  entry_list_.clear();
  compute_jobs_.clear();
  job_id_ptrs_.clear();
  job_id_ptrs_map_.clear();
  return true;
}


bool TemplateEntry::LoadBeforeSet(IDSet<job_id_t>* before_set,
                                  const size_t& index,
                                  const std::vector<job_id_t>& inner_job_ids,
                                  const std::vector<job_id_t>& outer_job_ids) {
  if (!finalized_) {
    dbg(DBG_ERROR, "ERROR: template has NOT been finalized and cannot get instantiated!\n");
    return false;
  }

  assert(entry_list_.size() == job_id_ptrs_.size());

  if (inner_job_ids.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided ids does not match the required ids!\n");
    return false;
  }

  assert(index < inner_job_ids.size());

  // Set the job_id pointers to the new values.
  size_t idx = 0;
  PtrList::iterator piter = job_id_ptrs_.begin();
  for (; piter != job_id_ptrs_.end(); ++piter) {
    *(*piter) = inner_job_ids[idx];
    ++idx;
  }

  idx = 0;
  EntryList::iterator iter = entry_list_.begin();
  for (; idx < index; ++idx) {
    ++iter;
  }

  {
    // TODO(omidm) Does accesing a field in class make a copy?
    PtrSet::iterator it = iter->before_set_ptrs_.begin();
    for (; it != iter->before_set_ptrs_.end(); ++it) {
      before_set->insert(*(*it));
    }
  }

  return true;
}







bool TemplateEntry::Instantiate(JobManager *job_manager,
                                const std::vector<job_id_t>& inner_job_ids,
                                const std::vector<job_id_t>& outer_job_ids,
                                const std::vector<Parameter>& parameters,
                                const job_id_t& parent_job_id) {
  if (!finalized_) {
    dbg(DBG_ERROR, "ERROR: template has NOT been finalized and cannot get instantiated!\n");
    return false;
  }

  assert(entry_list_.size() == job_id_ptrs_.size());

  if (inner_job_ids.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided ids does not match the required ids!\n");
    return false;
  }

  if (parameters.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided parameters does not match the required ids!\n");
    return false;
  }


  // Set the job_id pointers to the new values.
  size_t index = 0;
  PtrList::iterator piter = job_id_ptrs_.begin();
  for (; piter != job_id_ptrs_.end(); ++piter) {
    *(*piter) = inner_job_ids[index];
    ++index;
  }

  index = 0;
  EntryList::iterator iter = entry_list_.begin();
  for (; iter != entry_list_.end(); ++iter) {
    IDSet<job_id_t> before_set;
    {
      // TODO(omidm) Does accesing a field in class make a copy?
      PtrSet::iterator it = iter->before_set_ptrs_.begin();
      for (; it != iter->before_set_ptrs_.end(); ++it) {
        before_set.insert(*(*it));
      }
    }

    IDSet<job_id_t> after_set;
    {
      // TODO(omidm) Does accesing a field in class make a copy?
      PtrSet::iterator it = iter->after_set_ptrs_.begin();
      for (; it != iter->after_set_ptrs_.end(); ++it) {
        after_set.insert(*(*it));
      }
    }

    job_manager->AddComputeJobEntry(iter->job_name_,
                                    *(job_id_ptrs_[index]),
                                    iter->read_set_,
                                    iter->write_set_,
                                    before_set,
                                    after_set,
                                    parent_job_id,
                                    *(iter->future_job_id_ptr_),
                                    iter->sterile_,
                                    iter->region_,
                                    parameters[index]);
    ++index;
  }

  return true;
}

size_t TemplateEntry::GetParentJobIndices(std::list<size_t>* list) {
  assert(finalized_);
  list->clear();
  *list = parent_job_indices_;
  return list->size();
}


TemplateJobEntry* TemplateEntry::GetJobAtIndex(size_t index) {
  if (index >= compute_jobs_.size()) {
    return NULL;
  }

  return compute_jobs_[index];
}


bool TemplateEntry::GetComplexJobEntry(ComplexJobEntry*& complex_job,
                                       const job_id_t& job_id,
                                       const job_id_t& parent_job_id,
                                       const std::vector<job_id_t>& inner_job_ids,
                                       const std::vector<job_id_t>& outer_job_ids,
                                       const std::vector<Parameter>& parameters) {
  if (!finalized_) {
    dbg(DBG_ERROR, "ERROR: template has NOT been finalized and cannot get instantiated!\n");
    complex_job = NULL;
    return false;
  }

  if (inner_job_ids.size() != compute_jobs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided ids does not match the required ids!\n");
    complex_job = NULL;
    return false;
  }

  if (parameters.size() != compute_jobs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided parameters does not match the required ids!\n");
    complex_job = NULL;
    return false;
  }

  complex_job = new ComplexJobEntry(job_id,
                                    parent_job_id,
                                    this,
                                    inner_job_ids,
                                    outer_job_ids,
                                    parameters);

  return true;
}

TemplateJobEntry* TemplateEntry::AddComputeJob(const std::string& job_name,
                                               const job_id_t& job_id,
                                               const IDSet<logical_data_id_t>& read_set,
                                               const IDSet<logical_data_id_t>& write_set,
                                               const IDSet<job_id_t>& before_set,
                                               const IDSet<job_id_t>& after_set,
                                               const job_id_t& parent_job_id,
                                               const job_id_t& future_job_id,
                                               const bool& sterile,
                                               const GeometricRegion& region) {
  if (finalized_) {
    dbg(DBG_ERROR, "ERROR: template has been finalized and cannot add compute job!\n");
    return NULL;
  }

  TemplateJobEntry *job =
    new TemplateJobEntry(job_name,
                         job_id,
                         compute_jobs_.size(),
                         read_set,
                         write_set,
                         before_set,
                         after_set,
                         sterile,
                         region,
                         this);

  compute_jobs_.push_back(job);
  AddTemplateJobEntryToJobGraph(job);

  if (!ACTIVE_INSTANTIATION) {
    return job;
  } else {
    // TODO(omidm): the rest of theis function could be ignored if we don't need
    // the Instantiation function to work.
    boost::shared_ptr<job_id_t> job_id_ptr;
    {
      PtrMap::iterator iter = job_id_ptrs_map_.find(job_id);
      if (iter == job_id_ptrs_map_.end()) {
        job_id_ptr = boost::shared_ptr<job_id_t>(new job_id_t(job_id));
        job_id_ptrs_map_[job_id] = job_id_ptr;
      } else {
        job_id_ptr = iter->second;
      }
    }

    PtrSet before_set_ptrs;
    {
      IDSet<job_id_t>::IDSetIter it = before_set.begin();
      for (; it != before_set.end(); ++it) {
        boost::shared_ptr<job_id_t> ptr;
        PtrMap::iterator iter = job_id_ptrs_map_.find(*it);
        if (iter == job_id_ptrs_map_.end()) {
          ptr = boost::shared_ptr<job_id_t>(new job_id_t(*it));
          job_id_ptrs_map_[*it] = ptr;
        } else {
          ptr = iter->second;
        }
        before_set_ptrs.insert(ptr);
      }
    }

    PtrSet after_set_ptrs;
    {
      IDSet<job_id_t>::IDSetIter it = after_set.begin();
      for (; it != after_set.end(); ++it) {
        boost::shared_ptr<job_id_t> ptr;
        PtrMap::iterator iter = job_id_ptrs_map_.find(*it);
        if (iter == job_id_ptrs_map_.end()) {
          ptr = boost::shared_ptr<job_id_t>(new job_id_t(*it));
          job_id_ptrs_map_[*it] = ptr;
        } else {
          ptr = iter->second;
        }
        after_set_ptrs.insert(ptr);
      }
    }

    TemplateComputeJobEntry entry(job_name,
        job_id_ptr,
        read_set,
        write_set,
        before_set_ptrs,
        after_set_ptrs,
        future_job_id_ptr_,
        sterile,
        region);

    entry_list_.push_back(entry);
    job_id_ptrs_.push_back(job_id_ptr);
    // End of block
  }

  return job;
}

bool TemplateEntry::AddExplicitCopyJob() {
  dbg(DBG_ERROR, "ERROR: explicit copy jobs from application are not supported yet!.\n");
  exit(-1);
  return false;
}


bool TemplateEntry::AddTemplateJobEntryToJobGraph(TemplateJobEntry *job) {
  Vertex<TemplateJobEntry, job_id_t> *vertex;
  bool added_node = job_graph_.AddVertex(job->job_id(), job, &vertex);
  assert(added_node);
  std::cout << "OMID for: " << template_name_ << " ADD:  " << job->job_name() << " " << job->job_id() << std::endl; // NOLINT
  std::cout << "OMID BEFORESET: " << job->before_set_p()->ToNetworkData() << std::endl; // NOLINT

  if (job->before_set_p()->size() == 0) {
    traverse_queue_.push_back(vertex);
    assign_ordered_indices_.push_back(job->index());
    std::cout << "OMID for: " << template_name_ << " INIT:  " << job->job_name() << " " << job->job_id() << std::endl; // NOLINT
  } else {
    IDSet<job_id_t>::ConstIter it;
    for (it = job->before_set_p()->begin(); it != job->before_set_p()->end(); ++it) {
      bool added_edge = job_graph_.AddEdge(*it, job->job_id());
      assert(added_edge);
    }
  }

  return true;
}

bool TemplateEntry::InitializeCursor(ComplexJobEntry::Cursor* cursor) {
  assert(finalized_);
  size_t pivot = 0;
  size_t index = assign_ordered_indices_[pivot];

  cursor->set_index(index);
  cursor->set_pivot(pivot);
  if (index == last_assign_index_) {
    cursor->set_state(ComplexJobEntry::Cursor::END_ALL);
  } else if (assign_batch_mark_indices_.count(index)) {
    cursor->set_state(ComplexJobEntry::Cursor::END_BATCH);
  } else {
    cursor->set_state(ComplexJobEntry::Cursor::MID_BATCH);
  }

  return true;
}

bool TemplateEntry::AdvanceCursorForAssignment(ComplexJobEntry::Cursor* cursor) {
  assert(finalized_);
  size_t pivot = cursor->pivot() + 1;
  if (pivot >= assign_ordered_indices_.size()) {
    return false;
  }

  size_t index = assign_ordered_indices_[pivot];

  cursor->set_index(index);
  cursor->set_pivot(pivot);
  if (index == last_assign_index_) {
    cursor->set_state(ComplexJobEntry::Cursor::END_ALL);
  } else if (assign_batch_mark_indices_.count(index)) {
    cursor->set_state(ComplexJobEntry::Cursor::END_BATCH);
  } else {
    cursor->set_state(ComplexJobEntry::Cursor::MID_BATCH);
  }

  return true;
}

void TemplateEntry::AddToAccessPattern(const logical_data_id_t& ldid,
                                       const data_version_t& diff_version,
                                       const size_t& job_index) {
  boost::unique_lock<boost::mutex> lock(access_pattern_mutex_);

  AccessIndex::iterator iter = access_pattern_.find(ldid);
  if (iter != access_pattern_.end()) {
    VersionIndex::iterator it = iter->second->find(diff_version);
    if (it != iter->second->end()) {
      it->second->push_back(job_index);
    } else {
      Bucket *b = new Bucket();
      b->push_back(job_index);
      iter->second->insert(make_pair(diff_version, b));
    }
  } else {
    VersionIndex *vi = new VersionIndex();
    Bucket *b = new Bucket();
    b->push_back(job_index);
    vi->insert(make_pair(diff_version, b));
    access_pattern_.insert(make_pair(ldid, vi));
  }
}

size_t TemplateEntry::QueryAccessPattern(const logical_data_id_t& ldid,
                                         const data_version_t& diff_version,
                                         std::list<size_t>* indices) {
  boost::unique_lock<boost::mutex> lock(access_pattern_mutex_);

  size_t count = 0;
  indices->clear();

  AccessIndex::iterator iter = access_pattern_.find(ldid);
  if (iter != access_pattern_.end()) {
    VersionIndex::iterator it = iter->second->find(diff_version);
    if (it != iter->second->end()) {
      Bucket::iterator i = it->second->begin();
      for (; i != it->second->end(); ++i) {
        indices->push_back(*i);
        ++count;
      }
    }
  }

  return count;
}

std::string TemplateEntry::ProduceBindingRecordName(size_t load_balancing_tag,
                                                    const std::string& grand_parent_name) {
  std::string key = int2string(load_balancing_tag) + "-" + template_name_ + "-" + grand_parent_name;
  return key;
}

bool TemplateEntry::AddBindingRecord(size_t load_balancing_tag,
                                     const std::string& grand_parent_name,
                                     const std::vector<job_id_t>& compute_job_ids,
                                     BindingTemplate*& binding_template) {
  std::string record_name = ProduceBindingRecordName(load_balancing_tag, grand_parent_name);

  BindingMap::iterator iter = binding_records_.find(record_name);
  if (iter != binding_records_.end()) {
    binding_template = NULL;
    return false;
  }

  BindingTemplate* bt = new BindingTemplate(record_name,
                                            compute_job_ids,
                                            this);
  binding_records_[record_name] = bt;
  binding_template = bt;
  return true;
}

bool TemplateEntry::QueryBindingRecord(size_t load_balancing_tag,
                                       const std::string& grand_parent_name,
                                       BindingTemplate*& binding_template) {
  std::string record_name = ProduceBindingRecordName(load_balancing_tag, grand_parent_name);

  BindingMap::iterator iter = binding_records_.find(record_name);
  if (iter == binding_records_.end()) {
    binding_template = NULL;
    return false;
  }

  binding_template = iter->second;
  return true;
}





