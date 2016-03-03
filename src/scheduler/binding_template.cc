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
  * This is BindingTemplate module to hold and instantiate copy and compute
  * commands sent to workers to bind physical data to a batch of jobs in a
  * template complex job.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/scheduler/binding_template.h"
#include "src/scheduler/template_entry.h"
#include "src/scheduler/job_manager.h"

#define MEGA_RCR_COUNT_THREASHOLD 2

using namespace nimbus; // NOLINT

BindingTemplate::BindingTemplate(const std::string& record_name,
                                 const std::vector<job_id_t>& compute_job_ids,
                                 TemplateEntry *template_entry,
                                 bool worker_template_active,
                                 bool mega_rcr_job_active) {
  finalized_ = false;
  established_command_template_ = false;
  record_name_ = record_name;
  template_entry_ = template_entry;
  future_job_id_ptr_ = JobIdPtr(new job_id_t(0));
  worker_template_active_ = worker_template_active;
  mega_rcr_job_active_ = mega_rcr_job_active;

  std::vector<job_id_t>::const_iterator iter = compute_job_ids.begin();
  for (; iter != compute_job_ids.end(); ++iter) {
    JobIdPtr job_id_ptr = JobIdPtr(new job_id_t(*iter));
    compute_job_id_map_[*iter] = job_id_ptr;
  }
}

BindingTemplate::~BindingTemplate() {
}

bool BindingTemplate::finalized() const {
  boost::unique_lock<boost::mutex> lock(mutex_);
  return finalized_;
}

size_t BindingTemplate::copy_job_num() const {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  return copy_job_id_list_.size();
}

size_t BindingTemplate::compute_job_num() const {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  return compute_job_id_list_.size();
}

std::string BindingTemplate::record_name() const {
  return record_name_;
}

const BindingTemplate::PatternMetaData* BindingTemplate::patterns_meta_data_p() const {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  return &patterns_meta_data_;
}

const BindingTemplate::PatternList* BindingTemplate::entry_pattern_list_p() const {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  return &entry_pattern_list_;
}

const BindingTemplate::PatternList* BindingTemplate::end_pattern_list_p() const {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  return &end_pattern_list_;
}



void
BindingTemplate::GetRequiredUpdatesForCascading(boost::shared_ptr<VersionList> vlist_write_diff,
                                                ConstPatternList *patterns,
                                                std::vector<data_version_t> *versions_diff) {
  patterns->clear();
  versions_diff->clear();

  VersionMap vmap_diff;
  {
    VersionList::const_iterator iter = vlist_write_diff->begin();
    for (; iter != vlist_write_diff->end(); ++iter) {
      vmap_diff.set_entry(iter->first, iter->second);
    }
  }

  PatternMap::iterator iter = entry_pattern_map_.begin();
  for (; iter != entry_pattern_map_.end(); ++iter) {
    if (iter->second->version_type_ == SCRATCH) {
      dbg(DBG_ERROR, "ERROR: entry pattern cannot have scratch type!\n");
      assert(false);
    }
    if (iter->second->version_type_ == WILD_CARD) {
      continue;
    }

    if (iter->second->version_diff_from_base_ != 0) {
      dbg(DBG_ERROR, "ERROR: regular entry pattern cannot have non zero version diff!\n");
      assert(false);
    }

    logical_data_id_t ldid = iter->second->ldid_;
    physical_data_id_t pdid = iter->second->pdid_;

    data_version_t write_diff_version;
    if (!vmap_diff.query_entry(ldid, &write_diff_version)) {
      write_diff_version = 0;
    }

    PatternMap::iterator it = end_pattern_map_.find(pdid);
    assert(it != end_pattern_map_.end());
    if (it->second->version_type_ == WILD_CARD) {
      dbg(DBG_ERROR, "ERROR: ennd pattern cannot have wild_card type!\n");
      assert(false);
    }

    data_version_t v_diff_needed = iter->second->version_diff_from_base_ + write_diff_version;
    if (it->second->version_type_ == SCRATCH) {
      patterns->push_back(it->second);
      versions_diff->push_back(v_diff_needed);
    } else {
      data_version_t v_diff_now = it->second->version_diff_from_base_;

      if (v_diff_now != v_diff_needed) {
        assert(v_diff_now < v_diff_needed);
        patterns->push_back(it->second);
        versions_diff->push_back(v_diff_needed);
      }
    }
  }
}


bool BindingTemplate::Finalize(const std::vector<job_id_t>& compute_job_ids) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(!finalized_);
  assert(compute_job_id_map_.size() == template_entry_->compute_jobs_num());
  assert(compute_job_id_map_.size() == compute_job_ids.size());
  assert(compute_job_to_command_map_.size() == compute_job_ids.size());


  if (mega_rcr_job_active_) {
    std::cout << "OMID Template: "
              << template_entry_->template_name()
              << " RCR num: "
              << candidate_rcr_.size()
              << std::endl;

    // TODO(omidm): have to check for all jobs, not just compute!!
    std::map<job_id_t, job_id_t> rcr_compute_pairs;
    {
      CommandTemplateList::iterator iter = command_templates_.begin();
      for (; iter != command_templates_.end(); ++iter) {
        JobIdPtrSet bs = (*iter)->before_set_ptr_;
        JobIdPtrSet::iterator i = bs.begin();
        for (; i != bs.end(); ++i) {
          job_id_t rcr = *(*i);
          if (candidate_rcr_.contains(rcr)) {
            candidate_rcr_.remove(rcr);
            assert(rcr_compute_pairs.count(rcr) == 0);
            if ((*iter)->type_ == COMPUTE) {
              rcr_compute_pairs[rcr] = *((*iter)->job_id_ptr_);
            }
          } else if (rcr_compute_pairs.count(rcr)) {
            rcr_compute_pairs.erase(rcr);
          }
        }
      }
    }

    candidate_rcr_.clear();
    std::map<job_id_t, job_id_t> rcr_to_mega;
    std::map<job_id_t, job_id_t> compute_to_mega;
    std::map<job_id_t, IDSet<job_id_t> > compute_rcr_group;
    {
      std::map<job_id_t, job_id_t>::iterator iter = rcr_compute_pairs.begin();
      for (; iter != rcr_compute_pairs.end(); ++iter) {
        compute_rcr_group[iter->second].insert(iter->first);
      }
      std::map<job_id_t, IDSet<job_id_t> >::iterator it = compute_rcr_group.begin();
      for (; it != compute_rcr_group.end();) {
        IDSet<job_id_t> pool = it->second;
        if (pool.size() < MEGA_RCR_COUNT_THREASHOLD) {
          compute_rcr_group.erase(it++);
        } else {
          assert(pool.size() > 0);
          IDSet<job_id_t>::IDSetIter i = pool.begin();
          job_id_t mega_job_id = *i;
          compute_to_mega[it->first] = mega_job_id;
          for (; i != pool.end(); ++i) {
            candidate_rcr_.insert(*i);
            rcr_to_mega[*i] = mega_job_id;
          }
          ++it;
        }
      }
    }

    std::cout << "OMID Template: "
              << template_entry_->template_name()
              << "Mega Group Num: "
              << compute_rcr_group.size()
              << std::endl;

    // Checking that no ther job has candidates in before set!
    // If they do it could cause race!!!.
    {
      {
        std::map<job_id_t, CombineJobCommandTemplate*>::iterator iter =
          combine_job_to_command_map_.begin();
        for (; iter != combine_job_to_command_map_.end(); ++iter) {
          JobIdPtrSet bs = iter->second->before_set_ptr_;
          JobIdPtrSet::iterator i = bs.begin();
          for (; i != bs.end(); ++i) {
            job_id_t rcr = *(*i);
            if (candidate_rcr_.contains(rcr)) {
              assert(false);
              std::cout << "DANGER 0: "
                        << template_entry_->template_name()
                        << std::endl;
            }
          }
        }
      }
      {
        std::map<job_id_t, RemoteCopyReceiveCommandTemplate*>::iterator iter =
          rcr_job_to_command_map_.begin();
        for (; iter != rcr_job_to_command_map_.end(); ++iter) {
          JobIdPtrSet bs = iter->second->before_set_ptr_;
          JobIdPtrSet::iterator i = bs.begin();
          for (; i != bs.end(); ++i) {
            job_id_t rcr = *(*i);
            if (candidate_rcr_.contains(rcr)) {
              assert(false);
              std::cout << "DANGER 1: "
                        << template_entry_->template_name()
                        << std::endl;
            }
          }
        }
      }
      {
        std::map<job_id_t, RemoteCopySendCommandTemplate*>::iterator iter =
          rcs_job_to_command_map_.begin();
        for (; iter != rcs_job_to_command_map_.end(); ++iter) {
          JobIdPtrSet bs = iter->second->before_set_ptr_;
          JobIdPtrSet::iterator i = bs.begin();
          for (; i != bs.end(); ++i) {
            job_id_t rcr = *(*i);
            if (candidate_rcr_.contains(rcr)) {
              assert(false);
              std::cout << "DANGER 2: "
                        << template_entry_->template_name()
                        << std::endl;
            }
          }
        }
      }
      {
        std::map<job_id_t, LocalCopyCommandTemplate*>::iterator iter =
          lc_job_to_command_map_.begin();
        for (; iter != lc_job_to_command_map_.end(); ++iter) {
          JobIdPtrSet bs = iter->second->before_set_ptr_;
          JobIdPtrSet::iterator i = bs.begin();
          for (; i != bs.end(); ++i) {
            job_id_t rcr = *(*i);
            if (candidate_rcr_.contains(rcr)) {
              assert(false);
              std::cout << "DANGER 3: "
                        << template_entry_->template_name()
                        << std::endl;
            }
          }
        }
      }
    }




    // create the mega rcr commands
    {
      // TODO(omidm): complete
      std::map<job_id_t, IDSet<job_id_t> >::iterator iter = compute_rcr_group.begin();
      for (; iter != compute_rcr_group.end(); ++iter) {
        job_id_t mega_id = compute_to_mega[iter->first];
        JobIdPtr job_id_ptr = GetExistingCopyJobIdPtr(mega_id);
        JobIdPtrList rcr_job_id_ptrs;
        PhyIdPtrList to_physical_data_id_ptrs;

        std::map<job_id_t, ComputeJobCommandTemplate*>::iterator iter_c =
          compute_job_to_command_map_.find(iter->first);
        assert(iter_c != compute_job_to_command_map_.end());
        worker_id_t worker_id = iter_c->second->worker_id_;

        IDSet<job_id_t> rcr_group = iter->second;
        IDSet<job_id_t>::IDSetIter it = rcr_group.begin();
        for (; it != rcr_group.end(); ++it) {
          std::map<job_id_t, RemoteCopyReceiveCommandTemplate*>::iterator i =
            rcr_job_to_command_map_.find(*it);
          assert(i != rcr_job_to_command_map_.end());
          assert(i->second->worker_id_ == worker_id);
          rcr_job_id_ptrs.push_back(i->second->job_id_ptr_);
          to_physical_data_id_ptrs.push_back(i->second->to_physical_data_id_ptr_);
        }

        mega_rcr_job_to_command_map_[mega_id] =
          new MegaRCRCommandTemplate(job_id_ptr,
                                     rcr_job_id_ptrs,
                                     to_physical_data_id_ptrs,
                                     worker_id);
      }
    }

    // fix the before set of compute jobs with mega rcr and clean the obsolete ids
    {
      std::map<job_id_t, IDSet<job_id_t> >::iterator iter = compute_rcr_group.begin();
      for (; iter != compute_rcr_group.end(); ++iter) {
        std::map<job_id_t, ComputeJobCommandTemplate*>::iterator it =
          compute_job_to_command_map_.find(iter->first);
        assert(it != compute_job_to_command_map_.end());
        JobIdPtrSet bs = it->second->before_set_ptr_;
        JobIdPtrSet::iterator i = bs.begin();
        for (; i != bs.end();) {
          job_id_t rcr = *(*i);
          if (iter->second.contains(rcr)) {
            bs.erase(i++);
          } else {
            ++i;
          }
        }
        bs.insert(GetExistingCopyJobIdPtr(compute_to_mega[iter->first]));
        it->second->before_set_ptr_ = bs;
      }
    }

    // clean up obsolete rcr and insert mega rcr in command_templates_
    // Also set the meta_rcr_id for the sender job
    {
      // TODO(omidm): complete
      IDSet<job_id_t> inserted_mega_rcr;
      CommandTemplateList::iterator iter = command_templates_.begin();
      for (; iter != command_templates_.end();) {
        CommandTemplate *ct = *iter;
        job_id_t rcr, mega_id;
        RemoteCopySendCommandTemplate *rcsc;
        RemoteCopyReceiveCommandTemplate *rcrc;
        switch (ct->type_) {
          case COMPUTE:
          case COMBINE:
          case LC:
            ++iter;
            break;
          case RCR:
            rcrc = reinterpret_cast<RemoteCopyReceiveCommandTemplate*>(ct);
            rcr = *(rcrc->job_id_ptr_);
            if (candidate_rcr_.contains(rcr)) {
              delete rcrc;
              mega_id = rcr_to_mega[rcr];
              if (!inserted_mega_rcr.contains(mega_id)) {
                *iter = mega_rcr_job_to_command_map_[mega_id];
                inserted_mega_rcr.insert(mega_id);
                ++iter;
              } else {
                command_templates_.erase(iter++);
              }
            } else {
              ++iter;
            }
            break;
          case RCS:
            rcsc = reinterpret_cast<RemoteCopySendCommandTemplate*>(ct);
            rcr = *(rcsc->receive_job_id_ptr_);
            if (candidate_rcr_.contains(rcr)) {
              mega_id = rcr_to_mega[rcr];
              rcsc->mega_rcr_job_id_ptr_ = GetExistingCopyJobIdPtr(mega_id);
              if (worker_template_active_) {
                job_worker_map_[mega_id].insert(rcsc->worker_id_);
              }
            }
            ++iter;
            break;
          default:
            assert(false);
        }
      }
    }



    {
      std::map<job_id_t, IDSet<job_id_t> >::iterator iter = compute_rcr_group.begin();
      for (; iter != compute_rcr_group.end(); ++iter) {
        std::cout << "       OMID Group job_name: "
                  << compute_job_to_command_map_[iter->first]->job_name_
                  << " worker id: "
                  << compute_job_to_command_map_[iter->first]->worker_id_
                  << " RCR num: "
                  << iter->second.size()
                  << std::endl;
      }
    }

    //
  }

  // Set param index for the compute commands.
  {
    size_t idx = 0;
    std::vector<job_id_t>::const_iterator iter = compute_job_ids.begin();
    for (; iter != compute_job_ids.end(); ++iter) {
      std::map<job_id_t, ComputeJobCommandTemplate*>::iterator it =
        compute_job_to_command_map_.find(*iter);
      assert(it != compute_job_to_command_map_.end());
      it->second->param_index_ = idx;
      ++idx;
    }
  }

  compute_job_id_list_.clear();
  {
    std::vector<job_id_t>::const_iterator iter = compute_job_ids.begin();
    for (; iter != compute_job_ids.end(); ++iter) {
      JobIdPtrMap::iterator it = compute_job_id_map_.find(*iter);
      assert(it != compute_job_id_map_.end());
      compute_job_id_list_.push_back(it->second);
    }
  }

  if (worker_template_active_) {
    std::vector<job_id_t>::const_iterator iter = compute_job_ids.begin();
    for (; iter != compute_job_ids.end(); ++iter) {
      JobIdPtrMap::iterator it_p = compute_job_id_map_.find(*iter);
      assert(it_p != compute_job_id_map_.end());

      JobWorkerMap::iterator it = job_worker_map_.find(*iter);
      assert(it != job_worker_map_.end());
      std::set<worker_id_t>::iterator i = it->second.begin();
      for (; i != it->second.end(); ++i) {
        worker_job_ids_[*i].push_back(it_p->second);
      }
    }
  }


  copy_job_id_list_.clear();
  {
    JobIdPtrMap::iterator iter = copy_job_id_map_.begin();
    for (; iter != copy_job_id_map_.end(); ++iter) {
      copy_job_id_list_.push_back(iter->second);
    }
  }

  if (worker_template_active_) {
    JobIdPtrMap::iterator iter = copy_job_id_map_.begin();
    for (; iter != copy_job_id_map_.end(); ++iter) {
      JobWorkerMap::iterator it = job_worker_map_.find(iter->first);
      assert(it != job_worker_map_.end());
      std::set<worker_id_t>::iterator i = it->second.begin();
      for (; i != it->second.end(); ++i) {
        worker_job_ids_[*i].push_back(iter->second);
      }
    }
  }


  phy_id_list_.clear();
  end_pattern_list_.clear();
  entry_pattern_list_.clear();
  {
    PatternSorted::iterator iter = ordered_entry_patterns_.begin();
    for (; iter != ordered_entry_patterns_.end(); ++iter) {
      size_t regular_num = 0;
      size_t wild_card_num = 0;
      {
        PatternList::iterator it = iter->second.first->begin();
        for (; it != iter->second.first->end(); ++it) {
          ++regular_num;
          physical_data_id_t pdid = (*it)->pdid_;
          {
            PhyIdPtrMap::iterator i = phy_id_map_.find(pdid);
            assert(i != phy_id_map_.end());
            phy_id_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = end_pattern_map_.find(pdid);
            assert(i != end_pattern_map_.end());
            end_pattern_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = entry_pattern_map_.find(pdid);
            assert(i != entry_pattern_map_.end());
            entry_pattern_list_.push_back(i->second);
            assert(i->second == (*it));
          }
        }
      }
      {
        PatternList::iterator it = iter->second.second->begin();
        for (; it != iter->second.second->end(); ++it) {
          ++wild_card_num;
          physical_data_id_t pdid = (*it)->pdid_;
          {
            PhyIdPtrMap::iterator i = phy_id_map_.find(pdid);
            assert(i != phy_id_map_.end());
            phy_id_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = end_pattern_map_.find(pdid);
            assert(i != end_pattern_map_.end());
            end_pattern_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = entry_pattern_map_.find(pdid);
            assert(i != entry_pattern_map_.end());
            entry_pattern_list_.push_back(i->second);
            assert(i->second == (*it));
          }
        }
      }

      patterns_meta_data_.push_back(std::make_pair(regular_num, wild_card_num));
    }
  }

  if (worker_template_active_) {
    JobIdPtrMap::iterator iter = phy_id_map_.begin();
    for (; iter != phy_id_map_.end(); ++iter) {
      PhyWorkerMap::iterator it = phy_worker_map_.find(iter->first);
      assert(it != phy_worker_map_.end());
      std::set<worker_id_t>::iterator i = it->second.begin();
      for (; i != it->second.end(); ++i) {
        worker_phy_ids_[*i].push_back(iter->second);
      }
    }
  }

  if (worker_template_active_) {
    ComputeJobCommandTemplate *cc;
    CommandTemplateList::iterator iter = command_templates_.begin();
    for (; iter != command_templates_.end(); ++iter) {
      CommandTemplate *ct = *iter;
      if (ct->type_ == COMPUTE) {
        cc = reinterpret_cast<ComputeJobCommandTemplate*>(ct);
        worker_parameter_indices_[cc->worker_id_].push_back(cc->param_index_);
      }
    }
  }


  // Testing the integrity of computations!

  {
    PatternList::iterator iter = entry_pattern_list_.begin();
    for (; iter != entry_pattern_list_.end(); ++iter) {
      if ((*iter)->version_type_ == SCRATCH) {
        dbg(DBG_ERROR, "ERROR: entry pattern cannot have scratch type!\n");
        assert(false);
      }
      if ((*iter)->version_type_ == REGULAR) {
        if ((*iter)->version_diff_from_base_ != 0) {
          dbg(DBG_ERROR, "ERROR: regular entry pattern cannot have non zero version diff!\n");
          assert(false);
        }
      }
    }
  }


  {
    PatternList::iterator iter = end_pattern_list_.begin();
    for (; iter != end_pattern_list_.end(); ++iter) {
      if ((*iter)->version_type_ == WILD_CARD) {
        dbg(DBG_ERROR, "ERROR: end pattern cannot have wild-card type!\n");
        assert(false);
      }
    }
  }

  assert(copy_job_id_map_.size() == copy_job_id_list_.size());
  assert(compute_job_id_map_.size() == compute_job_id_list_.size());
  assert(phy_id_map_.size() == phy_id_list_.size());
  assert(phy_id_map_.size() == end_pattern_map_.size());
  assert(phy_id_map_.size() == entry_pattern_map_.size());
  assert(end_pattern_map_.size() == end_pattern_list_.size());
  assert(entry_pattern_map_.size() == entry_pattern_list_.size());
  assert(entry_pattern_list_.size() == end_pattern_list_.size());
  assert(ordered_entry_patterns_.size() == patterns_meta_data_.size());

  finalized_ = true;
  return true;
}

bool BindingTemplate::Instantiate(const std::vector<job_id_t>& compute_job_ids,
                                  const std::vector<Parameter>& parameters,
                                  const std::vector<job_id_t>& copy_job_ids,
                                  const std::vector<physical_data_id_t> *physical_ids,
                                  const ExtraDependency& extra_dependency,
                                  const template_id_t& template_generation_id,
                                  SchedulerServer *server) {
  Log log(Log::NO_FILE);
  log.log_StartTimer();

  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  assert(compute_job_ids.size() == compute_job_id_list_.size());
  assert(copy_job_ids.size() == copy_job_id_list_.size());
  assert(physical_ids->size() == phy_id_list_.size());

  {
    size_t idx = 0;
    JobIdPtrList::iterator iter = compute_job_id_list_.begin();
    for (; iter != compute_job_id_list_.end(); ++iter) {
      *(*iter) = compute_job_ids[idx];
      ++idx;
    }
  }

  {
    size_t idx = 0;
    JobIdPtrList::iterator iter = copy_job_id_list_.begin();
    for (; iter != copy_job_id_list_.end(); ++iter) {
      *(*iter) = copy_job_ids[idx];
      ++idx;
    }
  }

  {
    size_t idx = 0;
    PhyIdPtrList::iterator iter = phy_id_list_.begin();
    for (; iter != phy_id_list_.end(); ++iter) {
      *(*iter) = physical_ids->operator[](idx);
      ++idx;
    }
  }

  log.log_StopTimer();
  std::cout << "COMPLEX: Instantiate Loading: "
    << log.timer() << std::endl;


  if (established_command_template_) {
    SpawnCommandTemplateAtWorkers(parameters, extra_dependency, server, template_generation_id);
    return true;
  }


  if (worker_template_active_ && !established_command_template_) {
    SendCommandTemplateHeaderToWorkers(server);
  }

  CommandTemplateList::iterator iter = command_templates_.begin();
  for (; iter != command_templates_.end(); ++iter) {
    CommandTemplate *ct = *iter;
    ComputeJobCommandTemplate *cc = NULL;
    IDSet<job_id_t> ed;
    {
      ExtraDependency::const_iterator it =
        extra_dependency.find(ct->worker_id_);
      if (it != extra_dependency.end()) {
        ed = it->second;
      }
    }
    switch (ct->type_) {
      case COMPUTE:
        cc = reinterpret_cast<ComputeJobCommandTemplate*>(ct);
        SendComputeJobCommand(cc,
                              parameters[cc->param_index_],
                              ed,
                              server);
        break;
      case COMBINE:
        SendCombineJobCommand(reinterpret_cast<CombineJobCommandTemplate*>(ct),
                              ed,
                              server);
        break;
      case LC:
        SendLocalCopyCommand(reinterpret_cast<LocalCopyCommandTemplate*>(ct),
                             ed,
                             server);
        break;
      case RCS:
        SendRemoteCopySendCommand(reinterpret_cast<RemoteCopySendCommandTemplate*>(ct),
                                  ed,
                                  server);
        break;
      case RCR:
        SendRemoteCopyReceiveCommand(reinterpret_cast<RemoteCopyReceiveCommandTemplate*>(ct),
                                     ed,
                                     server);
        break;
      case MEGA_RCR:
        SendMegaRCRCommand(reinterpret_cast<MegaRCRCommandTemplate*>(ct),
                           ed,
                           server);
        break;
      default:
        assert(false);
    }
  }


  if (worker_template_active_ && !established_command_template_) {
    SendCommandTemplateFinalizeToWorkers(server);
    established_command_template_ = true;
  }


  return true;
}

void BindingTemplate::SendCommandTemplateHeaderToWorkers(SchedulerServer *server) {
  std::set<worker_id_t>::iterator iter = worker_ids_.begin();
  for (; iter != worker_ids_.end(); ++iter) {
    worker_id_t w_id = *iter;

    std::vector<job_id_t> outer_job_ids;
    std::vector<job_id_t> inner_job_ids;
    {
      std::map<worker_id_t, JobIdPtrList>::iterator it = worker_job_ids_.find(w_id);
      // Worker definitely has a job, copy or compute. -omidm
      assert(it != worker_job_ids_.end());
      JobIdPtrList::iterator i = it->second.begin();
      for (; i != it->second.end(); ++i) {
        inner_job_ids.push_back(*(*i));
      }
    }

    std::vector<physical_data_id_t> phy_ids;
    {
      std::map<worker_id_t, PhyIdPtrList>::iterator it = worker_phy_ids_.find(w_id);
      // a worker does not necessarily have physical data,
      // e.g. it has compute jobs with empty read/write set -omidm
      if (it != worker_phy_ids_.end()) {
        PhyIdPtrList::iterator i = it->second.begin();
        for (; i != it->second.end(); ++i) {
          phy_ids.push_back(*(*i));
        }
      }
    }

    StartCommandTemplateCommand cm(record_name_,
                                   inner_job_ids,
                                   outer_job_ids,
                                   phy_ids);

    SchedulerWorker *worker;
    if (!server->GetSchedulerWorkerById(worker, w_id)) {
      assert(false);
    }
    server->SendCommand(worker, &cm);
  }
}

void BindingTemplate::SendCommandTemplateFinalizeToWorkers(SchedulerServer *server) {
  std::set<worker_id_t>::iterator iter = worker_ids_.begin();
  for (; iter != worker_ids_.end(); ++iter) {
    worker_id_t w_id = *iter;

    EndCommandTemplateCommand cm(record_name_);

    SchedulerWorker *worker;
    if (!server->GetSchedulerWorkerById(worker, w_id)) {
      assert(false);
    }
    server->SendCommand(worker, &cm);
  }
}

void BindingTemplate::SpawnCommandTemplateAtWorkers(const std::vector<Parameter>& parameters,
                                                    const ExtraDependency& extra_dependency,
                                                    SchedulerServer *server,
                                                    const template_id_t& template_generation_id) {
  std::set<worker_id_t>::iterator iter = worker_ids_.begin();
  for (; iter != worker_ids_.end(); ++iter) {
    worker_id_t w_id = *iter;

    std::vector<job_id_t> outer_job_ids;
    std::vector<job_id_t> inner_job_ids;
    {
      std::map<worker_id_t, JobIdPtrList>::iterator it = worker_job_ids_.find(w_id);
      // Worker definitely has a job, copy or compute. -omidm
      assert(it != worker_job_ids_.end());
      JobIdPtrList::iterator i = it->second.begin();
      for (; i != it->second.end(); ++i) {
        inner_job_ids.push_back(*(*i));
      }
    }

    std::vector<Parameter> params;
    {
      std::map<worker_id_t, std::vector<size_t> >::iterator it =
        worker_parameter_indices_.find(w_id);
      // a worker does not necessarily have compute jobs, and so parameter -omidm
      if (it != worker_parameter_indices_.end()) {
        std::vector<size_t>::iterator i = it->second.begin();
        for (; i != it->second.end(); ++i) {
          assert((*i) <parameters.size());
          params.push_back(parameters[*i]);
        }
      }
    }

    std::vector<physical_data_id_t> phy_ids;
    {
      std::map<worker_id_t, PhyIdPtrList>::iterator it = worker_phy_ids_.find(w_id);
      // a worker does not necessarily have physical data,
      // e.g. it has compute jobs with empty read/write set -omidm
      if (it != worker_phy_ids_.end()) {
        PhyIdPtrList::iterator i = it->second.begin();
        for (; i != it->second.end(); ++i) {
          phy_ids.push_back(*(*i));
        }
      }
    }

    IDSet<job_id_t> ed;
    {
      ExtraDependency::const_iterator it =
        extra_dependency.find(w_id);
      if (it != extra_dependency.end()) {
        ed = it->second;
      }
    }

    SpawnCommandTemplateCommand cm(record_name_,
                                   inner_job_ids,
                                   outer_job_ids,
                                   ed,
                                   params,
                                   phy_ids,
                                   template_generation_id);

    SchedulerWorker *worker;
    if (!server->GetSchedulerWorkerById(worker, w_id)) {
      assert(false);
    }
    server->SendCommand(worker, &cm);
  }
}

void BindingTemplate::SendComputeJobCommand(ComputeJobCommandTemplate* command,
                                            const Parameter& parameter,
                                            const IDSet<job_id_t>& extra_dependency,
                                            SchedulerServer *server) {
  std::string job_name = command->job_name_;
  ID<job_id_t> job_id(*(command->job_id_ptr_));
  ID<job_id_t> future_job_id(*(command->future_job_id_ptr_));

  IDSet<physical_data_id_t> read_set, write_set, scratch_set, reduce_set;
  IDSet<job_id_t> before_set, after_set;

  {
    JobIdPtrSet::iterator it = command->before_set_ptr_.begin();
    for (; it != command->before_set_ptr_.end(); ++it) {
      before_set.insert(*(*it));
    }
  }
  {
    JobIdPtrSet::iterator it = command->after_set_ptr_.begin();
    for (; it != command->after_set_ptr_.end(); ++it) {
      after_set.insert(*(*it));
    }
  }
  {
    PhyIdPtrSet::iterator it = command->read_set_ptr_.begin();
    for (; it != command->read_set_ptr_.end(); ++it) {
      read_set.insert(*(*it));
    }
  }
  {
    PhyIdPtrSet::iterator it = command->write_set_ptr_.begin();
    for (; it != command->write_set_ptr_.end(); ++it) {
      write_set.insert(*(*it));
    }
  }
  {
    PhyIdPtrSet::iterator it = command->scratch_set_ptr_.begin();
    for (; it != command->scratch_set_ptr_.end(); ++it) {
      scratch_set.insert(*(*it));
    }
  }
  {
    PhyIdPtrSet::iterator it = command->reduce_set_ptr_.begin();
    for (; it != command->reduce_set_ptr_.end(); ++it) {
      reduce_set.insert(*(*it));
    }
  }

  ComputeJobCommand cm(job_name,
                       job_id,
                       read_set,
                       write_set,
                       scratch_set,
                       reduce_set,
                       before_set,
                       extra_dependency,
                       after_set,
                       future_job_id,
                       command->sterile_,
                       command->region_,
                       parameter);

  SchedulerWorker *worker;
  if (!server->GetSchedulerWorkerById(worker, command->worker_id_)) {
    assert(false);
  }
  server->SendCommand(worker, &cm);
}


void BindingTemplate::SendCombineJobCommand(CombineJobCommandTemplate* command,
                                            const IDSet<job_id_t>& extra_dependency,
                                            SchedulerServer *server) {
  std::string job_name = command->job_name_;
  ID<job_id_t> job_id(*(command->job_id_ptr_));

  IDSet<physical_data_id_t> scratch_set, reduce_set;
  IDSet<job_id_t> before_set;

  {
    JobIdPtrSet::iterator it = command->before_set_ptr_.begin();
    for (; it != command->before_set_ptr_.end(); ++it) {
      before_set.insert(*(*it));
    }
  }
  {
    PhyIdPtrSet::iterator it = command->scratch_set_ptr_.begin();
    for (; it != command->scratch_set_ptr_.end(); ++it) {
      scratch_set.insert(*(*it));
    }
  }
  {
    PhyIdPtrSet::iterator it = command->reduce_set_ptr_.begin();
    for (; it != command->reduce_set_ptr_.end(); ++it) {
      reduce_set.insert(*(*it));
    }
  }

  CombineJobCommand cm(job_name,
                       job_id,
                       scratch_set,
                       reduce_set,
                       before_set,
                       extra_dependency,
                       command->region_);

  SchedulerWorker *worker;
  if (!server->GetSchedulerWorkerById(worker, command->worker_id_)) {
    assert(false);
  }
  server->SendCommand(worker, &cm);
}


void BindingTemplate::SendLocalCopyCommand(LocalCopyCommandTemplate* command,
                                           const IDSet<job_id_t>& extra_dependency,
                                           SchedulerServer *server) {
  ID<job_id_t> job_id(*(command->job_id_ptr_));
  ID<physical_data_id_t> from_data_id(*(command->from_physical_data_id_ptr_));
  ID<physical_data_id_t> to_data_id(*(command->to_physical_data_id_ptr_));

  IDSet<job_id_t> before_set;

  {
    JobIdPtrSet::iterator it = command->before_set_ptr_.begin();
    for (; it != command->before_set_ptr_.end(); ++it) {
      before_set.insert(*(*it));
    }
  }

  LocalCopyCommand cm_c(job_id,
                        from_data_id,
                        to_data_id,
                        before_set,
                        extra_dependency);

  SchedulerWorker *worker;
  if (!server->GetSchedulerWorkerById(worker, command->worker_id_)) {
    assert(false);
  }
  server->SendCommand(worker, &cm_c);
}

void BindingTemplate::SendRemoteCopySendCommand(RemoteCopySendCommandTemplate* command,
                                                const IDSet<job_id_t>& extra_dependency,
                                                SchedulerServer *server) {
  ID<job_id_t> job_id(*(command->job_id_ptr_));
  ID<job_id_t> receive_job_id(*(command->receive_job_id_ptr_));
  ID<job_id_t> mega_rcr_job_id(*(command->mega_rcr_job_id_ptr_));
  ID<physical_data_id_t> from_data_id(*(command->from_physical_data_id_ptr_));
  ID<worker_id_t> to_worker_id(command->to_worker_id_);
  std::string to_ip = command->to_ip_;
  ID<port_t> to_port(command->to_port_);

  IDSet<job_id_t> before_set;

  {
    JobIdPtrSet::iterator it = command->before_set_ptr_.begin();
    for (; it != command->before_set_ptr_.end(); ++it) {
      before_set.insert(*(*it));
    }
  }

  RemoteCopySendCommand cm_s(job_id,
                             receive_job_id,
                             mega_rcr_job_id,
                             from_data_id,
                             to_worker_id,
                             to_ip,
                             to_port,
                             before_set,
                             extra_dependency);

  SchedulerWorker *worker;
  if (!server->GetSchedulerWorkerById(worker, command->worker_id_)) {
    assert(false);
  }
  server->SendCommand(worker, &cm_s);
}

void BindingTemplate::SendRemoteCopyReceiveCommand(RemoteCopyReceiveCommandTemplate* command,
                                                   const IDSet<job_id_t>& extra_dependency,
                                                   SchedulerServer *server) {
  ID<job_id_t> job_id(*(command->job_id_ptr_));
  ID<physical_data_id_t> to_data_id(*(command->to_physical_data_id_ptr_));

  IDSet<job_id_t> before_set;

  {
    JobIdPtrSet::iterator it = command->before_set_ptr_.begin();
    for (; it != command->before_set_ptr_.end(); ++it) {
      before_set.insert(*(*it));
    }
  }

  RemoteCopyReceiveCommand cm_r(job_id,
                                to_data_id,
                                before_set,
                                extra_dependency);

  SchedulerWorker *worker;
  if (!server->GetSchedulerWorkerById(worker, command->worker_id_)) {
    assert(false);
  }
  server->SendCommand(worker, &cm_r);
}

void BindingTemplate::SendMegaRCRCommand(MegaRCRCommandTemplate* command,
                                         const IDSet<job_id_t>& extra_dependency,
                                         SchedulerServer *server) {
  ID<job_id_t> job_id(*(command->job_id_ptr_));

  std::vector<job_id_t> receive_job_ids;
  {
    JobIdPtrList::iterator it = command->rcr_job_id_ptrs_.begin();
    for (; it != command->rcr_job_id_ptrs_.end(); ++it) {
      receive_job_ids.push_back(*(*it));
    }
  }

  std::vector<physical_data_id_t> to_physical_data_ids;
  {
    PhyIdPtrList::iterator it = command->to_physical_data_id_ptrs_.begin();
    for (; it != command->to_physical_data_id_ptrs_.end(); ++it) {
      to_physical_data_ids.push_back(*(*it));
    }
  }

  MegaRCRCommand cm_r(job_id,
                      receive_job_ids,
                      to_physical_data_ids,
                      extra_dependency);

  SchedulerWorker *worker;
  if (!server->GetSchedulerWorkerById(worker, command->worker_id_)) {
    assert(false);
  }
  server->SendCommand(worker, &cm_r);
}

bool BindingTemplate::TrackDataObject(const worker_id_t& worker_id,
                                      const logical_data_id_t& ldid,
                                      const physical_data_id_t& pdid,
                                      VersionType version_type,
                                      data_version_t version_diff_from_base) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  PhyIdPtrMap::iterator iter =  phy_id_map_.find(pdid);
  if (iter != phy_id_map_.end()) {
    return true;
  }

  if (version_type == SCRATCH) {
    dbg(DBG_ERROR, "ERROR: cannot track a scratch data as a new entry (pdid: %lu).\n", pdid);
    assert(false);
  }

  PhyIdPtr pdid_ptr = PhyIdPtr(new physical_data_id_t(pdid));
  phy_id_map_[pdid] = pdid_ptr;

  {
    PatternEntry *pattern =
      new PatternEntry(worker_id, ldid, pdid, version_type, version_diff_from_base);
    entry_pattern_map_[pdid] = pattern;

    PatternSorted::iterator iter = ordered_entry_patterns_.find(ldid);
    if (iter != ordered_entry_patterns_.end()) {
      if (version_type == REGULAR) {
        iter->second.first->push_back(pattern);
      } else {
        iter->second.second->push_back(pattern);
      }
    } else {
      PatternList *regular = new PatternList();
      PatternList *wild_card = new PatternList();
      if (version_type == REGULAR) {
        regular->push_back(pattern);
      } else {
        wild_card->push_back(pattern);
      }
      ordered_entry_patterns_[ldid] = std::make_pair(regular, wild_card);
    }
  }

  {
    PatternEntry *pattern =
      new PatternEntry(worker_id, ldid, pdid, version_type, version_diff_from_base);
    end_pattern_map_[pdid] = pattern;
  }

  return true;
}

bool BindingTemplate::UpdateDataObject(const physical_data_id_t& pdid,
                                       VersionType version_type,
                                       data_version_t version_diff_from_base) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  PatternMap::iterator iter =  end_pattern_map_.find(pdid);
  if (iter == end_pattern_map_.end()) {
    dbg(DBG_ERROR, "ERROR:, no record found to update the pdid: %lu.\n", pdid);
    assert(false);
    return false;
  }

  PatternEntry *pe = iter->second;
  pe->version_type_ = version_type;
  pe->version_diff_from_base_ = version_diff_from_base;

  return true;
}

bool BindingTemplate::AddComputeJobCommand(ComputeJobCommand* command,
                                           worker_id_t w_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  worker_ids_.insert(w_id);
  JobIdPtr compute_job_id_ptr = GetExistingComputeJobIdPtr(command->job_id().elem());

  if (worker_template_active_) {
    job_worker_map_[*compute_job_id_ptr].insert(w_id);
  }

  PhyIdPtrSet read_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->read_set_p()->begin();
    for (; iter != command->read_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      read_set.insert(phy_id_ptr);

      if (worker_template_active_) {
        phy_worker_map_[*phy_id_ptr].insert(w_id);
      }
    }
  }

  PhyIdPtrSet write_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->write_set_p()->begin();
    for (; iter != command->write_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      write_set.insert(phy_id_ptr);

      if (worker_template_active_) {
        phy_worker_map_[*phy_id_ptr].insert(w_id);
      }
    }
  }

  PhyIdPtrSet scratch_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->scratch_set_p()->begin();
    for (; iter != command->scratch_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      scratch_set.insert(phy_id_ptr);

      if (worker_template_active_) {
        phy_worker_map_[*phy_id_ptr].insert(w_id);
      }
    }
  }

  PhyIdPtrSet reduce_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->reduce_set_p()->begin();
    for (; iter != command->reduce_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      reduce_set.insert(phy_id_ptr);

      if (worker_template_active_) {
        phy_worker_map_[*phy_id_ptr].insert(w_id);
      }
    }
  }

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        if (!GetComputeJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      }
      before_set.insert(job_id_ptr);

      if (worker_template_active_) {
        job_worker_map_[*job_id_ptr].insert(w_id);
      }
    }
  }

  JobIdPtrSet after_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->after_set_p()->begin();
    for (; iter != command->after_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        if (!GetComputeJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      }
      after_set.insert(job_id_ptr);

      if (worker_template_active_) {
        job_worker_map_[*job_id_ptr].insert(w_id);
      }
    }
  }

  ComputeJobCommandTemplate *cm =
    new ComputeJobCommandTemplate(command->job_name(),
                                  compute_job_id_ptr,
                                  read_set,
                                  write_set,
                                  scratch_set,
                                  reduce_set,
                                  before_set,
                                  after_set,
                                  future_job_id_ptr_,
                                  command->sterile(),
                                  command->region(),
                                  w_id);

  // Keep this mapping to set the param_index in Finalize - omidm
  compute_job_to_command_map_[*compute_job_id_ptr] = cm;

  command_templates_.push_back(cm);

  return true;
}

bool BindingTemplate::AddCombineJobCommand(CombineJobCommand* command,
                                            worker_id_t w_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  worker_ids_.insert(w_id);
  JobIdPtr combine_job_id_ptr = GetCopyJobIdPtr(command->job_id().elem());

  if (worker_template_active_) {
    job_worker_map_[*combine_job_id_ptr].insert(w_id);
  }

  PhyIdPtrSet scratch_set;
  assert(command->scratch_set_p()->size() == 1);
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->scratch_set_p()->begin();
    for (; iter != command->scratch_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      scratch_set.insert(phy_id_ptr);

      if (worker_template_active_) {
        phy_worker_map_[*phy_id_ptr].insert(w_id);
      }
    }
  }

  PhyIdPtrSet reduce_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->reduce_set_p()->begin();
    for (; iter != command->reduce_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      reduce_set.insert(phy_id_ptr);

      if (worker_template_active_) {
        phy_worker_map_[*phy_id_ptr].insert(w_id);
      }
    }
  }

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        if (!GetComputeJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      }
      before_set.insert(job_id_ptr);

      if (worker_template_active_) {
        job_worker_map_[*job_id_ptr].insert(w_id);
      }
    }
  }

  CombineJobCommandTemplate *cm =
    new CombineJobCommandTemplate(command->job_name(),
                                  combine_job_id_ptr,
                                  scratch_set,
                                  reduce_set,
                                  before_set,
                                  command->region(),
                                  w_id);

  combine_job_to_command_map_[*combine_job_id_ptr] = cm;

  command_templates_.push_back(cm);

  return true;
}

bool BindingTemplate::AddLocalCopyCommand(LocalCopyCommand* command,
                                          worker_id_t w_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  worker_ids_.insert(w_id);
  JobIdPtr copy_job_id_ptr = GetCopyJobIdPtr(command->job_id().elem());

  if (worker_template_active_) {
    job_worker_map_[*copy_job_id_ptr].insert(w_id);
  }

  PhyIdPtr from_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->from_physical_data_id().elem());

  if (worker_template_active_) {
    phy_worker_map_[*from_physical_data_id_ptr].insert(w_id);
  }

  PhyIdPtr to_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->to_physical_data_id().elem());

  if (worker_template_active_) {
    phy_worker_map_[*to_physical_data_id_ptr].insert(w_id);
  }

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        if (!GetComputeJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      }
      before_set.insert(job_id_ptr);

      if (worker_template_active_) {
        job_worker_map_[*job_id_ptr].insert(w_id);
      }
    }
  }

  LocalCopyCommandTemplate *cm =
    new LocalCopyCommandTemplate(copy_job_id_ptr,
                                 from_physical_data_id_ptr,
                                 to_physical_data_id_ptr,
                                 before_set,
                                 w_id);

  lc_job_to_command_map_[*copy_job_id_ptr] = cm;

  command_templates_.push_back(cm);

  return true;
}

bool BindingTemplate::AddRemoteCopySendCommand(RemoteCopySendCommand* command,
                                               worker_id_t w_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  worker_ids_.insert(w_id);
  JobIdPtr copy_job_id_ptr = GetCopyJobIdPtr(command->job_id().elem());

  if (worker_template_active_) {
    job_worker_map_[*copy_job_id_ptr].insert(w_id);
  }

  JobIdPtr receive_job_id_ptr = GetCopyJobIdPtr(command->receive_job_id().elem());

  assert(command->mega_rcr_job_id().elem() == NIMBUS_KERNEL_JOB_ID);
  static const JobIdPtr default_mega_rcr_job_id_ptr =
    JobIdPtr(new job_id_t(NIMBUS_KERNEL_JOB_ID));

  if (worker_template_active_) {
    job_worker_map_[*receive_job_id_ptr].insert(w_id);
  }

  PhyIdPtr from_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->from_physical_data_id().elem());

  if (worker_template_active_) {
    phy_worker_map_[*from_physical_data_id_ptr].insert(w_id);
  }

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        if (!GetComputeJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      }
      before_set.insert(job_id_ptr);

      if (worker_template_active_) {
        job_worker_map_[*job_id_ptr].insert(w_id);
      }
    }
  }

  RemoteCopySendCommandTemplate *cm =
    new RemoteCopySendCommandTemplate(copy_job_id_ptr,
                                      receive_job_id_ptr,
                                      default_mega_rcr_job_id_ptr,
                                      from_physical_data_id_ptr,
                                      command->to_worker_id(),
                                      command->to_ip(),
                                      command->to_port(),
                                      before_set,
                                      w_id);

  rcs_job_to_command_map_[*copy_job_id_ptr] = cm;

  command_templates_.push_back(cm);

  return true;
}

bool BindingTemplate::AddRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* command,
                                                  worker_id_t w_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  worker_ids_.insert(w_id);
  JobIdPtr copy_job_id_ptr = GetCopyJobIdPtr(command->job_id().elem());

  if (worker_template_active_) {
    job_worker_map_[*copy_job_id_ptr].insert(w_id);
  }

  PhyIdPtr to_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->to_physical_data_id().elem());

  if (worker_template_active_) {
    phy_worker_map_[*to_physical_data_id_ptr].insert(w_id);
  }

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        if (!GetComputeJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      }
      before_set.insert(job_id_ptr);

      if (worker_template_active_) {
        job_worker_map_[*job_id_ptr].insert(w_id);
      }
    }
  }

  if (before_set.size() == 0) {
    candidate_rcr_.insert(*copy_job_id_ptr);
  }

  RemoteCopyReceiveCommandTemplate *cm =
    new RemoteCopyReceiveCommandTemplate(copy_job_id_ptr,
                                         to_physical_data_id_ptr,
                                         before_set,
                                         w_id);

  // Keep this mapping for mega rcr job computations.
  rcr_job_to_command_map_[*copy_job_id_ptr] = cm;

  command_templates_.push_back(cm);

  return true;
}

BindingTemplate::JobIdPtr BindingTemplate::GetCopyJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = copy_job_id_map_.find(job_id);
  if (iter != copy_job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    job_id_ptr = JobIdPtr(new job_id_t(job_id));
    copy_job_id_map_[job_id] = job_id_ptr;
  }

  return job_id_ptr;
}

bool BindingTemplate::GetCopyJobIdPtrIfExisted(job_id_t job_id, JobIdPtr *job_id_ptr) {
  JobIdPtrMap::iterator iter = copy_job_id_map_.find(job_id);
  if (iter != copy_job_id_map_.end()) {
    *job_id_ptr = iter->second;
    return true;
  }

  return false;
}

BindingTemplate::JobIdPtr BindingTemplate::GetExistingCopyJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = copy_job_id_map_.find(job_id);
  if (iter != compute_job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    assert(false);
  }

  return job_id_ptr;
}

BindingTemplate::JobIdPtr BindingTemplate::GetExistingComputeJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = compute_job_id_map_.find(job_id);
  if (iter != compute_job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    assert(false);
  }

  return job_id_ptr;
}

bool BindingTemplate::GetComputeJobIdPtrIfExisted(job_id_t job_id, JobIdPtr *job_id_ptr) {
  JobIdPtrMap::iterator iter = compute_job_id_map_.find(job_id);
  if (iter != compute_job_id_map_.end()) {
    *job_id_ptr = iter->second;
    return true;
  }

  return false;
}


BindingTemplate::PhyIdPtr BindingTemplate::GetExistingPhyIdPtr(physical_data_id_t pdid) {
  PhyIdPtr phy_id_ptr;

  PhyIdPtrMap::iterator iter = phy_id_map_.find(pdid);
  if (iter != phy_id_map_.end()) {
    phy_id_ptr = iter->second;
  } else {
    assert(false);
  }

  return phy_id_ptr;
}


