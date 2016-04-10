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
  * This is ExecutionTemplate module to hold and instantiate execution template
  * for a worker without building/destroying the execution graph each time.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/execution_template.h"
#include "src/shared/scheduler_client.h"

using namespace nimbus; // NOLINT

ExecutionTemplate::ExecutionTemplate(const std::string& execution_template_name,
                                     const std::vector<job_id_t>& inner_job_ids,
                                     const std::vector<job_id_t>& outer_job_ids,
                                     const std::vector<physical_data_id_t>& phy_ids) {
  finalized_ = false;
  mark_stat_ = false;
  copy_job_num_ = 0;
  compute_job_num_ = 0;
  job_done_counter_ = 0;
  ready_job_counter_ = 0;
  pending_instantiate_ = false;
  execution_template_name_ = execution_template_name;
  template_generation_id_ = NIMBUS_INIT_TEMPLATE_ID;
  future_job_id_ptr_ = JobIdPtr(new job_id_t(0));

  {
    std::vector<job_id_t>::const_iterator iter = inner_job_ids.begin();
    for (; iter != inner_job_ids.end(); ++iter) {
      JobIdPtr job_id_ptr = JobIdPtr(new job_id_t(*iter));
      inner_job_id_map_[*iter] = job_id_ptr;
      inner_job_id_list_.push_back(job_id_ptr);
    }
  }

  {
    std::vector<job_id_t>::const_iterator iter = outer_job_ids.begin();
    for (; iter != outer_job_ids.end(); ++iter) {
      JobIdPtr job_id_ptr = JobIdPtr(new job_id_t(*iter));
      outer_job_id_map_[*iter] = job_id_ptr;
      outer_job_id_list_.push_back(job_id_ptr);
    }
  }

  {
    std::vector<physical_data_id_t>::const_iterator iter = phy_ids.begin();
    for (; iter != phy_ids.end(); ++iter) {
      PhyIdPtr phy_id_ptr = PhyIdPtr(new physical_data_id_t(*iter));
      phy_id_map_[*iter] = phy_id_ptr;
      phy_id_list_.push_back(phy_id_ptr);
    }
  }
}

ExecutionTemplate::~ExecutionTemplate() {
}

bool ExecutionTemplate::finalized() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return finalized_;
}

size_t ExecutionTemplate::job_num() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(finalized_);
  return copy_job_num_ + compute_job_num_;
}

size_t ExecutionTemplate::copy_job_num() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(finalized_);
  return copy_job_num_;
}

size_t ExecutionTemplate::compute_job_num() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(finalized_);
  return compute_job_num_;
}

std::string ExecutionTemplate::execution_template_name() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return execution_template_name_;
}

template_id_t ExecutionTemplate::template_generation_id() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return template_generation_id_;
}

bool ExecutionTemplate::pending_instantiate() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return pending_instantiate_;
}

template_id_t ExecutionTemplate::pending_template_generation_id() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return pending_template_generation_id_;
}

size_t ExecutionTemplate::ready_job_counter() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return ready_job_counter_;
}

bool ExecutionTemplate::Finalize() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(!finalized_);
  assert(seed_job_templates_.size() == 0);
  assert(job_templates_list_.size() == 0);

  JobTemplateMap::iterator iter = job_templates_.begin();
  for (; iter != job_templates_.end(); ++iter) {
    JobTemplate *jt = iter->second;
    if (jt->type_ == COMPUTE) {
      compute_job_id_list_.push_back(jt->job_id_ptr_);
    }

    job_templates_list_.push_back(jt);
    if (jt->dependency_num_ == 0) {
      seed_job_templates_.push_back(jt);
      continue;
    }

    IDSet<job_id_t>::IDSetIter it = jt->before_set_.begin();
    for (; it != jt->before_set_.end(); ++it) {
      JobTemplateMap::iterator i = job_templates_.find(*it);
      assert(i != job_templates_.end());
      i->second->after_set_job_templates_.push_back(jt);
    }
  }

  finalized_ = true;

  return true;
}

bool ExecutionTemplate::Instantiate(const std::vector<job_id_t>& inner_job_ids,
                                    const std::vector<job_id_t>& outer_job_ids,
                                    const IDSet<job_id_t>& extra_dependency,
                                    const std::vector<Parameter>& parameters,
                                    const std::vector<physical_data_id_t>& physical_ids,
                                    const WorkerDataExchanger::EventList& pending_events,
                                    const template_id_t& template_generation_id,
                                    JobList *ready_jobs) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (job_done_counter_ != 0) {
    std::cout << "WARNING: Instantiating a new execution template "
              << "before finishing the previous one!\n"
              << "         Number of left jobs: " << job_done_counter_ << std::endl;

    assert(!pending_instantiate_);

    pending_instantiate_ = true;
    pending_inner_job_ids_ = inner_job_ids;
    pending_outer_job_ids_ = outer_job_ids;
    pending_extra_dependency_ = extra_dependency;
    pending_parameters_ = parameters;
    pending_physical_ids_ = physical_ids;
    pending_template_generation_id_ = template_generation_id;
    ready_jobs->clear();
    return false;
  }

  assert(ready_job_counter_ == 0);
  assert(!pending_instantiate_);
  assert(finalized_);
  assert(inner_job_ids.size() == inner_job_id_list_.size());
  assert(outer_job_ids.size() == outer_job_id_list_.size());
  assert(physical_ids.size() == phy_id_list_.size());
  assert(blocked_on_extra_dependency_.size() == 0);
  // assert(seed_job_templates_.size() > 0);

  mark_stat_ = false;

  extra_dependency_ = extra_dependency;
  if (extra_dependency_.size() != 0) {
    dbg(DBG_WARN, "WARNING: extra dependency size of %lu for %s!\n",
        extra_dependency_.size(), execution_template_name_.c_str());
  }

  template_generation_id_ = template_generation_id;
  job_done_counter_ = copy_job_num_ + compute_job_num_;

  {
    size_t idx = 0;
    JobIdPtrList::iterator iter = inner_job_id_list_.begin();
    for (; iter != inner_job_id_list_.end(); ++iter) {
      *(*iter) = inner_job_ids[idx];
      ++idx;
    }
  }

  {
    size_t idx = 0;
    JobIdPtrList::iterator iter = outer_job_id_list_.begin();
    for (; iter != outer_job_id_list_.end(); ++iter) {
      *(*iter) = outer_job_ids[idx];
      ++idx;
    }
  }

  {
    size_t idx = 0;
    PhyIdPtrList::iterator iter = phy_id_list_.begin();
    for (; iter != phy_id_list_.end(); ++iter) {
      *(*iter) = physical_ids[idx];
      ++idx;
    }
  }

  parameters_ = parameters;

  {
    JobTemplateVector::iterator iter = job_templates_list_.begin();
    for (; iter != job_templates_list_.end(); ++iter) {
      JobTemplate *jt = *iter;
      assert(jt->dependency_counter_ == 0);
      jt->dependency_counter_ = jt->dependency_num_;
    }
  }

  ready_jobs->clear();
  size_t counter = 0;
  {
    JobTemplateVector::iterator iter = seed_job_templates_.begin();
    for (; iter != seed_job_templates_.end(); ++iter) {
      JobTemplate *jt = *iter;
      jt->Refresh(parameters_, template_generation_id_);
      if (extra_dependency_.size() == 0) {
        ready_jobs->push_back((jt->job_));
        ++counter;
      } else {
        blocked_on_extra_dependency_.push_back((jt->job_));
      }
    }
  }

  ready_job_counter_ += counter;

  {
    WorkerDataExchanger::EventList::const_iterator iter = pending_events.begin();
    for (; iter != pending_events.end(); ++iter) {
      ProcessReceiveEvent(*iter, ready_jobs, true);
    }
  }


  return true;
}

bool ExecutionTemplate::InstantiatePending(const WorkerDataExchanger::EventList& pending_events,
                                           JobList *ready_jobs) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(pending_instantiate_);
  assert(job_done_counter_ == 0);
  assert(ready_job_counter_ == 0);
  assert(extra_dependency_.size() == 0);
  assert(blocked_on_extra_dependency_.size() == 0);
  pending_instantiate_ = false;

  Instantiate(pending_inner_job_ids_,
              pending_outer_job_ids_,
              pending_extra_dependency_,
              pending_parameters_,
              pending_physical_ids_,
              pending_events,
              pending_template_generation_id_,
              ready_jobs);

  return true;
}

bool ExecutionTemplate::MarkInnerJobDone(const job_id_t& shadow_job_id,
                                         JobList *ready_jobs,
                                         bool prepare_rewind_phase,
                                         bool mark_stat,
                                         bool append) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(extra_dependency_.size() == 0);
  assert(finalized_);
  if (!append) {
    ready_jobs->clear();
  }

  bool et_complete = false;

  if (mark_stat) {
    mark_stat_ = mark_stat;
  }

  assert(job_done_counter_ > 0);
  assert(ready_job_counter_ > 0);
  --job_done_counter_;
  --ready_job_counter_;
  if (job_done_counter_ == 0) {
    assert(ready_job_counter_ == 0);
    et_complete = true;
  }

  if (prepare_rewind_phase) {
    return et_complete;
  }

  JobTemplateMap::iterator it = job_templates_.find(shadow_job_id);
  assert(it != job_templates_.end());
  JobTemplateVector ready_list;
  it->second->ClearAfterSet(&ready_list);

  size_t counter = 0;
  JobTemplateVector::iterator iter = ready_list.begin();
  for (; iter != ready_list.end(); ++iter) {
    JobTemplate *jt = *iter;
    jt->Refresh(parameters_, template_generation_id_);
    ready_jobs->push_back((jt->job_));
    ++counter;
  }

  ready_job_counter_ += counter;

  return et_complete;
}


void ExecutionTemplate::NotifyJobDone(const job_id_t& job_id,
                                      JobList *ready_jobs,
                                      bool &prepare_rewind_phase,
                                      bool append) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(finalized_);
  if (!append) {
    ready_jobs->clear();
  }

  if (prepare_rewind_phase) {
    return;
  }

  pending_extra_dependency_.remove(job_id);

  if (extra_dependency_.size() == 0) {
    return;
  }

  size_t counter = 0;
  extra_dependency_.remove(job_id);
  if (extra_dependency_.size() == 0) {
    JobList::iterator iter = blocked_on_extra_dependency_.begin();
    for (; iter != blocked_on_extra_dependency_.end(); ++iter) {
      ready_jobs->push_back(*iter);
      ++counter;
    }
    blocked_on_extra_dependency_.clear();
  }

  ready_job_counter_ += counter;

  return;
}


bool ExecutionTemplate::GenerateMegaJobDoneCommand(MegaJobDoneCommand **cm) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(finalized_);
  std::vector<job_id_t> job_ids;

  JobIdPtrList::iterator iter = compute_job_id_list_.begin();
  for (; iter != compute_job_id_list_.end(); ++iter) {
    job_ids.push_back(*(*iter));
  }

  *cm = new MegaJobDoneCommand(job_ids, mark_stat_);
  return true;
}



void ExecutionTemplate::ProcessReceiveEvent(const WorkerDataExchanger::Event& e,
                                            JobList *ready_jobs,
                                            bool append) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(finalized_);
  if (!append) {
    ready_jobs->clear();
  }

  assert(e.template_generation_id_ == template_generation_id_);

  size_t counter = 0;
  JobTemplateMap::iterator iter;
  if (e.mega_rcr_job_id_ == NIMBUS_KERNEL_JOB_ID) {
    iter = job_templates_.find(e.receive_job_id_);
    assert(iter != job_templates_.end());
    JobTemplate *jt = iter->second;
    RemoteCopyReceiveJob *rcr = dynamic_cast<RemoteCopyReceiveJob*>(jt->job_); // NOLINT
    assert(rcr != NULL);
    rcr->set_data_version(e.version_);
    rcr->set_serialized_data(e.ser_data_);

    assert(jt->dependency_counter_ > 0);
    --(jt->dependency_counter_);
    if (jt->dependency_counter_ == 0) {
      jt->Refresh(parameters_, template_generation_id_);
      if (extra_dependency_.size() == 0) {
        ready_jobs->push_back(rcr);
        ++counter;
      } else {
        blocked_on_extra_dependency_.push_back(rcr);
      }
    }
  } else {
    iter = job_templates_.find(e.mega_rcr_job_id_);
    assert(iter != job_templates_.end());
    JobTemplate *jt = iter->second;
    MegaRCRJob *mrcr = dynamic_cast<MegaRCRJob*>(jt->job_); // NOLINT
    assert(mrcr != NULL);
    // TODO(omidm): fix the version setting for mega rcr job.
    mrcr->set_serialized_data(e.receive_job_id_, e.ser_data_);

    assert(jt->dependency_counter_ > 0);
    --(jt->dependency_counter_);
    if (jt->dependency_counter_ == 0) {
      jt->Refresh(parameters_, template_generation_id_);
      if (extra_dependency_.size() == 0) {
        ready_jobs->push_back(mrcr);
        ++counter;
      } else {
        blocked_on_extra_dependency_.push_back(mrcr);
      }
    }
  }

  ready_job_counter_ += counter;
}

bool ExecutionTemplate::AddComputeJobTemplate(ComputeJobCommand* cm,
                                              Application *app) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(!finalized_);

  job_id_t job_id = cm->job_id().elem();
  Job* job = app->CloneJob(cm->job_name());
  job->set_name("Compute:" + cm->job_name());
  job->set_sterile(cm->sterile());
  job->set_region(cm->region());
  job->set_shadow_job_id(job_id);
  job->set_execution_template(this);

  JobIdPtr compute_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  PhyIdPtrSet read_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = cm->read_set_p()->begin();
    for (; iter != cm->read_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      read_set.insert(phy_id_ptr);
    }
  }

  PhyIdPtrSet write_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = cm->write_set_p()->begin();
    for (; iter != cm->write_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      write_set.insert(phy_id_ptr);
    }
  }

  PhyIdPtrSet scratch_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = cm->scratch_set_p()->begin();
    for (; iter != cm->scratch_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      scratch_set.insert(phy_id_ptr);
    }
  }

  PhyIdPtrSet reduce_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = cm->reduce_set_p()->begin();
    for (; iter != cm->reduce_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      reduce_set.insert(phy_id_ptr);
    }
  }

  ComputeJobTemplate *jt =
    new ComputeJobTemplate(job,
                           compute_job_id_ptr,
                           read_set,
                           write_set,
                           scratch_set,
                           reduce_set,
                           cm->before_set(),
                           future_job_id_ptr_,
                           compute_job_num_++);

  job_templates_[job_id] = jt;

  return true;
}

void ExecutionTemplate::ComputeJobTemplate::Refresh(
    const std::vector<Parameter>& parameters,
    const template_id_t& template_generation_id) {
  job_->set_id(ID<job_id_t>(*job_id_ptr_));
  job_->set_future_job_id(ID<job_id_t>(*future_job_id_ptr_));

  job_->clear_template_variables();

  IDSet<physical_data_id_t> read_set, write_set, scratch_set, reduce_set;

  {
    PhyIdPtrSet::iterator it = read_set_ptr_.begin();
    for (; it != read_set_ptr_.end(); ++it) {
      read_set.insert(*(*it));
    }
  }
  job_->set_read_set(read_set);

  {
    PhyIdPtrSet::iterator it = write_set_ptr_.begin();
    for (; it != write_set_ptr_.end(); ++it) {
      write_set.insert(*(*it));
    }
  }
  job_->set_write_set(write_set);

  {
    PhyIdPtrSet::iterator it = scratch_set_ptr_.begin();
    for (; it != scratch_set_ptr_.end(); ++it) {
      scratch_set.insert(*(*it));
    }
  }
  job_->set_scratch_set(scratch_set);

  {
    PhyIdPtrSet::iterator it = reduce_set_ptr_.begin();
    for (; it != reduce_set_ptr_.end(); ++it) {
      reduce_set.insert(*(*it));
    }
  }
  job_->set_reduce_set(reduce_set);

  job_->set_parameters(parameters[param_index_]);
}

bool ExecutionTemplate::AddCombineJobTemplate(CombineJobCommand* cm,
                                              Application *app) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(!finalized_);

  job_id_t job_id = cm->job_id().elem();
  Job* job = app->CloneJob(cm->job_name());
  job->set_name("Combine:" + cm->job_name());
  job->set_region(cm->region());
  job->set_shadow_job_id(job_id);
  job->set_execution_template(this);

  JobIdPtr combine_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  PhyIdPtrSet scratch_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = cm->scratch_set_p()->begin();
    for (; iter != cm->scratch_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      scratch_set.insert(phy_id_ptr);
    }
  }

  PhyIdPtrSet reduce_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = cm->reduce_set_p()->begin();
    for (; iter != cm->reduce_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      reduce_set.insert(phy_id_ptr);
    }
  }

  CombineJobTemplate *jt =
    new CombineJobTemplate(job,
                           combine_job_id_ptr,
                           scratch_set,
                           reduce_set,
                           cm->before_set());

  ++copy_job_num_;
  job_templates_[job_id] = jt;

  return true;
}

void ExecutionTemplate::CombineJobTemplate::Refresh(
    const std::vector<Parameter>& parameters,
    const template_id_t& template_generation_id) {
  job_->set_id(ID<job_id_t>(*job_id_ptr_));

  job_->clear_template_variables();

  IDSet<physical_data_id_t> scratch_set, reduce_set;

  {
    PhyIdPtrSet::iterator it = scratch_set_ptr_.begin();
    for (; it != scratch_set_ptr_.end(); ++it) {
      scratch_set.insert(*(*it));
    }
  }
  job_->set_scratch_set(scratch_set);

  {
    PhyIdPtrSet::iterator it = reduce_set_ptr_.begin();
    for (; it != reduce_set_ptr_.end(); ++it) {
      reduce_set.insert(*(*it));
    }
  }
  job_->set_reduce_set(reduce_set);
}

bool ExecutionTemplate::AddLocalCopyJobTemplate(LocalCopyCommand* cm,
                                                Application *app) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(!finalized_);

  job_id_t job_id = cm->job_id().elem();
  LocalCopyJob * job = new LocalCopyJob(app);
  job->set_name("LocalCopy");
  job->set_shadow_job_id(job_id);
  job->set_execution_template(this);

  JobIdPtr copy_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  PhyIdPtr from_physical_data_id_ptr =
    GetExistingPhyIdPtr(cm->from_physical_data_id().elem());

  PhyIdPtr to_physical_data_id_ptr =
    GetExistingPhyIdPtr(cm->to_physical_data_id().elem());

  LocalCopyJobTemplate *jt =
    new LocalCopyJobTemplate(job,
                             copy_job_id_ptr,
                             from_physical_data_id_ptr,
                             to_physical_data_id_ptr,
                             cm->before_set());

  ++copy_job_num_;
  job_templates_[job_id] = jt;

  return true;
}

void ExecutionTemplate::LocalCopyJobTemplate::Refresh(
    const std::vector<Parameter>& parameters,
    const template_id_t& template_generation_id) {
  job_->set_id(ID<job_id_t>(*job_id_ptr_));

  IDSet<physical_data_id_t> read_set;
  read_set.insert(*from_physical_data_id_ptr_);
  job_->set_read_set(read_set);

  IDSet<physical_data_id_t> write_set;
  write_set.insert(*to_physical_data_id_ptr_);
  job_->set_write_set(write_set);
}

bool ExecutionTemplate::AddRemoteCopySendJobTemplate(RemoteCopySendCommand* cm,
                                                     Application *app,
                                                     WorkerDataExchanger *dx) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(!finalized_);

  job_id_t job_id = cm->job_id().elem();
  RemoteCopySendJob * job = new RemoteCopySendJob(dx, app);
  dx->AddContactInfo(cm->to_worker_id().elem(),
                     cm->to_ip(),
                     cm->to_port().elem());
  job->set_name("RemoteCopySend");
  job->set_to_worker_id(cm->to_worker_id());
  job->set_to_ip(cm->to_ip());
  job->set_to_port(cm->to_port());
  job->set_receive_job_id(cm->receive_job_id());
  job->set_mega_rcr_job_id(cm->mega_rcr_job_id());
  job->set_shadow_job_id(job_id);
  job->set_execution_template(this);

  JobIdPtr copy_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  // JobIdPtr receive_job_id_ptr = GetExistingInnerJobIdPtr(cm->receive_job_id().elem());

  // static const JobIdPtr default_mega_rcr_job_id_ptr =
  //   JobIdPtr(new job_id_t(NIMBUS_KERNEL_JOB_ID));

  // JobIdPtr mega_rcr_job_id_ptr;
  // if (cm->mega_rcr_job_id().elem() == NIMBUS_KERNEL_JOB_ID) {
  //   mega_rcr_job_id_ptr = default_mega_rcr_job_id_ptr;
  // } else {
  //   mega_rcr_job_id_ptr = GetExistingInnerJobIdPtr(cm->mega_rcr_job_id().elem());
  // }

  PhyIdPtr from_physical_data_id_ptr =
    GetExistingPhyIdPtr(cm->from_physical_data_id().elem());

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = cm->before_set_p()->begin();
    for (; iter != cm->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr = GetExistingInnerJobIdPtr(*iter);
      before_set.insert(job_id_ptr);
    }
  }

  RemoteCopySendJobTemplate *jt =
    new RemoteCopySendJobTemplate(job,
                                  copy_job_id_ptr,
                                  // receive_job_id_ptr,
                                  // mega_rcr_job_id_ptr,
                                  from_physical_data_id_ptr,
                                  cm->before_set());

  ++copy_job_num_;
  job_templates_[job_id] = jt;

  return true;
}

void ExecutionTemplate::RemoteCopySendJobTemplate::Refresh(
    const std::vector<Parameter>& parameters,
    const template_id_t& template_generation_id) {
  RemoteCopySendJob *rcs = reinterpret_cast<RemoteCopySendJob*>(job_);
  assert(rcs);
  rcs->set_id(ID<job_id_t>(*job_id_ptr_));
  // rcs->set_receive_job_id(ID<job_id_t>(*receive_job_id_ptr_));
  // rcs->set_mega_rcr_job_id(ID<job_id_t>(*mega_rcr_job_id_ptr_));
  rcs->set_template_generation_id(template_generation_id);

  IDSet<physical_data_id_t> read_set;
  read_set.insert(*from_physical_data_id_ptr_);
  rcs->set_read_set(read_set);
}

bool ExecutionTemplate::AddRemoteCopyReceiveJobTemplate(RemoteCopyReceiveCommand* cm,
                                                        Application *app) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(!finalized_);

  job_id_t job_id = cm->job_id().elem();
  RemoteCopyReceiveJob * job = new RemoteCopyReceiveJob(app);
  job->set_name("RemoteCopyReceive");
  job->set_id(cm->job_id());
  job->set_shadow_job_id(job_id);
  job->set_execution_template(this);

  JobIdPtr copy_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  PhyIdPtr to_physical_data_id_ptr =
    GetExistingPhyIdPtr(cm->to_physical_data_id().elem());

  RemoteCopyReceiveJobTemplate *jt =
    new RemoteCopyReceiveJobTemplate(job,
                                     copy_job_id_ptr,
                                     to_physical_data_id_ptr,
                                     cm->before_set());

  ++copy_job_num_;
  job_templates_[job_id] = jt;

  return true;
}

void ExecutionTemplate::RemoteCopyReceiveJobTemplate::Refresh(
    const std::vector<Parameter>& parameters,
    const template_id_t& template_generation_id) {
  job_->set_id(ID<job_id_t>(*job_id_ptr_));

  IDSet<physical_data_id_t> write_set;
  write_set.insert(*to_physical_data_id_ptr_);
  job_->set_write_set(write_set);
}

bool ExecutionTemplate::AddMegaRCRJobTemplate(MegaRCRCommand* cm,
                                              Application *app) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(!finalized_);

  job_id_t job_id = cm->job_id().elem();
  MegaRCRJob *job = new MegaRCRJob(app);
  job->set_name("MegaRCR");
  job->set_receive_job_ids(cm->receive_job_ids());
  job->set_shadow_job_id(job_id);
  job->set_execution_template(this);

  JobIdPtr copy_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  // JobIdPtrList receive_job_id_ptrs;
  // {
  //   std::vector<job_id_t>::const_iterator iter = cm->receive_job_ids_p()->begin();
  //   for (; iter != cm->receive_job_ids_p()->end(); ++iter) {
  //     JobIdPtr job_id_ptr = GetExistingInnerJobIdPtr(*iter);
  //     receive_job_id_ptrs.push_back(job_id_ptr);
  //   }
  // }

  PhyIdPtrList to_phy_id_ptrs;
  {
    std::vector<physical_data_id_t>::const_iterator iter =
      cm->to_physical_data_ids_p()->begin();
    for (; iter != cm->to_physical_data_ids_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      to_phy_id_ptrs.push_back(phy_id_ptr);
    }
  }

  IDSet<job_id_t> empty_before_set;

  MegaRCRJobTemplate *jt =
    new MegaRCRJobTemplate(job,
                           copy_job_id_ptr,
                           // receive_job_id_ptrs,
                           to_phy_id_ptrs,
                           empty_before_set);

  ++copy_job_num_;
  job_templates_[job_id] = jt;

  return true;
}

void ExecutionTemplate::MegaRCRJobTemplate::Refresh(
    const std::vector<Parameter>& parameters,
    const template_id_t& template_generation_id) {
  MegaRCRJob *mrcr = reinterpret_cast<MegaRCRJob*>(job_); // NOLINT
  assert(mrcr);

  mrcr->set_id(ID<job_id_t>(*job_id_ptr_));

  // std::vector<job_id_t> receive_job_ids;
  // {
  //   JobIdPtrList::iterator it = receive_job_id_ptrs_.begin();
  //   for (; it != receive_job_id_ptrs_.end(); ++it) {
  //     receive_job_ids.push_back(*(*it));
  //   }
  // }
  // mrcr->set_receive_job_ids(receive_job_ids);

  std::vector<physical_data_id_t> to_phy_ids;
  {
    PhyIdPtrList::iterator it = to_phy_id_ptrs_.begin();
    for (; it != to_phy_id_ptrs_.end(); ++it) {
      to_phy_ids.push_back(*(*it));
    }
  }
  mrcr->set_to_phy_ids(to_phy_ids);
}

ExecutionTemplate::JobIdPtr ExecutionTemplate::GetExistingInnerJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = inner_job_id_map_.find(job_id);
  if (iter != inner_job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    assert(false);
  }

  return job_id_ptr;
}

ExecutionTemplate::PhyIdPtr ExecutionTemplate::GetExistingPhyIdPtr(physical_data_id_t pdid) {
  PhyIdPtr phy_id_ptr;

  PhyIdPtrMap::iterator iter = phy_id_map_.find(pdid);
  if (iter != phy_id_map_.end()) {
    phy_id_ptr = iter->second;
  } else {
    assert(false);
  }

  return phy_id_ptr;
}

void ExecutionTemplate::JobTemplate::ClearAfterSet(JobTemplateVector *ready_list) {
  ready_list->clear();

  JobTemplateVector::iterator iter = after_set_job_templates_.begin();
  for (; iter != after_set_job_templates_.end(); ++iter) {
    assert((*iter)->dependency_counter_ > 0);
    --((*iter)->dependency_counter_);
    if ((*iter)->dependency_counter_ == 0) {
      ready_list->push_back(*iter);
    }
  }
}


