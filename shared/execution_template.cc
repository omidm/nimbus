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

#include "shared/execution_template.h"
#include "shared/scheduler_client.h"

using namespace nimbus; // NOLINT

ExecutionTemplate::ExecutionTemplate(const std::string& execution_template_name,
                                     const std::vector<job_id_t>& inner_job_ids,
                                     const std::vector<job_id_t>& outer_job_ids,
                                     const std::vector<physical_data_id_t>& phy_ids) {
  finalized_ = false;
  compute_job_num_ = 0;
  copy_job_num_ = 0;
  execution_template_name_ = execution_template_name;
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
  boost::unique_lock<boost::mutex> lock(mutex_);
  return finalized_;
}

size_t ExecutionTemplate::copy_job_num() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  return copy_job_num_;
}

size_t ExecutionTemplate::compute_job_num() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  return compute_job_num_;
}

std::string ExecutionTemplate::execution_template_name() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  return execution_template_name_;
}

bool ExecutionTemplate::Finalize() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(!finalized_);
  finalized_ = true;
  return true;
}

bool ExecutionTemplate::Instantiate(const std::vector<job_id_t>& inner_job_ids,
                                  const std::vector<job_id_t>& outer_job_ids,
                                  const std::vector<Parameter>& parameters,
                                  const std::vector<physical_data_id_t>& physical_ids,
                                  JobList *seed_jobs) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(finalized_);
  assert(inner_job_ids.size() == inner_job_id_list_.size());
  assert(outer_job_ids.size() == outer_job_id_list_.size());
  assert(physical_ids.size() == phy_id_list_.size());

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

  ComputeJobTemplate *cc;
  JobTemplateMap::iterator iter = job_templates_.begin();
  for (; iter != job_templates_.end(); ++iter) {
    BaseJobTemplate *ct = iter->second;
    switch (ct->type_) {
      case COMPUTE:
        cc = reinterpret_cast<ComputeJobTemplate*>(ct);
        assert(cc->param_index_ < parameters.size());
        break;
      case LC:
        break;
      case RCS:
        break;
      case RCR:
        break;
      case MEGA_RCR:
        break;
      default:
        assert(false);
    }
  }

  return true;
}

bool ExecutionTemplate::MarkJobDone(const job_id_t& shadow_job_id,
                                    JobList *ready_jobs) {
  return false;
}

bool ExecutionTemplate::AddComputeJobTemplate(ComputeJobCommand* cm,
                                              Application *app) {
  std::string job_name = cm->job_name();
  job_id_t job_id = cm->job_id().elem();
  Job* job = app->CloneJob(cm->job_name());
  job->set_name("Compute:" + cm->job_name());
  job->set_sterile(cm->sterile());
  job->set_region(cm->region());
  // job->set_shadow_id(job_id);
  // job->set_execution_template(this);

  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(!finalized_);

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

  ComputeJobTemplate *jt =
    new ComputeJobTemplate(job,
                           compute_job_id_ptr,
                           read_set,
                           write_set,
                           cm->before_set(),
                           future_job_id_ptr_,
                           compute_job_num_++);

  job_templates_[job_id] = jt;

  return true;
}


void ExecutionTemplate::BaseJobTemplate::RemoveDependences(JobTemplateVector *ready_list) {
  ready_list->clear();

  JobTemplateVector::iterator iter = after_set_job_templates_.begin();
  for (; iter != after_set_job_templates_.end(); ++iter) {
    assert((*iter)->dependency_counter_ > 0);
    --(*iter)->dependency_counter_;
    if ((*iter)->dependency_counter_ == 0) {
      ready_list->push_back(*iter);
    }
  }
}

void ExecutionTemplate::ComputeJobTemplate::Refresh(
    const std::vector<Parameter>& parameters) {
  job_->set_id(ID<job_id_t>(*job_id_ptr_));
  job_->set_future_job_id(ID<job_id_t>(*future_job_id_ptr_));

  IDSet<physical_data_id_t> read_set, write_set;

  assert(dependency_counter_ == 0);
  dependency_counter_ = before_set_count_;

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

  job_->set_parameters(parameters[param_index_]);
}




bool ExecutionTemplate::AddLocalCopyJobTemplate(LocalCopyCommand* cm,
                                                Application *app) {
  job_id_t job_id = cm->job_id().elem();
  LocalCopyJob * job = new LocalCopyJob(app);
  job->set_name("LocalCopy");
  // job->set_shadow_id(job_id);
  // job->set_execution_template(this);

  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(!finalized_);

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
    const std::vector<Parameter>& parameters) {
  job_->set_id(ID<job_id_t>(*job_id_ptr_));

  IDSet<physical_data_id_t> read_set;
  read_set.insert(*from_physical_data_id_ptr_);
  job_->set_read_set(read_set);

  IDSet<physical_data_id_t> write_set;
  write_set.insert(*to_physical_data_id_ptr_);
  job_->set_write_set(write_set);
}



bool ExecutionTemplate::AddRemoteCopySendJobTemplate(RemoteCopySendCommand* command) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(!finalized_);
  job_id_t job_id = command->job_id().elem();
  JobIdPtr copy_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  JobIdPtr receive_job_id_ptr = GetExistingInnerJobIdPtr(command->receive_job_id().elem());

  static const JobIdPtr default_mega_rcr_job_id_ptr =
    JobIdPtr(new job_id_t(NIMBUS_KERNEL_JOB_ID));

  JobIdPtr mega_rcr_job_id_ptr;
  if (command->mega_rcr_job_id().elem() == NIMBUS_KERNEL_JOB_ID) {
    mega_rcr_job_id_ptr = default_mega_rcr_job_id_ptr;
  } else {
    mega_rcr_job_id_ptr = GetExistingInnerJobIdPtr(command->mega_rcr_job_id().elem());
  }

  PhyIdPtr from_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->from_physical_data_id().elem());

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr = GetExistingInnerJobIdPtr(*iter);
      before_set.insert(job_id_ptr);
    }
  }

  RemoteCopySendJobTemplate *cm =
    new RemoteCopySendJobTemplate(copy_job_id_ptr,
                                      receive_job_id_ptr,
                                      mega_rcr_job_id_ptr,
                                      from_physical_data_id_ptr,
                                      command->to_worker_id(),
                                      command->to_ip(),
                                      command->to_port(),
                                      before_set,
                                      0);

  ++copy_job_num_;
  job_templates_[job_id] = cm;

  return true;
}

bool ExecutionTemplate::AddRemoteCopyReceiveJobTemplate(RemoteCopyReceiveCommand* command) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(!finalized_);
  job_id_t job_id = command->job_id().elem();
  JobIdPtr copy_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  PhyIdPtr to_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->to_physical_data_id().elem());

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr = GetExistingInnerJobIdPtr(*iter);
      before_set.insert(job_id_ptr);
    }
  }

  RemoteCopyReceiveJobTemplate *cm =
    new RemoteCopyReceiveJobTemplate(copy_job_id_ptr,
                                         to_physical_data_id_ptr,
                                         before_set,
                                         0);

  ++copy_job_num_;
  job_templates_[job_id] = cm;

  return true;
}

bool ExecutionTemplate::AddMegaRCRJobTemplate(MegaRCRCommand* command) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(!finalized_);
  job_id_t job_id = command->job_id().elem();
  JobIdPtr copy_job_id_ptr = GetExistingInnerJobIdPtr(job_id);

  JobIdPtrList receive_job_id_ptrs;
  {
    std::vector<job_id_t>::const_iterator iter = command->receive_job_ids_p()->begin();
    for (; iter != command->receive_job_ids_p()->end(); ++iter) {
      JobIdPtr job_id_ptr = GetExistingInnerJobIdPtr(*iter);
      receive_job_id_ptrs.push_back(job_id_ptr);
    }
  }

  PhyIdPtrList to_phy_id_ptrs;
  {
    std::vector<physical_data_id_t>::const_iterator iter =
      command->to_physical_data_ids_p()->begin();
    for (; iter != command->to_physical_data_ids_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      to_phy_id_ptrs.push_back(phy_id_ptr);
    }
  }


  MegaRCRJobTemplate *cm =
    new MegaRCRJobTemplate(copy_job_id_ptr,
                               receive_job_id_ptrs,
                               to_phy_id_ptrs,
                               0);

  copy_job_num_ += command->receive_job_ids_p()->size();
  job_templates_[job_id] = cm;

  return true;
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


