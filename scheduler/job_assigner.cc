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
  * This is the base job assigner module. It provides methods required for
  * binding a job to a worker and submitting the compute and copy jobs to the
  * worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "scheduler/job_assigner.h"
#include "scheduler/load_balancer.h"

#define LB_UPDATE_RATE 100
#define JOB_ASSIGNER_THREAD_NUM 1

#define _NIMBUS_NO_LOG

namespace nimbus {

JobAssigner::JobAssigner() {
  Initialize();
}

void JobAssigner::Initialize() {
  thread_num_ = 0;
  checkpoint_id_ = NIMBUS_INIT_CHECKPOINT_ID;
  pending_assignment_ = 0;
  server_ = NULL;
  id_maker_ = NULL;
  job_manager_ = NULL;
  data_manager_ = NULL;
  load_balancer_ = NULL;
  log_.set_file_name("log_job_assigner");
  log_before_set_.set_file_name("log_before_set");
  log_assign_stamp_.set_file_name("log_assign_stamp");
}

JobAssigner::~JobAssigner() {
}

void JobAssigner::set_id_maker(IDMaker *id_maker) {
  id_maker_ = id_maker;
}

void JobAssigner::set_server(SchedulerServer *server) {
  server_ = server;
}

void JobAssigner::set_job_manager(JobManager *job_manager) {
  job_manager_ = job_manager;
}

void JobAssigner::set_data_manager(DataManager *data_manager) {
  data_manager_ = data_manager;
}

void JobAssigner::set_load_balancer(LoadBalancer *load_balancer) {
  load_balancer_ = load_balancer;
}

void JobAssigner::set_thread_num(size_t thread_num) {
  thread_num_ = thread_num;
}

void JobAssigner::set_checkpoint_id(checkpoint_id_t checkpoint_id) {
  checkpoint_id_ = checkpoint_id;
}

void JobAssigner::Run() {
  for (size_t i = 0; i < thread_num_; ++i) {
    job_assigner_threads_.push_back(
        new boost::thread(boost::bind(&JobAssigner::JobAssignerThread, this)));
  }
}


void JobAssigner::JobAssignerThread() {
  while (true) {
    JobEntry *job;
    {
      boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);

      while (job_queue_.size() == 0) {
        job_queue_cond_.wait(job_queue_lock);
      }

      JobEntryList::iterator iter = job_queue_.begin();
      job = *iter;
      job_queue_.erase(iter);
      ++pending_assignment_;
    }

    if (!AssignJob(job)) {
      dbg(DBG_ERROR, "ERROR: JobAssigner: could not assign job %lu.\n", job->job_id());
      exit(-1);
    }

    {
      boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);
      --pending_assignment_;
      job_queue_cond_.notify_all();
    }
  }
}

void JobAssigner::AssignJobs(const JobEntryList& list) {
  if (thread_num_ == 0) {
    boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);
    assert(job_queue_.size() == 0);
    job_queue_ = list;
    JobEntry *job;
    JobEntryList::iterator iter = job_queue_.begin();
    for (; iter != job_queue_.end();) {
      job = *iter;
      if (!AssignJob(job)) {
        dbg(DBG_ERROR, "ERROR: JobAssigner: could not assign job %lu.\n", job->job_id());
        exit(-1);
      }
      job_queue_.erase(iter++);
    }
  } else {
    boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);
    assert(job_queue_.size() == 0);
    job_queue_ = list;
    job_queue_cond_.notify_all();

    while ((job_queue_.size() > 0) || (pending_assignment_ > 0)) {
      job_queue_cond_.wait(job_queue_lock);
    }
  }
}

bool JobAssigner::QueryDataManagerForPatterns(
                      ComplexJobEntry* job,
                      const BindingTemplate::PatternList* patterns,
                      const BindingTemplate::PatternMetaData* patterns_meta_data,
                      std::vector<physical_data_id_t>* physical_ids) {
  physical_ids->clear();

  if (!job_manager_->ResolveJobDataVersionsForPattern(job, patterns)) {
    assert(false);
  }

  BindingTemplate::PatternMetaData::const_iterator it =
    patterns_meta_data->begin();
  BindingTemplate::PatternList::const_iterator iter;
  for (iter = patterns->begin(); iter != patterns->end();) {
    size_t ldid_count = it->first + it->second;

    logical_data_id_t ldid = (*iter)->ldid_;
    LogicalDataObject* ldo =
      const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(ldid));

    PhysicalDataList instances;
    data_manager_->AllInstances(ldo, &instances);

    // For now we assume that binding for template does not require creating
    // new instances -omidm
    assert(instances.size() >= ldid_count);

    data_version_t base_version;
    if (!job->vmap_read()->query_entry(ldid, &base_version)) {
      assert(false);
    }

    for (size_t i = 0; i < ldid_count; ++i) {
      BindingTemplate::PatternEntry *pattern = (*iter);
      bool found = false;

      data_version_t version = pattern->version_diff_from_base_ + base_version;

      PhysicalDataList::iterator pi = instances.begin();
      for (; pi != instances.end(); ++pi) {
        if ((pi->worker() == pattern->worker_id_) &&
            ((pi->version() == version) ||
             (pattern->version_type_ == BindingTemplate::WILD_CARD))) {
          found = true;
          physical_ids->push_back(pi->id());
          instances.erase(pi);
          break;
        }
      }

      assert(found);
      ++iter;
    }
    ++it;
  }

  return true;
}

bool JobAssigner::UpdateDataManagerByPatterns(
                      ComplexJobEntry* job,
                      const BindingTemplate::PatternList* patterns,
                      const std::vector<physical_data_id_t>* physical_ids) {
  assert(patterns->size() == physical_ids->size());

  IDSet<job_id_t> list_job_read;
  list_job_read.insert(job->job_id());
  job_id_t last_job_write = job->job_id();

  size_t physical_ids_idx = 0;
  BindingTemplate::PatternList::const_iterator iter = patterns->begin();
  for (; iter != patterns->end(); ++iter) {
    BindingTemplate::PatternEntry *pe = *iter;
    data_version_t base_version;
    if (!job->vmap_read()->query_entry(pe->ldid_, &base_version)) {
      assert(false);
    }
    data_version_t version = base_version + pe->version_diff_from_base_;

    if (!data_manager_->UpdateVersionAndAccessRecord(pe->ldid_,
                                                     physical_ids->operator[](physical_ids_idx),
                                                     version,
                                                     list_job_read,
                                                     last_job_write)) {
      assert(false);
    }

    ++physical_ids_idx;
  }

  return true;
}

bool JobAssigner::AssignComplexJob(ComplexJobEntry *job) {
  BindingTemplate *bt = job->binding_template();
  size_t copy_job_num = bt->copy_job_num();
  size_t compute_job_num = bt->compute_job_num();

  assert(compute_job_num == job->inner_job_ids_p()->size());

  std::vector<job_id_t> copy_job_ids;
  id_maker_->GetNewJobID(&copy_job_ids, copy_job_num);

  std::vector<physical_data_id_t> physical_ids;

  QueryDataManagerForPatterns(job,
                              bt->entry_pattern_list_p(),
                              bt->patterns_meta_data_p(),
                              &physical_ids);

  bt->Instantiate(job->inner_job_ids(),
                  job->parameters(),
                  copy_job_ids,
                  physical_ids,
                  server_);

  UpdateDataManagerByPatterns(job,
                              bt->end_pattern_list_p(),
                              &physical_ids);

  return true;
}

bool JobAssigner::AssignJob(JobEntry *job) {
  if (job->job_type() == JOB_CMPX) {
    ComplexJobEntry *xj = reinterpret_cast<ComplexJobEntry*>(job);
    return AssignComplexJob(xj);
  }

  SchedulerWorker* worker = job->assigned_worker();
  log_assign_.log_StartTimer();
  log_prepare_.log_ResetTimer();
  log_job_manager_.log_ResetTimer();
  log_data_manager_.log_ResetTimer();
  log_version_manager_.log_ResetTimer();

  log_version_manager_.log_ResumeTimer();
  job_manager_->ResolveJobDataVersions(job);
  assert(job->versioned());
  if (job->memoize()) {
    job_manager_->MemoizeVersionsForTemplate(job);
  }
  log_version_manager_.log_StopTimer();


  if (job->checkpoint_id() > NIMBUS_INIT_CHECKPOINT_ID) {
    SaveJobContextForCheckpoint(job);
  }


  log_prepare_.log_ResumeTimer();
  bool prepared_data = true;
  IDSet<logical_data_id_t>::ConstIter it;
  for (it = job->union_set_p()->begin(); it != job->union_set_p()->end(); ++it) {
    if (!PrepareDataForJobAtWorker(job, worker, *it)) {
      prepared_data = false;
      break;
    }
  }
  log_prepare_.log_StopTimer();


  if (prepared_data) {
    PrintLog(job);

    // BINDING MEMOIZE - omidm
    if (job->memoize_binding()) {
      assert(job->job_type() == JOB_SHDW);
      ShadowJobEntry *sj = reinterpret_cast<ShadowJobEntry*>(job);
      BindingTemplate *bt = sj->binding_template();

      ID<job_id_t> job_id(job->job_id());
      ID<job_id_t> future_job_id(job->future_job_id());
      IDSet<physical_data_id_t> read_set, write_set;
      // TODO(omidm): check the return value of the following methods.
      job->GetPhysicalReadSet(&read_set);
      job->GetPhysicalWriteSet(&write_set);
      ComputeJobCommand cm(job->job_name(),
                           job_id,
                           read_set,
                           write_set,
                           job->before_set(),
                           job->after_set(),
                           future_job_id,
                           job->sterile(),
                           job->region(),
                           job->params());

      bt->AddComputeJobCommand(&cm, worker->worker_id());
    }
    // BINDING MEMOIZE - omidm


    log_job_manager_.log_ResumeTimer();
    job_manager_->UpdateJobBeforeSet(job);
    log_job_manager_.log_StopTimer();
    SendComputeJobToWorker(worker, job);

    // BINDING MEMOIZE - omidm
    if (job->memoize_binding()) {
      assert(job->job_type() == JOB_SHDW);
      ShadowJobEntry *sj = reinterpret_cast<ShadowJobEntry*>(job);
      ComplexJobEntry *xj = sj->complex_job();

      if (job->to_finalize_binding_template()) {
        job->binding_template()->Finalize(xj->inner_job_ids());
      }
    }
    // BINDING MEMOIZE - omidm

    log_job_manager_.log_ResumeTimer();
    job_manager_->NotifyJobAssignment(job);
    load_balancer_->NotifyJobAssignment(job);
    log_job_manager_.log_StopTimer();

#ifndef _NIMBUS_NO_LOG
    log_assign_.log_StopTimer();
    std::cout << "jobmanager: " << log_job_manager_.timer() << " n: " << job->job_name() << std::endl; // NOLINT
    std::cout << "datamanager: " << log_data_manager_.timer() << " n: " << job->job_name() << std::endl; // NOLINT
    std::cout << "versionmanager: " << log_version_manager_.timer() << " n: " << job->job_name() << std::endl; // NOLINT
    std::cout << "prepare: " << log_prepare_.timer() << " n: " << job->job_name() << std::endl; // NOLINT
    std::cout << "assign: " << log_assign_.timer() << " n: " << job->job_name() << std::endl; // NOLINT
#endif

    return true;
  }

  return false;
}

void JobAssigner::PrintLog(JobEntry *job) {
  char buff[LOG_MAX_BUFF_SIZE];

  snprintf(buff, sizeof(buff), "id: %lu bs: %s.",
      job->job_id(), job->before_set_p()->ToString().c_str());
  log_before_set_.log_WriteToFile(std::string(buff), LOG_INFO);

  snprintf(buff, sizeof(buff), "%10.9f id: %lu n: %s.",
      Log::GetRawTime(), job->job_id(), job->job_name().c_str());
  log_assign_stamp_.log_WriteToFile(std::string(buff), LOG_INFO);
}

bool JobAssigner::PrepareDataForJobAtWorker(JobEntry* job,
                                            SchedulerWorker* worker,
                                            logical_data_id_t l_id) {
  bool reading = job->read_set_p()->contains(l_id);
  bool writing = job->write_set_p()->contains(l_id);
  assert(reading || writing);

  log_data_manager_.log_ResumeTimer();
  LogicalDataObject* ldo =
    const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(l_id));
  log_data_manager_.log_StopTimer();

  boost::unique_lock<boost::mutex> lock(ldo->mutex());

  data_version_t version;
  if (reading) {
    if (!job->vmap_read()->query_entry(l_id, &version)) {
      dbg(DBG_ERROR, "ERROR: logical id %lu is not versioned in the read context of %s.\n",
          l_id, job->job_name().c_str());
      assert(false);
    }
  }

  // BINDING MEMOIZE - omidm
  BindingTemplate *bt = NULL;
  data_version_t version_diff;
  bool memoize_binding = false;
  if (job->memoize_binding()) {
    assert(job->job_type() == JOB_SHDW);
    ShadowJobEntry *sj = reinterpret_cast<ShadowJobEntry*>(job);
    bt = sj->binding_template();

    if (!sj->vmap_read_diff()->query_entry(l_id, &version_diff)) {
      assert(false);
    }

    memoize_binding = true;
  }
  // BINDING MEMOIZE - omidm

  // Just for checking
  data_version_t unused_version;
  if (writing) {
    if (!job->vmap_write()->query_entry(l_id, &unused_version)) {
      dbg(DBG_ERROR, "ERROR: logical id %lu is not versioned in the write context of %s.\n",
          l_id, job->job_name().c_str());
      assert(false);
    }
  }

  if (!reading) {
    PhysicalData target_instance;
    GetFreeDataAtWorker(worker, ldo, &target_instance);

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      bt->TrackDataObject(worker->worker_id(),
                          l_id,
                          target_instance.id(),
                          BindingTemplate::WILD_CARD,
                          0);
    }
    // BINDING MEMOIZE - omidm

    log_job_manager_.log_ResumeTimer();
    if (job_manager_->CausingUnwantedSerialization(job, l_id, target_instance, memoize_binding)) {
      dbg(DBG_ERROR, "Why serializing!!\n");
      dbg(DBG_SCHED, "Causing unwanted serialization for data %lu.\n", l_id);
      assert(false);
    }
    log_job_manager_.log_StopTimer();

    AllocateLdoInstanceToJob(job, ldo, target_instance);

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      assert(writing);
      bt->UpdateDataObject(target_instance.id(),
                           version_diff + 1);
    }
    // BINDING MEMOIZE - omidm

    return true;
  }

  PhysicalDataList instances_at_worker;
  log_data_manager_.log_ResumeTimer();
  data_manager_->InstancesByWorkerAndVersion(
      ldo, worker->worker_id(), version, &instances_at_worker);
  log_data_manager_.log_StopTimer();

  JobEntryList list;
  VersionedLogicalData vld(l_id, version);
  log_version_manager_.log_ResumeTimer();
  job_manager_->GetJobsNeedDataVersion(&list, vld);
  log_version_manager_.log_StopTimer();
  assert(list.size() >= 1);
  bool writing_needed_version = (list.size() > 1) && writing;


  if (instances_at_worker.size() > 1) {
    PhysicalData target_instance;

    bool found = false;
    PhysicalDataList::iterator iter;
    for (iter = instances_at_worker.begin(); iter != instances_at_worker.end(); iter++) {
      log_job_manager_.log_ResumeTimer();
      if (!job_manager_->CausingUnwantedSerialization(job, l_id, *iter, memoize_binding)) {
        log_job_manager_.log_StopTimer();
        target_instance = *iter;
        found = true;
        break;
      }
      log_job_manager_.log_StopTimer();
    }

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      if (found) {
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            target_instance.id(),
                            BindingTemplate::REGULAR,
                            version_diff);
      }
    }
    // BINDING MEMOIZE - omidm

    if (!found) {
      dbg(DBG_SCHED, "Avoiding unwanted serialization for data %lu (1).\n", l_id);
      GetFreeDataAtWorker(worker, ldo, &target_instance);

      // BINDING MEMOIZE - omidm
      if (memoize_binding) {
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            target_instance.id(),
                            BindingTemplate::WILD_CARD,
                            0);
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            (*instances_at_worker.begin()).id(),
                            BindingTemplate::REGULAR,
                            version_diff);
        bt->UpdateDataObject(target_instance.id(),
                             version_diff);
      }
      // BINDING MEMOIZE - omidm

      // BINDING MEMOIZE - omidm
      LocalCopyData(worker, ldo, &(*instances_at_worker.begin()), &target_instance, bt);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      if (writing) {
        bt->UpdateDataObject(target_instance.id(),
                             version_diff + 1);
      }
    }
    // BINDING MEMOIZE - omidm

    return true;
  }


  if ((instances_at_worker.size() == 1) && !writing_needed_version) {
    PhysicalData target_instance;

    log_job_manager_.log_ResumeTimer();
    if (!job_manager_->CausingUnwantedSerialization(job, l_id, *instances_at_worker.begin(), memoize_binding)) { //NOLINT
      log_job_manager_.log_StopTimer();
      target_instance = *instances_at_worker.begin();

      // BINDING MEMOIZE - omidm
      if (memoize_binding) {
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            target_instance.id(),
                            BindingTemplate::REGULAR,
                            version_diff);
      }
      // BINDING MEMOIZE - omidm
    } else {
      log_job_manager_.log_StopTimer();
      dbg(DBG_SCHED, "Avoiding unwanted serialization for data %lu (2).\n", l_id);
      GetFreeDataAtWorker(worker, ldo, &target_instance);

      // BINDING MEMOIZE - omidm
      if (memoize_binding) {
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            target_instance.id(),
                            BindingTemplate::WILD_CARD,
                            0);
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            (*instances_at_worker.begin()).id(),
                            BindingTemplate::REGULAR,
                            version_diff);
        bt->UpdateDataObject(target_instance.id(),
                             version_diff);
      }
      // BINDING MEMOIZE - omidm

      // BINDING MEMOIZE - omidm
      LocalCopyData(worker, ldo, &(*instances_at_worker.begin()), &target_instance, bt);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      if (writing) {
        bt->UpdateDataObject(target_instance.id(),
                             version_diff + 1);
      }
    }
    // BINDING MEMOIZE - omidm

    return true;
  }


  if ((instances_at_worker.size() == 1) && writing_needed_version) {
    PhysicalData target_instance;

    log_job_manager_.log_ResumeTimer();
    if (!job_manager_->CausingUnwantedSerialization(job, l_id, *instances_at_worker.begin(), memoize_binding)) { // NOLINT
      log_job_manager_.log_StopTimer();
      target_instance = *instances_at_worker.begin();
      PhysicalData copy_data;
      GetFreeDataAtWorker(worker, ldo, &copy_data);

      // BINDING MEMOIZE - omidm
      if (memoize_binding) {
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            target_instance.id(),
                            BindingTemplate::REGULAR,
                            version_diff);
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            copy_data.id(),
                            BindingTemplate::WILD_CARD,
                            0);
        bt->UpdateDataObject(copy_data.id(),
                             version_diff);
      }
      // BINDING MEMOIZE - omidm

      // BINDING MEMOIZE - omidm
      LocalCopyData(worker, ldo, &target_instance, &copy_data, bt);
    } else {
      log_job_manager_.log_StopTimer();
      dbg(DBG_SCHED, "Avoiding unwanted serialization for data %lu (3).\n", l_id);
      GetFreeDataAtWorker(worker, ldo, &target_instance);

      // BINDING MEMOIZE - omidm
      if (memoize_binding) {
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            target_instance.id(),
                            BindingTemplate::WILD_CARD,
                            0);
        bt->TrackDataObject(worker->worker_id(),
                            l_id,
                            (*instances_at_worker.begin()).id(),
                            BindingTemplate::REGULAR,
                            version_diff);
        bt->UpdateDataObject(target_instance.id(),
                             version_diff);
      }
      // BINDING MEMOIZE - omidm

      // BINDING MEMOIZE - omidm
      LocalCopyData(worker, ldo, &(*instances_at_worker.begin()), &target_instance, bt);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      if (writing) {
        bt->UpdateDataObject(target_instance.id(),
                             version_diff + 1);
      }
    }
    // BINDING MEMOIZE - omidm

    return true;
  }


  if ((instances_at_worker.size() == 0) && (version == NIMBUS_INIT_DATA_VERSION)) {
    // BINDING MEMOIZE - omidm
    assert(!memoize_binding);
    // BINDING MEMOIZE - omidm

    PhysicalData created_data;
    CreateDataAtWorker(worker, ldo, &created_data);

    AllocateLdoInstanceToJob(job, ldo, created_data);
    return true;
  }

  PhysicalDataList instances_in_system;
  log_data_manager_.log_ResumeTimer();
  data_manager_->InstancesByVersion(ldo, version, &instances_in_system);
  log_data_manager_.log_StopTimer();

  if ((instances_at_worker.size() == 0) && (instances_in_system.size() >= 1)) {
    PhysicalData from_instance = *instances_in_system.begin();
    worker_id_t sender_id = from_instance.worker();
    SchedulerWorker* worker_sender;
    if (!server_->GetSchedulerWorkerById(worker_sender, sender_id)) {
      dbg(DBG_ERROR, "ERROR: could not find worker with id %lu.\n", sender_id);
      exit(-1);
    }

    PhysicalData target_instance;
    GetFreeDataAtWorker(worker, ldo, &target_instance);

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      bt->TrackDataObject(worker->worker_id(),
                          l_id,
                          target_instance.id(),
                          BindingTemplate::WILD_CARD,
                          0);
      bt->TrackDataObject(sender_id,
                          l_id,
                          from_instance.id(),
                          BindingTemplate::REGULAR,
                          version_diff);
      bt->UpdateDataObject(target_instance.id(),
                           version_diff);
    }
    // BINDING MEMOIZE - omidm

    // BINDING MEMOIZE - omidm
    RemoteCopyData(worker_sender, worker, ldo, &from_instance, &target_instance, bt);

    AllocateLdoInstanceToJob(job, ldo, target_instance);

    // BINDING MEMOIZE - omidm
    if (memoize_binding) {
      if (writing) {
        bt->UpdateDataObject(target_instance.id(),
                             version_diff + 1);
      }
    }
    // BINDING MEMOIZE - omidm

    return true;
  }

  assert(NIMBUS_FAULT_TOLERANCE_ACTIVE);

  // If you are here, you may find the version in checkpoint.
  dbg(DBG_WARN, "WARNING: looking in to checkpoint to load the data.\n");

  WorkerHandleList handles;
  job_manager_->GetHandleToLoadData(checkpoint_id_, l_id, version, &handles);

  if (handles.size() > 0) {
    std::string handle = handles.begin()->second;
    WorkerHandleList::iterator iter = handles.begin();
    for (; iter != handles.end(); ++iter) {
      if (iter->first == worker->worker_id()) {
        handle = iter->second;
        break;
      }
    }

    PhysicalData target_instance;
    GetFreeDataAtWorker(worker, ldo, &target_instance);
    LoadData(worker, ldo, version, &target_instance, handle);

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }

  dbg(DBG_ERROR, "ERROR: the version (%lu) of logical data %s (%lu) needed for job %s (%lu) does not exist.\n", // NOLINT
      version, ldo->variable().c_str(), l_id, job->job_name().c_str(), job->job_id());
  assert(instances_in_system.size() >= 1);

  return false;
}


bool JobAssigner::AllocateLdoInstanceToJob(JobEntry* job,
                                           LogicalDataObject* ldo,
                                           PhysicalData pd) {
  assert(job->versioned());
  PhysicalData pd_new = pd;

  if (job->write_set_p()->contains(ldo->id())) {
    data_version_t v_out;
    job->vmap_write()->query_entry(ldo->id(), &v_out);
    pd_new.set_version(v_out);
    pd_new.set_last_job_write(job->job_id());
    pd_new.clear_list_job_read();
    job->before_set_p()->insert(pd.list_job_read());
    job->before_set_p()->insert(pd.last_job_write());
  }

  if (job->read_set_p()->contains(ldo->id())) {
    data_version_t v_in;
    job->vmap_read()->query_entry(ldo->id(), &v_in);
    assert(v_in == pd.version());
    pd_new.add_to_list_job_read(job->job_id());
    job->before_set_p()->insert(pd.last_job_write());
  }

  job->set_physical_table_entry(ldo->id(), pd.id());

  log_data_manager_.log_ResumeTimer();
  data_manager_->UpdatePhysicalInstance(ldo, pd, pd_new);
  log_data_manager_.log_StopTimer();

  return true;
}


bool JobAssigner::SaveJobContextForCheckpoint(JobEntry *job) {
  assert(job->checkpoint_id() > NIMBUS_INIT_CHECKPOINT_ID);

  job_manager_->ResolveEntireContextForJob(job);
  job_manager_->CompleteJobForCheckpoint(job->checkpoint_id(), job);
  VersionMap::ConstIter iter = job->vmap_read()->content_p()->begin();
  for (; iter != job->vmap_read()->content_p()->end(); ++iter) {
    logical_data_id_t ldid = iter->first;
    data_version_t version = iter->second;
    LogicalDataObject* ldo =
      const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(ldid));

    PhysicalDataList instances_in_system;
    data_manager_->InstancesByVersion(ldo, version, &instances_in_system);


    if (instances_in_system.size() >= 1) {
      PhysicalDataList::iterator it = instances_in_system.begin();
      for (; it != instances_in_system.end(); ++it) {
        worker_id_t worker_id = it->worker();
        SchedulerWorker* worker;
        if (!server_->GetSchedulerWorkerById(worker, worker_id)) {
          dbg(DBG_ERROR, "ERROR: could not find worker with id %lu.\n", worker_id);
          exit(-1);
        }
        SaveData(worker, ldo, &(*it), job->checkpoint_id());
      }

      continue;
    }

    // If you are here, you may find the version in checkpoint.
    dbg(DBG_WARN, "WARNING: looking in to checkpoint to load the data for next checkpoint!.\n");

    WorkerHandleList handles;
    job_manager_->GetHandleToLoadData(checkpoint_id_, ldid, version, &handles);

    if (handles.size() > 0) {
      WorkerHandleList::iterator iter = handles.begin();
      for (; iter != handles.end(); ++iter) {
        worker_id_t worker_id = iter->first;
        std::string handle = iter->second;
        SchedulerWorker* worker;
        if (!server_->GetSchedulerWorkerById(worker, worker_id)) {
          // Worker could be down so you may not find it there
          worker = *(server_->workers()->begin());
        }

        PhysicalData target_instance;
        GetFreeDataAtWorker(worker, ldo, &target_instance);
        LoadData(worker, ldo, version, &target_instance, handle);
        SaveData(worker, ldo, &target_instance, job->checkpoint_id());
      }

      continue;
    }


    dbg(DBG_ERROR, "ERROR: version (%lu) of logical data %s (%lu) needed for CHECKPOINT job %s (%lu) does not exist.\n", // NOLINT
        version, ldo->variable().c_str(), ldid, job->job_name().c_str(), job->job_id());
    assert(instances_in_system.size() >= 1);
  }

  return true;
}

bool JobAssigner::LoadData(SchedulerWorker* worker,
                           LogicalDataObject* ldo,
                           data_version_t version,
                           PhysicalData* to_data,
                           std::string handle) {
  assert(worker->worker_id() == to_data->worker());

  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 1);
  IDSet<job_id_t> before;

  // Update data table.
  PhysicalData to_data_new = *to_data;
  to_data_new.set_version(version);
  to_data_new.set_last_job_write(j[0]);
  to_data_new.clear_list_job_read();
  data_manager_->UpdatePhysicalInstance(ldo, *to_data, to_data_new);

  // find the before set.
  before.insert(to_data->list_job_read());
  before.insert(to_data->last_job_write());
  job_manager_->UpdateBeforeSet(&before);

  // send the load command to worker.
  LoadDataCommand cm(ID<job_id_t>(j[0]),
                     ID<physical_data_id_t>(to_data->id()),
                     handle,
                     before);
  server_->SendCommand(worker, &cm);

  *to_data = to_data_new;

  return true;
}


bool JobAssigner::SaveData(SchedulerWorker* worker,
                           LogicalDataObject* ldo,
                           PhysicalData* from_data,
                           checkpoint_id_t checkpoint_id) {
  assert(worker->worker_id() == from_data->worker());

  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 1);
  IDSet<job_id_t> before;

  // Update data table.
  PhysicalData from_data_new = *from_data;
  from_data_new.add_to_list_job_read(j[0]);
  data_manager_->UpdatePhysicalInstance(ldo, *from_data, from_data_new);

  // find the before set.
  before.insert(from_data->last_job_write());
  job_manager_->UpdateBeforeSet(&before);

  // Add the entry for the checkpoint
  job_manager_->AddSaveDataJobToCheckpoint(checkpoint_id,
                                           j[0],
                                           ldo->id(),
                                           from_data->version(),
                                           worker->worker_id());

  // send save data command to worker.
  SaveDataCommand cm(ID<job_id_t>(j[0]),
                     ID<physical_data_id_t>(from_data->id()),
                     ID<checkpoint_id_t>(checkpoint_id),
                     before);
  server_->SendCommand(worker, &cm);

  *from_data = from_data_new;

  return true;
}

size_t JobAssigner::GetObsoleteLdoInstancesAtWorker(SchedulerWorker* worker,
                                                    LogicalDataObject* ldo,
                                                    PhysicalDataList* dest) {
  size_t count = 0;
  dest->clear();
  PhysicalDataList pv;
  log_data_manager_.log_ResumeTimer();
  data_manager_->InstancesByWorker(ldo, worker->worker_id(), &pv);
  log_data_manager_.log_StopTimer();
  PhysicalDataList::iterator iter = pv.begin();
  for (; iter != pv.end(); ++iter) {
    JobEntryList list;
    VersionedLogicalData vld(ldo->id(), iter->version());
    log_version_manager_.log_ResumeTimer();
    if (job_manager_->GetJobsNeedDataVersion(&list, vld) == 0) {
      dest->push_back(*iter);
      ++count;
    }
    log_version_manager_.log_StopTimer();
  }
  return count;
}

bool JobAssigner::CreateDataAtWorker(SchedulerWorker* worker,
                                     LogicalDataObject* ldo,
                                     PhysicalData* created_data) {
  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 1);
  std::vector<physical_data_id_t> d;
  id_maker_->GetNewPhysicalDataID(&d, 1);
  IDSet<job_id_t> before;

  // Update the job table.
  // JobEntry *job = job_manager_->AddCreateDataJobEntry(j[0]);

  // Update data table.
  IDSet<job_id_t> list_job_read;
  list_job_read.insert(j[0]);  // if other job wants to write, waits for creation.
  PhysicalData p(d[0], worker->worker_id(), NIMBUS_INIT_DATA_VERSION, list_job_read, j[0]);
  log_data_manager_.log_ResumeTimer();
  data_manager_->AddPhysicalInstance(ldo, p);
  log_data_manager_.log_StopTimer();

  // send the create command to worker.
  log_job_manager_.log_ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.log_StopTimer();
  CreateDataCommand cm(ID<job_id_t>(j[0]),
                       ldo->variable(),
                       ID<logical_data_id_t>(ldo->id()),
                       ID<physical_data_id_t>(d[0]),
                       before);
  server_->SendCommand(worker, &cm);

  // Notify assignment to job manager.
  // job->set_assigned_worker(worker);
  // job->set_before_set(before);
  // job_manager_->NotifyJobAssignment(job);

  *created_data = p;

  return true;
}

bool JobAssigner::RemoteCopyData(SchedulerWorker* from_worker,
                                 SchedulerWorker* to_worker,
                                 LogicalDataObject* ldo,
                                 PhysicalData* from_data,
                                 PhysicalData* to_data,
                                 BindingTemplate *bt) {
  assert(from_worker->worker_id() == from_data->worker());
  assert(to_worker->worker_id() == to_data->worker());

  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 2);
  job_id_t receive_id = j[0];
  job_id_t send_id = j[1];
  IDSet<job_id_t> before;
  // JobEntry *job;

  // Receive part

  // Update the job table.
  // job = job_manager_->AddRemoteCopyReceiveJobEntry(receive_id);

  // Update data table.
  PhysicalData to_data_new = *to_data;
  to_data_new.set_version(from_data->version());
  to_data_new.set_last_job_write(receive_id);
  to_data_new.clear_list_job_read();
  log_data_manager_.log_ResumeTimer();
  data_manager_->UpdatePhysicalInstance(ldo, *to_data, to_data_new);
  log_data_manager_.log_StopTimer();

  // send remote copy receive job to worker.
  before.clear();
  before.insert(to_data->list_job_read());
  before.insert(to_data->last_job_write());

  // BINDING MEMOIZE - omidm
  if (bt) {
    RemoteCopyReceiveCommand cm_r(ID<job_id_t>(receive_id),
                                  ID<physical_data_id_t>(to_data->id()),
                                  before);
    bt->AddRemoteCopyReceiveCommand(&cm_r, to_worker->worker_id());
  }
  // BINDING MEMOIZE - omidm

  log_job_manager_.log_ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.log_StopTimer();
  RemoteCopyReceiveCommand cm_r(ID<job_id_t>(receive_id),
                                ID<physical_data_id_t>(to_data->id()),
                                before);
  server_->SendCommand(to_worker, &cm_r);

  // Notify assignment to job manager.
  // job->set_assigned_worker(to_worker);
  // job->set_before_set(before);
  // job_manager_->NotifyJobAssignment(job);


  // Send Part.

  // Update the job table.
  // job = job_manager_->AddRemoteCopySendJobEntry(send_id);

  // Update data table.
  PhysicalData from_data_new = *from_data;
  from_data_new.add_to_list_job_read(send_id);
  log_data_manager_.log_ResumeTimer();
  data_manager_->UpdatePhysicalInstance(ldo, *from_data, from_data_new);
  log_data_manager_.log_StopTimer();

  // send remote copy send command to worker.
  before.clear();
  before.insert(from_data->last_job_write());

  // BINDING MEMOIZE - omidm
  if (bt) {
    RemoteCopySendCommand cm_s(ID<job_id_t>(send_id),
                               ID<job_id_t>(receive_id),
                               ID<physical_data_id_t>(from_data->id()),
                               ID<worker_id_t>(to_worker->worker_id()),
                               to_worker->ip(),
                               ID<port_t>(to_worker->port()),
                               before);
    bt->AddRemoteCopySendCommand(&cm_s, from_worker->worker_id());
  }
  // BINDING MEMOIZE - omidm

  log_job_manager_.log_ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.log_StopTimer();
  RemoteCopySendCommand cm_s(ID<job_id_t>(send_id),
                             ID<job_id_t>(receive_id),
                             ID<physical_data_id_t>(from_data->id()),
                             ID<worker_id_t>(to_worker->worker_id()),
                             to_worker->ip(),
                             ID<port_t>(to_worker->port()),
                             before);
  server_->SendCommand(from_worker, &cm_s);

  // Notify assignment to job manager.
  // job->set_assigned_worker(from_worker);
  // job->set_before_set(before);
  // job_manager_->NotifyJobAssignment(job);


  *from_data = from_data_new;
  *to_data = to_data_new;

  return true;
}

bool JobAssigner::LocalCopyData(SchedulerWorker* worker,
                                LogicalDataObject* ldo,
                                PhysicalData* from_data,
                                PhysicalData* to_data,
                                BindingTemplate *bt) {
  assert(worker->worker_id() == from_data->worker());
  assert(worker->worker_id() == to_data->worker());

  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 1);
  IDSet<job_id_t> before;

  // Update the job table.
  // JobEntry *job = job_manager_->AddLocalCopyJobEntry(j[0]);

  // Update data table.
  PhysicalData from_data_new = *from_data;
  from_data_new.add_to_list_job_read(j[0]);
  log_data_manager_.log_ResumeTimer();
  data_manager_->UpdatePhysicalInstance(ldo, *from_data, from_data_new);
  log_data_manager_.log_StopTimer();

  PhysicalData to_data_new = *to_data;
  to_data_new.set_version(from_data->version());
  to_data_new.set_last_job_write(j[0]);
  to_data_new.clear_list_job_read();
  log_data_manager_.log_ResumeTimer();
  data_manager_->UpdatePhysicalInstance(ldo, *to_data, to_data_new);
  log_data_manager_.log_StopTimer();

  // send local copy command to worker.
  before.insert(to_data->list_job_read());
  before.insert(to_data->last_job_write());
  before.insert(from_data->last_job_write());

  // BINDING MEMOIZE - omidm
  if (bt) {
    LocalCopyCommand cm_c(ID<job_id_t>(j[0]),
                          ID<physical_data_id_t>(from_data->id()),
                          ID<physical_data_id_t>(to_data->id()),
                          before);
    bt->AddLocalCopyCommand(&cm_c, worker->worker_id());
  }
  // BINDING MEMOIZE - omidm

  log_job_manager_.log_ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.log_StopTimer();
  LocalCopyCommand cm_c(ID<job_id_t>(j[0]),
                        ID<physical_data_id_t>(from_data->id()),
                        ID<physical_data_id_t>(to_data->id()),
                        before);
  server_->SendCommand(worker, &cm_c);

  // Notify assignment to job manager.
  // job->set_assigned_worker(worker);
  // job->set_before_set(before);
  // job_manager_->NotifyJobAssignment(job);

  *from_data = from_data_new;
  *to_data = to_data_new;

  return true;
}

bool JobAssigner::GetFreeDataAtWorker(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalData* free_data) {
  PhysicalDataList obsolete_instances;
  if (GetObsoleteLdoInstancesAtWorker(worker, ldo, &obsolete_instances) > 0) {
    *free_data = *obsolete_instances.begin();
    return true;
  }

  // Since there are no obsoletes, go ahead and create a new one.
  return CreateDataAtWorker(worker, ldo, free_data);
}

bool JobAssigner::SendComputeJobToWorker(SchedulerWorker* worker,
                                         JobEntry* job) {
  if ((job->job_type() == JOB_COMP) || (job->job_type() == JOB_SHDW)) {
    ID<job_id_t> job_id(job->job_id());
    ID<job_id_t> future_job_id(job->future_job_id());
    IDSet<physical_data_id_t> read_set, write_set;
    // TODO(omidm): check the return value of the following methods.
    job->GetPhysicalReadSet(&read_set);
    job->GetPhysicalWriteSet(&write_set);
    ComputeJobCommand cm(job->job_name(),
                         job_id,
                         read_set,
                         write_set,
                         job->before_set(),
                         job->after_set(),
                         future_job_id,
                         job->sterile(),
                         job->region(),
                         job->params());
    dbg(DBG_SCHED, "Sending compute job %lu to worker %lu.\n", job->job_id(), worker->worker_id());
    server_->SendCommand(worker, &cm);

    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: Job with id %lu is not a compute job.\n", job->job_id());
    return false;
  }
}

}  // namespace nimbus
