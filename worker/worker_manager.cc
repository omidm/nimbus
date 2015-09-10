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
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <pthread.h>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "shared/nimbus.h"
#include "shared/profiler_malloc.h"
#include "worker/core_config.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_computation.h"
#include "worker/worker_thread_auxiliary.h"
#include "worker/worker_thread_monitor.h"

#include "worker/worker_manager.h"

namespace nimbus {

uint64_t WorkerManager::inside_job_parallism = 0;
uint64_t WorkerManager::across_job_parallism = 0;

WorkerManager::WorkerManager() {
  pthread_mutex_init(&scheduling_needed_lock_, NULL);
  pthread_cond_init(&scheduling_needed_cond_, NULL);
  scheduling_needed_ = false;

  pthread_mutex_init(&scheduling_critical_section_lock_, NULL);

  pthread_mutex_init(&computation_job_queue_lock_, NULL);
  pthread_mutex_init(&auxiliary_job_queue_lock_, NULL);

  pthread_mutex_init(&local_job_done_list_lock_, NULL);

  if (across_job_parallism <= 0) {
    computation_thread_num = 1;
  } else {
    computation_thread_num = across_job_parallism;
  }
  dispatched_computation_job_count_ = 0;
  ready_jobs_count_ = 0;
  ongoing_parallelism_ = 0;
}

WorkerManager::~WorkerManager() {
  pthread_mutex_destroy(&scheduling_needed_lock_);
  pthread_cond_destroy(&scheduling_needed_cond_);

  pthread_mutex_destroy(&scheduling_critical_section_lock_);

  pthread_mutex_destroy(&computation_job_queue_lock_);
  pthread_mutex_destroy(&auxiliary_job_queue_lock_);
}


bool WorkerManager::ClassifyAndAddJob(Job* job) {
  bool result;
  if (use_auxiliary_thread_ &&
      (dynamic_cast<RemoteCopySendJob*>(job) ||  // NOLINT
      dynamic_cast<RemoteCopyReceiveJob*>(job) ||  // NOLINT
      dynamic_cast<MegaRCRJob*>(job) ||  // NOLINT
      dynamic_cast<LocalCopyJob*>(job) ||  // NOLINT
      dynamic_cast<CreateDataJob*>(job))) {  // NOLINT
    pthread_mutex_lock(&auxiliary_job_queue_lock_);
    auxiliary_job_list_.push_back(job);
    pthread_mutex_unlock(&auxiliary_job_queue_lock_);
    result = false;
  } else {
    pthread_mutex_lock(&computation_job_queue_lock_);
    PRIORITY_TYPE priority = LOW_PRIORITY;
    if (dynamic_cast<RemoteCopySendJob*>(job) ||     // NOLINT
        dynamic_cast<RemoteCopyReceiveJob*>(job) ||  // NOLINT
        dynamic_cast<MegaRCRJob*>(job) ||  // NOLINT
        dynamic_cast<LocalCopyJob*>(job)) {          // NOLINT
      priority = HIGH_PRIORITY;
    }
    computation_job_list_.push_back(job, priority);
    ++ready_jobs_count_;
    pthread_mutex_unlock(&computation_job_queue_lock_);
    result = true;
  }
  return result;
}

bool WorkerManager::PushJob(Job* job) {
  // TODO(quhang): when a job is dispatched.
  worker_->StatDispatchJob(1);
  if (ClassifyAndAddJob(job)) {
    TriggerScheduling();
  }
  return true;
}

bool WorkerManager::PushJobList(std::list<Job*>* job_list) {
  bool result = false;
  for (std::list<Job*>::iterator iter = job_list->begin();
       iter != job_list->end();
       ++iter) {
    if (ClassifyAndAddJob(*iter)) {
      result = true;
    }
  }
  worker_->StatDispatchJob(job_list->size());
  if (result) {
    TriggerScheduling();
  }
  return true;
}

bool WorkerManager::PullAuxiliaryJobs(WorkerThreadAuxiliary* worker_thread,
                       std::list<Job*>* job_list) {
  pthread_mutex_lock(&auxiliary_job_queue_lock_);
  job_list->splice(job_list->end(), auxiliary_job_list_);
  pthread_mutex_unlock(&auxiliary_job_queue_lock_);
  return true;
}

bool WorkerManager::FinishJob(Job* job) {
  // pthread_mutex_lock(&local_job_done_list_lock_);
  // local_job_done_list_.push_back(job);
  // pthread_mutex_unlock(&local_job_done_list_lock_);

  // worker_->StatEndJob(1);
  worker_->NotifyLocalJobDone(job);
  return true;
}

bool WorkerManager::GetLocalJobDoneList(JobList* buffer) {
  // TODO(omidm): this function should never get called, remove it.
  assert(false);
  pthread_mutex_lock(&local_job_done_list_lock_);
  int len = local_job_done_list_.size();
  buffer->splice(buffer->end(), local_job_done_list_);
  pthread_mutex_unlock(&local_job_done_list_lock_);
  // TODO(quhang): when a job is done.
  if (len != 0) {
    worker_->StatEndJob(len);
  }
  return true;
}

Job* WorkerManager::NextComputationJobToRun(
    WorkerThreadComputation* worker_thread) {
  pthread_mutex_lock(&scheduling_critical_section_lock_);
  ongoing_parallelism_ -= worker_thread->used_parallelism;
  busy_worker_thread_computation_list_.remove(worker_thread);
  idle_worker_thread_computation_list_.push_back(worker_thread);
  worker_thread->job_assigned = false;
  pthread_mutex_unlock(&scheduling_critical_section_lock_);

  TriggerScheduling();

  pthread_mutex_lock(&scheduling_critical_section_lock_);
  while (!worker_thread->job_assigned) {
    pthread_cond_wait(&worker_thread->thread_can_start,
                      &scheduling_critical_section_lock_);
  }
  worker_thread->job_assigned = false;
  Job* temp = worker_thread->next_job_to_run;
  assert(temp != NULL);
  worker_thread->next_job_to_run = NULL;
  pthread_mutex_unlock(&scheduling_critical_section_lock_);

  return temp;
}

bool WorkerManager::SendCommand(SchedulerCommand* command) {
  assert(worker_ != NULL);
  worker_->SendCommand(command);
  return true;
}

bool WorkerManager::LaunchThread(WorkerThread* worker_thread) {
  int error_code =
      pthread_create(&worker_thread->thread_id, NULL,
                     ThreadEntryPoint, worker_thread);
  assert(error_code == 0);
  return true;
}

bool WorkerManager::StartWorkerThreads() {
  // LaunchThread(new WorkerThreadMonitor(this));
  pthread_mutex_lock(&scheduling_critical_section_lock_);
  for (int i = 0; i < computation_thread_num; ++i) {
    busy_worker_thread_computation_list_.push_back(
        new WorkerThreadComputation(this));
    LaunchThread(busy_worker_thread_computation_list_.back());
  }
  if (use_auxiliary_thread_) {
    worker_thread_auxiliary_ =
        new WorkerThreadAuxiliary(this, &task_thread_pool_);
    LaunchThread(worker_thread_auxiliary_);
    worker_thread_auxiliary_->SetThreadNum(4);
  }
  int error_code = pthread_create(
      &scheduling_id_, NULL, SchedulingEntryPoint, this);
  assert(error_code == 0);
  pthread_mutex_unlock(&scheduling_critical_section_lock_);
  return true;
}

void WorkerManager::ScheduleComputationJobs() {
  pthread_mutex_lock(&computation_job_queue_lock_);
  pthread_mutex_lock(&scheduling_critical_section_lock_);
  while (!idle_worker_thread_computation_list_.empty() &&
         !computation_job_list_.empty()) {
    WorkerThreadComputation* worker_thread =
        idle_worker_thread_computation_list_.front();

    idle_worker_thread_computation_list_.pop_front();
    busy_worker_thread_computation_list_.push_back(worker_thread);

    assert(!worker_thread->job_assigned);
    worker_thread->next_job_to_run =
        computation_job_list_.front();
    computation_job_list_.pop_front();
    worker_thread->job_assigned = true;

    task_thread_pool_.FreeTaskThreads(worker_thread->allocated_threads);
    worker_thread->allocated_threads.clear();
    if (inside_job_parallism >= 2 ||
        worker_thread->next_job_to_run->SupportMultiThread()) {
      bool status =
          task_thread_pool_.AllocateTaskThreads(
              inside_job_parallism, &worker_thread->allocated_threads);
      assert(status);
      ongoing_parallelism_ += inside_job_parallism;
      worker_thread->used_parallelism = inside_job_parallism;
    } else {
      ++ongoing_parallelism_;
      worker_thread->used_parallelism = 1;
    }
    ++dispatched_computation_job_count_;
    --ready_jobs_count_;
    pthread_cond_signal(&worker_thread->thread_can_start);
  }
  if (use_auxiliary_thread_) {
    int physical_core_for_auxiliary_jobs =
        PHYSICAL_CORE_NUM - (ongoing_parallelism_ + 1) / 2;
    if (physical_core_for_auxiliary_jobs < 1) {
      physical_core_for_auxiliary_jobs = 1;
    }
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    for (int i = 0; i < physical_core_for_auxiliary_jobs; ++i) {
      CPU_SET(LOGICAL_CORE_X[i], &cpuset);
      CPU_SET(LOGICAL_CORE_Y[i], &cpuset);
    }
    worker_thread_auxiliary_->SetThreadAffinity(&cpuset);
    worker_thread_auxiliary_->SetThreadNum(
        physical_core_for_auxiliary_jobs * 8);
  }
  pthread_mutex_unlock(&scheduling_critical_section_lock_);
  pthread_mutex_unlock(&computation_job_queue_lock_);
}

void* WorkerManager::ThreadEntryPoint(void* parameters) {
  {
    struct sched_param param;
    param.sched_priority = 0;
    int st = pthread_setschedparam(pthread_self(), SCHED_BATCH, &param);
    if (st != 0) {
      // Scheduling setting goes wrong.
      std::exit(1);
    }
  }
  WorkerThread* worker_thread = reinterpret_cast<WorkerThread*>(parameters);
  worker_thread->Run();
  assert(false);
}

void* WorkerManager::SchedulingEntryPoint(
    void* parameters) {
  WorkerManager* worker_manager = reinterpret_cast<WorkerManager*>(parameters);
  while (true) {
    pthread_mutex_lock(&worker_manager->scheduling_needed_lock_);
    while (!worker_manager->scheduling_needed_) {
      pthread_cond_wait(&worker_manager->scheduling_needed_cond_,
                        &worker_manager->scheduling_needed_lock_);
    }
    worker_manager->scheduling_needed_ = false;
    pthread_mutex_unlock(&worker_manager->scheduling_needed_lock_);

    // Thread scheduling algorithm is invoked when new job is added or new
    // threads turn to idle.
    worker_manager->ScheduleComputationJobs();
  }
  assert(false);
}

void WorkerManager::TriggerScheduling() {
  pthread_mutex_lock(&scheduling_needed_lock_);
  scheduling_needed_ = true;
  pthread_cond_signal(&scheduling_needed_cond_);
  pthread_mutex_unlock(&scheduling_needed_lock_);
}

}  // namespace nimbus
