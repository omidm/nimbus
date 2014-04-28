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

#include "shared/nimbus.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_computation.h"
#include "worker/worker_thread_fast.h"
#include "worker/worker_thread_finish.h"
#include "worker/worker_thread_monitor.h"

#include "worker/worker_manager.h"

namespace nimbus {

WorkerManager::WorkerManager(bool multi_threaded) {
  pthread_mutex_init(&scheduling_needed_lock_, NULL);
  pthread_cond_init(&scheduling_needed_cond_, NULL);
  scheduling_needed_ = false;

  pthread_mutex_init(&scheduling_critical_section_lock_, NULL);

  log_ready_ = false;

  pthread_mutex_init(&computation_job_queue_lock_, NULL);

  pthread_mutex_init(&finish_job_queue_lock_, NULL);
  pthread_cond_init(&finish_job_queue_any_cond_, NULL);

  pthread_mutex_init(&fast_job_queue_lock_, NULL);
  pthread_cond_init(&fast_job_queue_any_cond_, NULL);

  if (multi_threaded) {
    computation_thread_num = 10;
    fast_thread_num = 1;
  } else {
    computation_thread_num = 1;
    fast_thread_num = 0;
  }
  idle_computation_threads_ = 0;
  dispatched_computation_job_count_ = 0;
  dispatched_finish_job_count_ = 0;
  dispatched_fast_job_count_= 0;
  ready_jobs_count_ = 0;
}

WorkerManager::~WorkerManager() {
  pthread_mutex_destroy(&scheduling_needed_lock_);
  pthread_cond_destroy(&scheduling_needed_cond_);

  pthread_mutex_destroy(&scheduling_critical_section_lock_);

  pthread_mutex_destroy(&computation_job_queue_lock_);

  pthread_mutex_destroy(&finish_job_queue_lock_);
  pthread_cond_destroy(&finish_job_queue_any_cond_);

  pthread_mutex_destroy(&fast_job_queue_lock_);
  pthread_cond_destroy(&fast_job_queue_any_cond_);
}

void WorkerManager::SetLoggingInterface(
    Log* log, Log* version_log, Log* data_hash_log,
    HighResolutionTimer* timer) {
  log_ = log;
  version_log_ = version_log;
  data_hash_log_ = data_hash_log_;
  timer_ = timer;
  log_ready_ = true;
}

bool WorkerManager::PushJob(Job* job) {
  if ((fast_thread_num != 0) && (
      (dynamic_cast<CreateDataJob*>(job) || // NOLINT
      dynamic_cast<LocalCopyJob*>(job) || // NOLINT
      dynamic_cast<RemoteCopySendJob*>(job) || // NOLINT
      dynamic_cast<RemoteCopyReceiveJob*>(job)))) { // NOLINT
    pthread_mutex_lock(&fast_job_queue_lock_);
    fast_job_list_.push_back(job);
    pthread_cond_signal(&fast_job_queue_any_cond_);
    pthread_mutex_unlock(&fast_job_queue_lock_);
  } else {
    pthread_mutex_lock(&computation_job_queue_lock_);
    computation_job_list_.push_back(job);
    ++ready_jobs_count_;
    pthread_mutex_unlock(&computation_job_queue_lock_);
    TriggerScheduling();
  }
  return true;
}

bool WorkerManager::PushFinishJob(Job* job) {
  pthread_mutex_lock(&finish_job_queue_lock_);
  finish_job_list_.push_back(job);
  pthread_cond_signal(&finish_job_queue_any_cond_);
  pthread_mutex_unlock(&finish_job_queue_lock_);
  return true;
}

Job* WorkerManager::NextComputationJobToRun(WorkerThread* worker_thread) {
  pthread_mutex_lock(&scheduling_critical_section_lock_);
  worker_thread->idle = true;
  worker_thread->job_assigned = false;
  ++idle_computation_threads_;
  pthread_mutex_unlock(&scheduling_critical_section_lock_);

  TriggerScheduling();

  pthread_mutex_lock(&scheduling_critical_section_lock_);
  while (!worker_thread->job_assigned) {
    pthread_cond_wait(&worker_thread->thread_can_start,
                      &scheduling_critical_section_lock_);
  }
  worker_thread->idle = false;
  --idle_computation_threads_;
  worker_thread->job_assigned = false;
  pthread_mutex_unlock(&scheduling_critical_section_lock_);

  Job* temp = worker_thread->next_job_to_run;
  assert(temp != NULL);
  worker_thread->next_job_to_run = NULL;

  return temp;
}

bool WorkerManager::PullFinishJobs(WorkerThread* worker_thread,
                                  std::list<Job*>* list_buffer) {
  pthread_mutex_lock(&finish_job_queue_lock_);
  worker_thread->idle = true;
  while (finish_job_list_.empty()) {
    pthread_cond_wait(&finish_job_queue_any_cond_,
                      &finish_job_queue_lock_);
  }
  list_buffer->clear();
  list_buffer->swap(finish_job_list_);
  dispatched_finish_job_count_ += list_buffer->size();
  worker_thread->idle = false;
  pthread_mutex_unlock(&finish_job_queue_lock_);
  return true;
}

bool WorkerManager::PullFastJobs(WorkerThread* worker_thread,
                                 std::list<Job*>* list_buffer) {
  pthread_mutex_lock(&fast_job_queue_lock_);
  worker_thread->idle = true;
  while (fast_job_list_.empty()) {
    pthread_cond_wait(&fast_job_queue_any_cond_,
                      &fast_job_queue_lock_);
  }
  list_buffer->clear();
  // TODO(quhang) hard-code.
  int count = 10;
  while (!fast_job_list_.empty() && count > 0) {
    list_buffer->push_back(fast_job_list_.front());
    fast_job_list_.pop_front();
    --count;
    ++dispatched_fast_job_count_;
  }
  worker_thread->idle = false;
  pthread_mutex_unlock(&fast_job_queue_lock_);
  return true;
}

bool WorkerManager::SendCommand(SchedulerCommand* command) {
  assert(worker_ != NULL);
  worker_->SendCommand(command);
  return true;
}

bool WorkerManager::LaunchThread(WorkerThread* worker_thread) {
  worker_thread_list_.push_back(worker_thread);
  assert(log_ready_);
  worker_thread->SetLoggingInterface(
      log_, version_log_, data_hash_log_, timer_);
  int error_code =
      pthread_create(&worker_thread->thread_id, NULL,
                     ThreadEntryPoint, worker_thread);
  assert(error_code == 0);
  return true;
}

bool WorkerManager::StartWorkerThreads() {
  LaunchThread(new WorkerThreadMonitor(this));
  LaunchThread(new WorkerThreadFinish(this, worker_->data_map()));
  for (int i = 0; i < fast_thread_num; ++i) {
    LaunchThread(new WorkerThreadFast(this));
  }
  for (int i = 0; i < computation_thread_num; ++i) {
    LaunchThread(new WorkerThreadComputation(this));
  }
  int error_code = pthread_create(
      &scheduling_id_, NULL, SchedulingEntryPoint, this);
  assert(error_code == 0);
  return true;
}

void WorkerManager::ScheduleComputationJobs() {
  // Try not to block the worker core thread when scheduling algorithm is
  // running.
  pthread_mutex_lock(&computation_job_queue_lock_);
  computation_job_to_schedule_list_.insert(
      computation_job_to_schedule_list_.end(),
      computation_job_list_.begin(),
      computation_job_list_.end());
  computation_job_list_.clear();
  pthread_mutex_unlock(&computation_job_queue_lock_);

  pthread_mutex_lock(&scheduling_critical_section_lock_);
  for (std::list<WorkerThread*>::iterator index = worker_thread_list_.begin();
       index != worker_thread_list_.end();
       ++index) {
    if (computation_job_to_schedule_list_.empty()) {
      break;
    }
    // TODO(quhang) RTTI is not good.
    WorkerThreadComputation* worker_thread =
        dynamic_cast<WorkerThreadComputation*>(*index);  // NOLINT
    if (worker_thread != NULL
        && worker_thread->idle && !worker_thread->job_assigned) {
      worker_thread->next_job_to_run =
          computation_job_to_schedule_list_.front();
      computation_job_to_schedule_list_.pop_front();
      worker_thread->job_assigned = true;
      ++dispatched_computation_job_count_;
      --ready_jobs_count_;
      pthread_cond_signal(&worker_thread->thread_can_start);
    }
  }
  pthread_mutex_unlock(&scheduling_critical_section_lock_);
}

void* WorkerManager::ThreadEntryPoint(void* parameters) {
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
