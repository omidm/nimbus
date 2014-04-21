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
#include "worker/worker_thread_finish.h"

#include "worker/worker_manager.h"

namespace nimbus {

WorkerManager::WorkerManager() {
  log_ready_ = false;
  pthread_mutex_init(&computation_job_queue_lock_, NULL);
  pthread_mutex_init(&finish_job_queue_lock_, NULL);
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

WorkerManager::~WorkerManager() {
  pthread_mutex_destroy(&computation_job_queue_lock_);
  pthread_mutex_destroy(&finish_job_queue_lock_);
}

Job* WorkerManager::PullComputationJob() {
  pthread_mutex_lock(&computation_job_queue_lock_);
  if (computation_job_list_.empty()) {
    pthread_mutex_unlock(&computation_job_queue_lock_);
    return NULL;
  } else {
    Job* temp = computation_job_list_.front();
    computation_job_list_.pop_front();
    pthread_mutex_unlock(&computation_job_queue_lock_);
    return temp;
  }
}

bool WorkerManager::PushComputationJob(Job* job) {
  pthread_mutex_lock(&computation_job_queue_lock_);
  computation_job_list_.push_back(job);
  pthread_mutex_unlock(&computation_job_queue_lock_);
  return true;
}

Job* WorkerManager::PullFinishJob() {
  pthread_mutex_lock(&finish_job_queue_lock_);
  if (finish_job_list_.empty()) {
    pthread_mutex_unlock(&finish_job_queue_lock_);
    return NULL;
  } else {
    Job* temp = finish_job_list_.front();
    finish_job_list_.pop_front();
    pthread_mutex_unlock(&finish_job_queue_lock_);
    return temp;
  }
}

bool WorkerManager::PushFinishJob(Job* job) {
  pthread_mutex_lock(&finish_job_queue_lock_);
  finish_job_list_.push_back(job);
  pthread_mutex_unlock(&finish_job_queue_lock_);
  return true;
}

// TODO(quhang) sychronization effort needed?
bool WorkerManager::SendCommand(SchedulerCommand* command) {
  assert(worker_ != NULL);
  worker_->SendCommand(command);
  return true;
}

bool WorkerManager::StartWorkerThreads(int thread_number) {
  WorkerThread* worker_thread =
      new WorkerThreadFinish(this, worker_->data_map());
  worker_thread_list.push_back(worker_thread);
  assert(log_ready_);
  worker_thread->SetLoggingInterface(
      log_, version_log_, data_hash_log_, timer_);
  int error_code =
      pthread_create(&worker_thread->thread_id, NULL,
                     ThreadEntryPoint, worker_thread);
  assert(error_code == 0);

  for (int i = 0; i < thread_number; ++i) {
    WorkerThread* worker_thread = new WorkerThreadComputation(this);
    worker_thread_list.push_back(worker_thread);
    assert(log_ready_);
    worker_thread->SetLoggingInterface(
        log_, version_log_, data_hash_log_, timer_);
    int error_code =
        pthread_create(&worker_thread->thread_id, NULL,
                       ThreadEntryPoint, worker_thread);
    assert(error_code == 0);
  }
  return true;
}

void* WorkerManager::ThreadEntryPoint(void* parameters) {
  WorkerThread* worker_thread = reinterpret_cast<WorkerThread*>(parameters);
  worker_thread->Run();
  assert(false);
}

}  // namespace nimbus
