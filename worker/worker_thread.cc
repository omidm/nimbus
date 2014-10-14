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

#include <list>

#include "worker/worker_manager.h"
#include "worker/worker_thread.h"

namespace nimbus {

WorkerThread::WorkerThread(WorkerManager* worker_manager) {
  worker_manager_ = worker_manager;
  pthread_cond_init(&thread_can_start, NULL);
  next_job_to_run = NULL;
  idle = true;
  job_assigned = false;
  pthread_mutex_init(&internal_thread_list_lock_, NULL);
  ClearRegisterThreads();
}

WorkerThread::~WorkerThread() {
  pthread_cond_destroy(&thread_can_start);
  pthread_mutex_destroy(&internal_thread_list_lock_);
}

void WorkerThread::SetLoggingInterface(
    Log* log, Log* version_log, Log* data_hash_log, Log* cache_log,
    HighResolutionTimer* timer) {
  log_ = log;
  version_log_ = version_log;
  data_hash_log_ = data_hash_log;
  cache_log_ = cache_log;
  timer_ = timer;
}

bool WorkerThread::RegisterThread(const pthread_t& thread_handle) {
  pthread_mutex_lock(&internal_thread_list_lock_);
  internal_thread_list_.push_back(thread_handle);
  pthread_mutex_unlock(&internal_thread_list_lock_);
  return false;
}

bool WorkerThread::DeregisterThread(const pthread_t& thread_handle) {
  pthread_mutex_lock(&internal_thread_list_lock_);
  for (std::list<pthread_t>::iterator iter = internal_thread_list_.begin();
       iter != internal_thread_list_.end();
       ++iter) {
    if (pthread_equal(*iter, thread_handle) != 0) {
      internal_thread_list_.erase(iter);
      pthread_mutex_unlock(&internal_thread_list_lock_);
      return true;
    }
  }
  pthread_mutex_unlock(&internal_thread_list_lock_);
  return false;
}

void WorkerThread::ClearRegisterThreads() {
  pthread_mutex_lock(&internal_thread_list_lock_);
  internal_thread_list_.clear();
  pthread_mutex_unlock(&internal_thread_list_lock_);
}

void WorkerThread::SetThreadAffinity(
    const size_t cpusetsize, const cpu_set_t* cpuset) {
  pthread_mutex_lock(&internal_thread_list_lock_);
  for (std::list<pthread_t>::iterator iter = internal_thread_list_.begin();
       iter != internal_thread_list_.end();
       ++iter) {
    pthread_setaffinity_np(*iter, cpusetsize, cpuset);
  }
  pthread_mutex_unlock(&internal_thread_list_lock_);
}

}  // namespace nimbus
