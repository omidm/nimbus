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
#include <sys/syscall.h>
#include <string>

#include "worker/worker.h"
#include "worker/worker_manager.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_auxiliary.h"
#include "worker/util_dumping.h"

namespace nimbus {

WorkerThreadAuxiliary::WorkerThreadAuxiliary(
    WorkerManager* worker_manager, TaskThreadPool* task_thread_pool)
    : WorkerThread(worker_manager) {
  task_thread_pool_ = task_thread_pool;
  active_thread_num_ = 0;
  thread_num_ = 0;
  pthread_mutex_init(&internal_lock_, NULL);
  pthread_cond_init(&internal_cond_, NULL);
}

WorkerThreadAuxiliary::~WorkerThreadAuxiliary() {
  // TODO(quhang): not implemented yet.
}

void WorkerThreadAuxiliary::Run() {
  struct sched_param param;
  param.sched_priority = 0;
  int st = pthread_setschedparam(pthread_self(), SCHED_OTHER, &param);
  if (st != 0) {
    // Scheduling setting goes wrong.
    exit(1);
  }
  while (true) {
    std::list<Job*> temp_job_list;
    bool success_flag =
        worker_manager_->PullAuxiliaryJobs(this, &temp_job_list);
    assert(success_flag);
    assert(worker_manager_ != NULL);
    pthread_mutex_lock(&internal_lock_);
    ready_job_list_.splice(ready_job_list_.end(), temp_job_list);
    pthread_cond_broadcast(&internal_cond_);
    pthread_mutex_unlock(&internal_lock_);
  }
}

void WorkerThreadAuxiliary::SetThreadNum(const int thread_num) {
  assert(thread_num >= 0);
  pthread_mutex_lock(&internal_lock_);
  if (thread_num > thread_num_) {
    for (int i = 0; i < thread_num - thread_num_; ++i) {
      WorkerTaskThreadAuxiliary* task_thread = NULL;
      if (stopped_task_threads_.empty()) {
        task_thread = new WorkerTaskThreadAuxiliary();
        task_thread->task_thread_wrapper =
            task_thread_pool_->AllocateTaskThread();
        assert(task_thread->task_thread_wrapper);
      } else {
        task_thread = stopped_task_threads_.front();
      }
      task_thread->worker_thread_auxiliary = this;
      task_thread->task_thread_wrapper->Run(
          WorkerTaskThreadAuxiliary::TaskThreadEntryPoint, this);
    }
  } else {
    for (int i = 0; i < thread_num_ - thread_num; ++i) {
      ready_job_list_.push_back(NULL);
    }
  }
  thread_num_ = thread_num;
  // Allocate more threads or deallocate threads.
  pthread_cond_broadcast(&internal_cond_);
  pthread_mutex_unlock(&internal_lock_);
}

void WorkerThreadAuxiliary::MainLoop(WorkerTaskThreadAuxiliary* task_thread) {
  pthread_mutex_lock(&internal_lock_);
  while (true) {
    while (ready_job_list_.empty()) {
      pthread_cond_wait(&internal_cond_, &internal_lock_);
    }
    Job* temp = ready_job_list_.front();
    ready_job_list_.pop_front();
    ++active_thread_num_;
    pthread_mutex_unlock(&internal_lock_);

    if (temp == NULL) {
      pthread_mutex_lock(&internal_lock_);
      --active_thread_num_;
      stopped_task_threads_.push_back(task_thread);
      pthread_mutex_unlock(&internal_lock_);
      throw ExceptionTaskThreadExit(NULL);
      // Never reached.
    }

    temp->set_worker_thread(this);
    ProcessJob(temp);
    temp->set_worker_thread(NULL);
    assert(worker_manager_ != NULL);
    bool success_flag = worker_manager_->FinishJob(temp);
    assert(success_flag);

    pthread_mutex_lock(&internal_lock_);
    --active_thread_num_;
  }
}


void WorkerThreadAuxiliary::ProcessJob(Job* job) {
#ifdef CACHE_LOG
  std::string jname = job->name();
  bool print_clog = false;
  if (cache_log_) {
    if (jname.find("Copy") != std::string::npos)
      print_clog = true;
    if (print_clog) {
      std::stringstream msg;
      pid_t tid = syscall(SYS_gettid);
      msg << "~~~ TID: " << tid << " App copy job start : " << jname << " " <<
        cache_log_->GetTime();
      cache_log_->WriteToFile(msg.str());
    }
  }
#endif
  dbg(DBG_WORKER, "[WORKER_THREAD] Execute auxiliary job, name=%s, id=%lld. \n",
      job->name().c_str(), job->id().elem());
  job->Execute(job->parameters(), job->data_array);
  dbg(DBG_WORKER, "[WORKER_THREAD] Finish executing auxiliary job, "
      "name=%s, id=%lld. \n", job->name().c_str(), job->id().elem());
#ifdef CACHE_LOG
  if (cache_log_) {
    if (print_clog) {
      std::stringstream msg;
      pid_t tid = syscall(SYS_gettid);
      msg << "~~~ TID: " << tid << " App copy job end : " << jname << " " <<
        cache_log_->GetTime();
      cache_log_->WriteToFile(msg.str());
    }
  }
#endif
}

void* WorkerTaskThreadAuxiliary::TaskThreadEntryPoint(void* parameter) {
  WorkerTaskThreadAuxiliary* stub =
      reinterpret_cast<WorkerTaskThreadAuxiliary*>(parameter);
  stub->worker_thread_auxiliary->MainLoop(stub);
  return NULL;
}

}  // namespace nimbus
