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

#include <sys/syscall.h>
#include <unistd.h>
#include <string>

#include "shared/profiler_malloc.h"
#include "worker/worker.h"
#include "worker/worker_manager.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_computation.h"

namespace nimbus {

WorkerThreadComputation::WorkerThreadComputation(WorkerManager* worker_manager)
    : WorkerThread(worker_manager) {
  use_threading_ = false;
  core_quota_ = 1;
  thread_queue = NULL;
}

WorkerThreadComputation::~WorkerThreadComputation() {
}

void WorkerThreadComputation::Run() {
  Job* job;
  while (true) {
    job = worker_manager_->NextComputationJobToRun(this);
    assert(job != NULL);
    job->set_use_threading(use_threading_);
    job->set_core_quota(core_quota_);
    thread_queue = NULL;
    job->set_thread_queue_hook(&thread_queue);
    ExecuteJob(job);
    assert(thread_queue == NULL);
    assert(worker_manager_ != NULL);
    bool success_flag = worker_manager_->FinishJob(job);
    assert(success_flag);
  }
}

void WorkerThreadComputation::ExecuteJob(Job* job) {
#ifndef MUTE_LOG
  log_->StartTimer();
  timer_->Start(job->id().elem());
#endif  // MUTE_LOG
#ifdef CACHE_LOG
  std::string jname = job->name();
  bool print_clog = false;
  bool print_cclog = false;
  if (jname.substr(0, 7) == "Compute") {
    print_clog = true;
  } else if (jname.find("Copy") != std::string::npos) {
    print_cclog = true;
  }
  if (print_clog) {
    std::stringstream msg;
    pid_t tid = syscall(SYS_gettid);
    msg << "~~~ TID: " << tid << " App compute job start : " << jname << " " << cache_log_->GetTime(); // NOLINT
    cache_log_->WriteToFile(msg.str());
  }
  if (print_cclog) {
    std::stringstream msg;
    pid_t tid = syscall(SYS_gettid);
    msg << "~~~ TID: " << tid << " App copy job start : " << jname << " " << cache_log_->GetTime();
    cache_log_->WriteToFile(msg.str());
  }
#endif
  // ProfilerMalloc::ResetThreadStatisticsByTid(pthread_self());
  dbg(DBG_WORKER, "[WORKER_THREAD] Execute job, name=%s, id=%lld. \n",
      job->name().c_str(), job->id().elem());
  job->Execute(job->parameters(), job->data_array);
  dbg(DBG_WORKER, "[WORKER_THREAD] Finish executing job, name=%s, id=%lld. \n",
      job->name().c_str(), job->id().elem());
  // size_t max_alloc = ProfilerMalloc::AllocMaxTid(pthread_self());
  // job->set_max_alloc(max_alloc);

#ifdef CACHE_LOG
  if (print_clog) {
    std::stringstream msg;
    pid_t tid = syscall(SYS_gettid);
    msg << "~~~ TID: " << tid << " App compute job end : " << jname << " " << cache_log_->GetTime();
    cache_log_->WriteToFile(msg.str());
  }
  if (print_cclog) {
    std::stringstream msg;
    pid_t tid = syscall(SYS_gettid);
    msg << "~~~ TID: " << tid << " App copy job end : " << jname << " " << cache_log_->GetTime();
    cache_log_->WriteToFile(msg.str());
  }
#endif
#ifndef MUTE_LOG
  double run_time = timer_->Stop(job->id().elem());
  log_->StopTimer();

  job->set_run_time(run_time);

  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff),
      "Execute Job, name: %35s  id: %6llu  length(s): %2.3lf  time(s): %6.3lf",
           job->name().c_str(), job->id().elem(),
           log_->timer(), log_->GetTime());
  log_->WriteToOutputStream(std::string(buff), LOG_INFO);

  char time_buff[LOG_MAX_BUFF_SIZE];
  snprintf(time_buff, sizeof(time_buff),
      "Queue Time: %2.9lf, Run Time: %2.9lf",
      job->wait_time(), job->run_time());
  log_->WriteToOutputStream(std::string(time_buff), LOG_INFO);
#endif  // MUTE_LOG
}

}  // namespace nimbus
