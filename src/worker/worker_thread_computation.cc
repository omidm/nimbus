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

#include "src/shared/profiler_malloc.h"
#include "src/shared/fast_log.h"
#include "src/worker/worker.h"
#include "src/worker/worker_manager.h"
#include "src/worker/worker_thread.h"
#include "src/worker/worker_thread_computation.h"

namespace nimbus {

WorkerThreadComputation::WorkerThreadComputation(WorkerManager* worker_manager)
    : WorkerThread(worker_manager) {
  pthread_cond_init(&thread_can_start, NULL);
  // next_job_to_run = NULL;
  // job_assigned = false;
  job_list_assigned_ = false;
  used_parallelism = 0;
}

WorkerThreadComputation::~WorkerThreadComputation() {
  pthread_cond_destroy(&thread_can_start);
}

void WorkerThreadComputation::Run() {
  timer::InitializeTimers();
  timer::StartTimer(timer::kTotal);
  Job* job;
  while (true) {
    // job = worker_manager_->NextComputationJobToRun(this);
    // assert(job != NULL);
    // job->set_worker_thread(this);
    // ExecuteJob(job);
    // job->set_worker_thread(NULL);
    // assert(worker_manager_ != NULL);
    // bool success_flag = worker_manager_->FinishJob(job);
    // assert(success_flag);
    worker_manager_->NextComputationJobToRun(this);
    assert(job_list_to_run_.size() > 0);
    JobList::iterator iter = job_list_to_run_.begin();
    for (; iter != job_list_to_run_.end(); ++iter) {
      Job *job = *iter;
      job->set_worker_thread(this);
      ExecuteJob(job);
      job->set_worker_thread(NULL);
      assert(worker_manager_ != NULL);
      bool success_flag = worker_manager_->FinishJob(job);
      assert(success_flag);
    }
    job_list_to_run_.clear();
  }
}

void WorkerThreadComputation::ExecuteJob(Job* job) {
  ProfilerMalloc::ResetBaseAlloc();
  dbg(DBG_WORKER, "[WORKER_THREAD] Execute job, name=%s, id=%lld. \n",
      job->name().c_str(), job->id().elem());
  if (dynamic_cast<RemoteCopySendJob*>(job) ||  // NOLINT
      dynamic_cast<RemoteCopyReceiveJob*>(job) ||  // NOLINT
      dynamic_cast<MegaRCRJob*>(job) ||  // NOLINT
      dynamic_cast<LocalCopyJob*>(job) ||  // NOLINT
      dynamic_cast<CreateDataJob*>(job)) {  // NOLINT
    timer::StartTimer(timer::kExecuteCopyJob);
    job->Execute(job->parameters(), job->data_array);
    timer::StopTimer(timer::kExecuteCopyJob);
  } else {
    timer::StartTimer(timer::kExecuteComputationJob);
    bool parent = !job->sterile();
    if (parent) {
      timer::StartTimer(timer::kExecuteParentJob);
    }
    job->Execute(job->parameters(), job->data_array);
    if (parent) {
      timer::StopTimer(timer::kExecuteParentJob);
    }
    timer::StopTimer(timer::kExecuteComputationJob);
  }
  dbg(DBG_WORKER, "[WORKER_THREAD] Finish executing job, name=%s, id=%lld. \n",
      job->name().c_str(), job->id().elem());
  // size_t max_alloc = ProfilerMalloc::AllocMaxTid(pthread_self());
  // job->set_max_alloc(max_alloc);
}
}  // namespace nimbus
