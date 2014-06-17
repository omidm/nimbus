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

#include <string>

#include "worker/worker.h"
#include "worker/worker_manager.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_fast.h"
#include "worker/util_dumping.h"

namespace nimbus {

WorkerThreadFast::WorkerThreadFast(WorkerManager* worker_manager)
    : WorkerThread(worker_manager) {
}

WorkerThreadFast::~WorkerThreadFast() {
}

void WorkerThreadFast::Run() {
  std::list<Job*> fast_job_list;
  while (true) {
    bool success_flag = worker_manager_->PullFastJobs(this, &fast_job_list);
    assert(success_flag);
    assert(worker_manager_ != NULL);
    for (std::list<Job*>::iterator index = fast_job_list.begin();
         index != fast_job_list.end();
         index++) {
      ProcessJob(*index);
      bool success_flag = worker_manager_->FinishJob(*index);
      assert(success_flag);
    }
    fast_job_list.clear();
  }
}

void WorkerThreadFast::ProcessJob(Job* job) {
#ifdef CACHE_LOG
  std::string jname = job->name();
  bool print_clog = false;
  if (jname.find("Copy") != std::string::npos)
    print_clog = true;
  if (print_clog) {
    std::stringstream msg;
    msg << "~~~ App copy job start : " << jname << " " << cache_log_->GetTime();
    cache_log_->WriteToFile(msg.str());
  }
#endif
  dbg(DBG_WORKER, "[WORKER_THREAD] Execute fast job, name=%s, id=%lld. \n",
      job->name().c_str(), job->id().elem());
  job->Execute(job->parameters(), job->data_array);
  dbg(DBG_WORKER, "[WORKER_THREAD] Finish executing fast job, "
      "name=%s, id=%lld. \n", job->name().c_str(), job->id().elem());
#ifdef CACHE_LOG
  if (print_clog) {
    std::stringstream msg;
    msg << "~~~ App copy job end : " << jname << " " << cache_log_->GetTime();
    cache_log_->WriteToFile(msg.str());
  }
#endif
}

}  // namespace nimbus
