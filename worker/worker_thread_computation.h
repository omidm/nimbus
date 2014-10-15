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
  * A worker thread that executes computation jobs that is dispatched by the
  * worker manager.
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_WORKER_WORKER_THREAD_COMPUTATION_H_
#define NIMBUS_WORKER_WORKER_THREAD_COMPUTATION_H_

#include <string>
#include "shared/nimbus.h"
#include "worker/task_thread_pool.h"

namespace nimbus {
class WorkerThreadComputation : public WorkerThread {
 public:
  explicit WorkerThreadComputation(WorkerManager* worker_manager);
  virtual ~WorkerThreadComputation();
  virtual void Run();
  void set_core_quota(int core_quota) {
    core_quota_ = core_quota;
  }
  int core_quota() {
    return core_quota_;
  }
  void set_use_threading(bool use_threading) {
    use_threading_ = use_threading;
  }
  bool use_threading() {
    return use_threading_;
  }
  TaskThreadPool::TaskThreadList allocated_threads;

 private:
  int core_quota_;
  bool use_threading_;
  void ExecuteJob(Job* job);
  uint64_t ParseLine(std::string line);
};
}  // namespace nimbus

#endif  // NIMBUS_WORKER_WORKER_THREAD_COMPUTATION_H_
