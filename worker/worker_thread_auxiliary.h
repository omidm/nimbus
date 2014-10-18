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

#ifndef NIMBUS_WORKER_WORKER_THREAD_AUXILIARY_H_
#define NIMBUS_WORKER_WORKER_THREAD_AUXILIARY_H_

#include <list>
#include "shared/nimbus.h"
#include "worker/task_thread_pool.h"

namespace nimbus {

class WorkerThreadAuxiliary;

class WorkerTaskThreadAuxiliary {
 public:
  WorkerTaskThreadAuxiliary() {
    task_thread_wrapper = NULL;
    worker_thread_auxiliary = NULL;
  }
  ~WorkerTaskThreadAuxiliary() {}
  TaskThreadWrapper* task_thread_wrapper;
  WorkerThreadAuxiliary* worker_thread_auxiliary;
  static void* TaskThreadEntryPoint(void* parameter);
};

class WorkerThreadAuxiliary : public WorkerThread {
 public:
  explicit WorkerThreadAuxiliary(
      WorkerManager* worker_manager, TaskThreadPool* task_thread_pool);
  virtual ~WorkerThreadAuxiliary();
  virtual void Run();
  void MainLoop(WorkerTaskThreadAuxiliary* task_thread);
  void SetThreadNum(const int thread_num);
 private:
  pthread_mutex_t internal_lock_;
  pthread_cond_t internal_cond_;
  Job* PullAuxiliaryJob();
  void ProcessJob(Job* job);
  TaskThreadPool* task_thread_pool_;
  int thread_num_;
  int active_thread_num_;
  std::list<Job*> ready_job_list_;
  std::list<WorkerTaskThreadAuxiliary*> stopped_task_threads_;
};
}  // namespace nimbus

#endif  // NIMBUS_WORKER_WORKER_THREAD_AUXILIARY_H_
