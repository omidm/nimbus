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
  * The worker manager extracts jobs from the ready queue, and dispatches jobs
  * to different worker threads.
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_WORKER_WORKER_MANAGER_H_
#define NIMBUS_WORKER_WORKER_MANAGER_H_

#include <cstdio>
#include <list>
#include <string>
#include "shared/high_resolution_timer.h"
#include "shared/log.h"
#include "shared/nimbus.h"
#include "worker/task_thread_pool.h"

namespace nimbus {
class SchedulerCommand;
class WorkerThreadMonitor;
class WorkerThread;
class WorkerThreadComputation;
class Worker;

class WorkerManager {
  friend class WorkerThreadMonitor;
 public:
  static int inside_job_parallism;
  static int across_job_parallism;
  explicit WorkerManager();
  ~WorkerManager();

  void SetLoggingInterface(
      Log* log, Log* version_log, Log* data_hash_log, Log* cache_log,
      HighResolutionTimer* timer);

  // Computation job: computation-intensive job.
  // Fast job: non-computation-intensive job.

  Job* NextComputationJobToRun(WorkerThreadComputation* worker_thread);

  bool PushJob(Job* job);
  bool PushJobList(std::list<Job*>* job_list);
  bool FinishJob(Job* job);

  bool SendCommand(SchedulerCommand* command);

  bool GetLocalJobDoneList(JobList* buffer);

  // Starts all the worker threads and the scheduling thread.
  bool StartWorkerThreads();

  // Configuration for the number of threads used.
  int computation_thread_num;

 public:
  // TODO(quhang) Not sure if maintaining such a pointer is good or not.
  Worker* worker_;
  /*
  static void RegisterThread(const pthread_t* child_thread,
                             const pthread_t* parenter_thread = NULL);
  struct ThreadStateStub {
    size_t allocated_memory;
    size_t freed_memory;
    bool is_active;
    ThreadStateStub() {
      allocated_memory = 0;
      freed_memory = 0;
      is_active = false;
    }
  };
  static void RegisterThreadStateStub(
      const pthread_t* thread, const ThreadStateStub* thread_stub);
  static void DeregisterThreadStateStub(const pthread_t* thread);
  static pthread_key_t* thread_state_stub_key;
  static void AllocateThreadStateStub();
  static ThreadStateStub* DetachThreadStateStub();
  */

 private:
  // Thread scheduling algorithm.
  void ScheduleComputationJobs();

  pthread_mutex_t scheduling_needed_lock_;
  pthread_cond_t scheduling_needed_cond_;
  // Protected by scheduling_needed_lock_.
  bool scheduling_needed_;
  // Triggers the scheduling algorithm.
  void TriggerScheduling();

  std::list<WorkerThread*> worker_thread_list_;
  bool LaunchThread(WorkerThread* worker_thread);
  // Entry point for each worker thread.
  static void* ThreadEntryPoint(void* parameters);

  pthread_t scheduling_id_;
  // Entry point for the scheduling thread.
  static void* SchedulingEntryPoint(void* parameters);

  pthread_mutex_t scheduling_critical_section_lock_;

  pthread_mutex_t computation_job_queue_lock_;
  // Protected by computation_job_queue_lock.
  std::list<Job*> computation_job_list_;

  std::list<Job*> computation_job_to_schedule_list_;

  pthread_mutex_t local_job_done_list_lock_;
  // Protected by local_job_done_queue_lock_.
  JobList local_job_done_list_;

  // Measures running states of the worker.
  int64_t dispatched_computation_job_count_;
  int idle_computation_threads_;
  int64_t ready_jobs_count_;

  // Logging data structures.
  void PrintTimeStamp(const char* event, const char* s, const uint64_t d);
  FILE* event_log;
  bool log_ready_;
  Log* log_;
  Log* version_log_;
  Log* data_hash_log_;
  Log* cache_log_;
  HighResolutionTimer* timer_;

  TaskThreadPool task_thread_pool_;
};
}  // namespace nimbus

#endif  // NIMBUS_WORKER_WORKER_MANAGER_H_
