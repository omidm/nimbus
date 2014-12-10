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
class WorkerThreadAuxiliary;
class Worker;

class WorkerManager {
  friend class WorkerThreadMonitor;
 public:
  // Configuration: used parallism level.
  static int inside_job_parallism;
  static int across_job_parallism;
  // Configuration: the number of computation threads used.
  int computation_thread_num;

  // Constructor and deconstructor
  explicit WorkerManager();
  ~WorkerManager();
  // Starts all the worker threads and the scheduling thread.
  bool StartWorkerThreads();

  // Interfaces for worker.
  // 1. Push a job.
  bool PushJob(Job* job);
  // 2. Push a list of jobs.
  bool PushJobList(std::list<Job*>* job_list);
  // 3. Retrieve the job done list.
  bool GetLocalJobDoneList(JobList* buffer);
  // 4. Trigger scheduling.
  void TriggerScheduling();

  // Interfaces for worker threads.
  // 1. Retrieve next computation job.
  Job* NextComputationJobToRun(WorkerThreadComputation* worker_thread);
  // 2. Send a command.
  bool SendCommand(SchedulerCommand* command);
  // 3. Finish a job
  bool FinishJob(Job* job);
  // 4. Retrieve next auxiliary job list.
  bool PullAuxiliaryJobs(WorkerThreadAuxiliary* worker_thread,
                         std::list<Job*>* job_list);

 public:
  Worker* worker_;

 private:
  // True when a scheduling algorithms is needed to be triggered.
  bool ClassifyAndAddJob(Job* job);
  // Thread scheduling algorithm.
  void ScheduleComputationJobs();
  // The thread pool.
  TaskThreadPool task_thread_pool_;

  // Threads.
  // 1. The thread to run scheduling.
  pthread_t scheduling_id_;
  static void* SchedulingEntryPoint(void* parameters);
  // 2. The worker threads.
  bool LaunchThread(WorkerThread* worker_thread);
  // Entry point for each worker thread.
  static void* ThreadEntryPoint(void* parameters);

  // The scheduling triggering mechanism.
  pthread_mutex_t scheduling_needed_lock_;
  pthread_cond_t scheduling_needed_cond_;
  bool scheduling_needed_;

  // Data structures and locks used for scheduling.
  // 1. Lock to protect schedulign critical session.
  pthread_mutex_t scheduling_critical_section_lock_;
  // 2. The list of ready computation jobs.
  pthread_mutex_t computation_job_queue_lock_;
  std::list<Job*> computation_job_list_;
  // 3. The list of ready computation jobs.
  pthread_mutex_t auxiliary_job_queue_lock_;
  std::list<Job*> auxiliary_job_list_;
  // 4. The list of finished jobs.
  pthread_mutex_t local_job_done_list_lock_;
  JobList local_job_done_list_;

 public:
  // Logging: set up logging interface.
  void SetLoggingInterface(
      Log* log, Log* version_log, Log* data_hash_log, Log* cache_log,
      HighResolutionTimer* timer);

  void SetEventLog(std::string wid_str);

 private:
  // Internal logging facility.
  void PrintTimeStamp(const char* event, const char* s, const uint64_t d);
  FILE* event_log;
  bool log_ready_;
  Log* log_;
  Log* version_log_;
  Log* data_hash_log_;
  Log* cache_log_;
  HighResolutionTimer* timer_;
  // Performance measurement.
  int64_t dispatched_computation_job_count_;
  int64_t ready_jobs_count_;

  static const bool use_auxiliary_thread_ = false;
  WorkerThreadAuxiliary* worker_thread_auxiliary_;
  std::list<WorkerThreadComputation*> idle_worker_thread_computation_list_;
  std::list<WorkerThreadComputation*> busy_worker_thread_computation_list_;
  int32_t ongoing_parallelism_;
};
}  // namespace nimbus

#endif  // NIMBUS_WORKER_WORKER_MANAGER_H_
