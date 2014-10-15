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
#include <map>
#include <string>
#include <vector>

#include "shared/nimbus.h"
#include "shared/profiler_malloc.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_computation.h"
#include "worker/worker_thread_monitor.h"

#include "worker/worker_manager.h"

namespace nimbus {

int WorkerManager::inside_job_parallism = 0;
int WorkerManager::across_job_parallism = 0;

WorkerManager::WorkerManager() {
  /*
  thread_state_stub_key = new pthread_key_t;
  int success_flag = pthread_key_create(thread_state_stub_key, NULL);
  assert(success_flag == 0);
  */
  event_log = fopen("event_be.txt", "w");
  pthread_mutex_init(&scheduling_needed_lock_, NULL);
  pthread_cond_init(&scheduling_needed_cond_, NULL);
  scheduling_needed_ = false;

  pthread_mutex_init(&scheduling_critical_section_lock_, NULL);

  log_ready_ = false;

  pthread_mutex_init(&computation_job_queue_lock_, NULL);

  pthread_mutex_init(&local_job_done_list_lock_, NULL);

  if (across_job_parallism <= 0) {
    computation_thread_num = 1;
  } else {
    computation_thread_num = across_job_parallism;
  }
  idle_computation_threads_ = 0;
  dispatched_computation_job_count_ = 0;
  ready_jobs_count_ = 0;
}

WorkerManager::~WorkerManager() {
  pthread_mutex_destroy(&scheduling_needed_lock_);
  pthread_cond_destroy(&scheduling_needed_cond_);

  pthread_mutex_destroy(&scheduling_critical_section_lock_);

  pthread_mutex_destroy(&computation_job_queue_lock_);
}

void WorkerManager::PrintTimeStamp(
    const char* event, const char* s, const uint64_t d) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(event_log, "%f %s %s %lu\n", time_sum, event, s, d);
}

void WorkerManager::SetLoggingInterface(
    Log* log, Log* version_log, Log* data_hash_log, Log* cache_log,
    HighResolutionTimer* timer) {
  log_ = log;
  version_log_ = version_log;
  data_hash_log_ = data_hash_log;
  cache_log_ = cache_log;
  timer_ = timer;
  log_ready_ = true;
}

bool WorkerManager::PushJob(Job* job) {
  pthread_mutex_lock(&computation_job_queue_lock_);
  computation_job_list_.push_back(job);
  ++ready_jobs_count_;
  pthread_mutex_unlock(&computation_job_queue_lock_);
  TriggerScheduling();
  return true;
}

bool WorkerManager::PushJobList(std::list<Job*>* job_list) {
  pthread_mutex_lock(&computation_job_queue_lock_);
  ready_jobs_count_ += job_list->size();
  computation_job_list_.splice(computation_job_list_.end(), *job_list);
  pthread_mutex_unlock(&computation_job_queue_lock_);
  TriggerScheduling();
  return true;
}

bool WorkerManager::FinishJob(Job* job) {
  pthread_mutex_lock(&local_job_done_list_lock_);
  PrintTimeStamp("f",
                 job->name().c_str(),
                 job->id().elem());
  local_job_done_list_.push_back(job);
  pthread_mutex_unlock(&local_job_done_list_lock_);
  return true;
}

bool WorkerManager::GetLocalJobDoneList(JobList* buffer) {
  pthread_mutex_lock(&local_job_done_list_lock_);
  buffer->splice(buffer->end(), local_job_done_list_);
  pthread_mutex_unlock(&local_job_done_list_lock_);
  return true;
}

Job* WorkerManager::NextComputationJobToRun(WorkerThread* worker_thread) {
  pthread_mutex_lock(&scheduling_critical_section_lock_);
  worker_thread->idle = true;
  worker_thread->job_assigned = false;
  ++idle_computation_threads_;
  pthread_mutex_unlock(&scheduling_critical_section_lock_);

  TriggerScheduling();

  pthread_mutex_lock(&scheduling_critical_section_lock_);
  while (!worker_thread->job_assigned) {
    pthread_cond_wait(&worker_thread->thread_can_start,
                      &scheduling_critical_section_lock_);
  }
  worker_thread->idle = false;
  --idle_computation_threads_;
  worker_thread->job_assigned = false;
  Job* temp = worker_thread->next_job_to_run;
  assert(temp != NULL);
  PrintTimeStamp("r",
                 worker_thread->next_job_to_run->name().c_str(),
                 worker_thread->next_job_to_run->id().elem());
  worker_thread->next_job_to_run = NULL;
  pthread_mutex_unlock(&scheduling_critical_section_lock_);

  return temp;
}

bool WorkerManager::SendCommand(SchedulerCommand* command) {
  assert(worker_ != NULL);
  worker_->SendCommand(command);
  return true;
}

bool WorkerManager::LaunchThread(WorkerThread* worker_thread) {
  worker_thread_list_.push_back(worker_thread);
  assert(log_ready_);
  worker_thread->SetLoggingInterface(
      log_, version_log_, data_hash_log_, cache_log_, timer_);
  int error_code =
      pthread_create(&worker_thread->thread_id, NULL,
                     ThreadEntryPoint, worker_thread);
  // WorkerManager::RegisterThread(&worker_thread->thread_id);
  assert(error_code == 0);
  return true;
}

bool WorkerManager::StartWorkerThreads() {
  // LaunchThread(new WorkerThreadMonitor(this));
  for (int i = 0; i < computation_thread_num; ++i) {
    LaunchThread(new WorkerThreadComputation(this));
  }
  int error_code = pthread_create(
      &scheduling_id_, NULL, SchedulingEntryPoint, this);
  assert(error_code == 0);
  return true;
}

void WorkerManager::ScheduleComputationJobs() {
  pthread_mutex_lock(&computation_job_queue_lock_);
  pthread_mutex_lock(&scheduling_critical_section_lock_);
  for (std::list<WorkerThread*>::iterator index = worker_thread_list_.begin();
       index != worker_thread_list_.end();
       ++index) {
    if (computation_job_list_.empty()) {
      break;
    }
    // TODO(quhang) RTTI is not good.
    WorkerThreadComputation* worker_thread =
        dynamic_cast<WorkerThreadComputation*>(*index);  // NOLINT
    if (worker_thread != NULL
        && worker_thread->idle && !worker_thread->job_assigned) {
      worker_thread->next_job_to_run =
          computation_job_list_.front();
      computation_job_list_.pop_front();
      // PrintTimeStamp("i",
      //                worker_thread->next_job_to_run->name().c_str(),
      //                worker_thread->next_job_to_run->id().elem());
      // dbg(DBG_WORKER_BD, DBG_WORKER_BD_S"Job(name %s, #%d) dispatched.\n",
      //     worker_thread->next_job_to_run->name().c_str(),
      //     worker_thread->next_job_to_run->id().elem());
      worker_thread->job_assigned = true;
      task_thread_pool_.FreeTaskThreads(worker_thread->allocated_threads);
      worker_thread->allocated_threads.clear();
      if (inside_job_parallism >= 2 ||
          worker_thread->next_job_to_run->SupportMultiThread()) {
        bool status =
            task_thread_pool_.AllocateTaskThreads(
                inside_job_parallism, &worker_thread->allocated_threads);
        assert(status);
      }
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      for (int i = 0; i < 2; ++i) {
        CPU_SET(i, &cpuset);
      }
      worker_thread->SetThreadAffinity(&cpuset);
      ++dispatched_computation_job_count_;
      --ready_jobs_count_;
      pthread_cond_signal(&worker_thread->thread_can_start);
    }
  }
  pthread_mutex_unlock(&scheduling_critical_section_lock_);
  pthread_mutex_unlock(&computation_job_queue_lock_);
}

void* WorkerManager::ThreadEntryPoint(void* parameters) {
  {
    struct sched_param param;
    param.sched_priority = 0;
    int st = pthread_setschedparam(pthread_self(), SCHED_BATCH, &param);
    if (st != 0) {
      // Scheduling setting goes wrong.
      std::exit(1);
    }
  }
  WorkerThread* worker_thread = reinterpret_cast<WorkerThread*>(parameters);
  worker_thread->Run();
  assert(false);
}

void* WorkerManager::SchedulingEntryPoint(
    void* parameters) {
  WorkerManager* worker_manager = reinterpret_cast<WorkerManager*>(parameters);
  while (true) {
    pthread_mutex_lock(&worker_manager->scheduling_needed_lock_);
    while (!worker_manager->scheduling_needed_) {
      pthread_cond_wait(&worker_manager->scheduling_needed_cond_,
                        &worker_manager->scheduling_needed_lock_);
    }
    worker_manager->scheduling_needed_ = false;
    pthread_mutex_unlock(&worker_manager->scheduling_needed_lock_);

    // Thread scheduling algorithm is invoked when new job is added or new
    // threads turn to idle.
    worker_manager->ScheduleComputationJobs();
  }
  assert(false);
}

void WorkerManager::TriggerScheduling() {
  pthread_mutex_lock(&scheduling_needed_lock_);
  scheduling_needed_ = true;
  pthread_cond_signal(&scheduling_needed_cond_);
  pthread_mutex_unlock(&scheduling_needed_lock_);
}

/*
void WorkerManager::RegisterThread(
    const pthread_t* child_thread, const pthread_t* parenter_thread) {
}

void WorkerManager::RegisterThreadStateStub(
    const pthread_t* thread, const ThreadStateStub* thread_stub) {
}

void WorkerManager::DeregisterThreadStateStub(const pthread_t* thread) {
}

void WorkerManager::AllocateThreadStateStub() {
  pthread_setspecific(*thread_state_stub_key, new ThreadStateStub);
}

WorkerManager::ThreadStateStub* WorkerManager::DetachThreadStateStub() {
  ThreadStateStub* temp =
      static_cast<ThreadStateStub*>(pthread_getspecific(*thread_state_stub_key));
  pthread_setspecific(*thread_state_stub_key, NULL);
  return temp;
}
*/

}  // namespace nimbus
