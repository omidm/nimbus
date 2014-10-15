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
 * Author: Hang Qu<quhang@stanford.edu>
 */

#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "worker/task_thread_pool.h"

namespace nimbus {

ExceptionTaskThreadExit:: ExceptionTaskThreadExit(void* user_result) {
  user_result_ = user_result;
}
void* ExceptionTaskThreadExit::user_result() {
  return user_result_;
}

TaskThreadWrapper::TaskThreadWrapper() {
  pthread_mutex_init(&internal_lock_, NULL);
  pthread_cond_init(&internal_cond_, NULL);
  user_result_ = NULL;
  user_parameter_ = NULL;
  user_function_ = NULL;
  has_work_to_do_ = false;
}
TaskThreadWrapper::~TaskThreadWrapper() {
  pthread_mutex_destroy(&internal_lock_);
  pthread_cond_destroy(&internal_cond_);
}
void* TaskThreadWrapper::ThreadRoutine(void* parameter) {
  {
    struct sched_param param;
    param.sched_priority = 0;
    int st = pthread_setschedparam(pthread_self(), SCHED_BATCH, &param);
    if (st != 0) {
      // Scheduling setting goes wrong.
      std::exit(1);
    }
  }
  TaskThreadWrapper* task_thread =
      static_cast<TaskThreadWrapper*>(parameter);
  task_thread->MainLoop();
  return NULL;
}
void TaskThreadWrapper::Run(
    void* (*user_function)(void* parameter), void* user_parameter) {
  pthread_mutex_lock(&internal_lock_);
  assert(!has_work_to_do_);
  has_work_to_do_ = true;
  user_function_ = user_function;
  user_parameter_ = user_parameter;
  pthread_cond_broadcast(&internal_cond_);
  pthread_mutex_unlock(&internal_lock_);
}
void* TaskThreadWrapper::Join() {
  void* result = NULL;
  pthread_mutex_lock(&internal_lock_);
  while (has_work_to_do_) {
    pthread_cond_wait(&internal_cond_, &internal_lock_);
  }
  result = user_result_;
  pthread_mutex_unlock(&internal_lock_);
  return result;
}
bool TaskThreadWrapper::Initialize() {
  int status = pthread_create(&thread_handle_, NULL,
                              ThreadRoutine, reinterpret_cast<void*>(this));
  return status == 0;
}
pthread_t TaskThreadWrapper::thread_handle() const {
  return thread_handle_;
}
void TaskThreadWrapper::MainLoop() {
  while (true) {
    pthread_mutex_lock(&internal_lock_);
    while (!has_work_to_do_) {
      pthread_cond_wait(&internal_cond_, &internal_lock_);
    }
    user_result_ = NULL;
    pthread_mutex_unlock(&internal_lock_);
    try {
      user_result_ = user_function_(user_parameter_);
    } catch(ExceptionTaskThreadExit finish_exception) {
      user_result_ = finish_exception.user_result();
    }
    pthread_mutex_lock(&internal_lock_);
    has_work_to_do_ = false;
    pthread_cond_broadcast(&internal_cond_);
    pthread_mutex_unlock(&internal_lock_);
  }  // Finishes loop.
}

TaskThreadPool::TaskThreadPool() {
  pthread_mutex_init(&list_lock_, NULL);
}
TaskThreadPool::~TaskThreadPool() {
  pthread_mutex_destroy(&list_lock_);
}
bool TaskThreadPool::AllocateTaskThreads(
    const int thread_num, TaskThreadList* thread_list) {
  pthread_mutex_lock(&list_lock_);
  for (int i = 0; i < thread_num; ++i) {
    if (free_list_.empty()) {
      TaskThreadWrapper* temp_thread = new TaskThreadWrapper;
      bool status = temp_thread->Initialize();
      assert(status);
      free_list_.push_back(temp_thread);
    }
    thread_list->push_back(free_list_.front());
    free_list_.pop_front();
  }
  pthread_mutex_unlock(&list_lock_);
  return true;
}
void TaskThreadPool::FreeTaskThreads(const TaskThreadList& thread_list) {
  pthread_mutex_lock(&list_lock_);
  for (TaskThreadList::const_iterator iter = thread_list.begin();
       iter != thread_list.end();
       ++iter) {
    free_list_.push_back(*iter);
  }
  pthread_mutex_unlock(&list_lock_);
}

}  // namespace nimbus
