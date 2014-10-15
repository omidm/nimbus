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

#ifndef NIMBUS_WORKER_TASK_THREAD_POOL_H_
#define NIMBUS_WORKER_TASK_THREAD_POOL_H_

#include <pthread.h>

#include <exception>
#include <list>

namespace nimbus {

// Exception to indicate an task thread finishes.
class ExceptionTaskThreadExit : public std::exception {
 public:
  explicit ExceptionTaskThreadExit(void* user_result);
  void* user_result();
 private:
  void* user_result_;
};

// A wrapper over pthread threads.
class TaskThreadWrapper {
 public:
  TaskThreadWrapper();
  ~TaskThreadWrapper();
  // Runs the user_function with user_parameter passed as parameter.
  void Run(void* (*user_function)(void* parameter), void* user_parameter);
  // Joins the task thread.
  void* Join();
  // Initializes the task thread.
  bool Initialize();
  pthread_t thread_handle() const;

 private:
  pthread_mutex_t internal_lock_;
  pthread_cond_t internal_cond_;
  pthread_t thread_handle_;
  void* user_parameter_;
  void* (*user_function_)(void* parameter);
  void* user_result_;
  bool has_work_to_do_;
  static void* ThreadRoutine(void* parameter);
  void MainLoop();
};

class TaskThreadPool {
 public:
  typedef std::list<TaskThreadWrapper*> TaskThreadList;
  TaskThreadPool();
  ~TaskThreadPool();
  bool AllocateTaskThreads(const int thread_num, TaskThreadList* thread_list);
  void FreeTaskThreads(const TaskThreadList& thread_list);

 private:
  pthread_mutex_t list_lock_;
  std::list<TaskThreadWrapper*> free_list_;
};

}  // namespace nimbus
#endif  // NIMBUS_WORKER_TASK_THREAD_POOL_H_
