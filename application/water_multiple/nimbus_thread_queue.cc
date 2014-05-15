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

#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include "pthread.h"
#include "shared/nimbus.h"

#include "application/water_multiple/nimbus_thread_queue.h"

namespace nimbus {

NimbusThreadQueue::NimbusThreadQueue(int thread_count) : THREAD_QUEUE(-1) {
  pthread_attr_init(&attr);
  pthread_cond_init(&done_condition, 0);
  pthread_cond_init(&todo_condition, 0);
  pthread_mutex_init(&queue_lock, 0);
  // Includes the parent thread.
  thread_quota = thread_count;
  active_threads = thread_count + 1;
  inactive_threads = 0;
  blocked_threads = 0;
  queue_length_ = 0;
  task_threads.resize(thread_count);

  for (int i = 0; i < (int)task_threads.size(); i++){
    task_threads[i] = new NimbusTaskThread(this, i);
    pthread_create(
        &(task_threads[i]->pthread_control_handler), NULL,
        NimbusTaskThread::ThreadRoutine, task_threads[i]);
  }
}
NimbusThreadQueue::~NimbusThreadQueue() {
  EXITER* exiter = new EXITER[task_threads.size()];
  --active_threads;
  ++blocked_threads;
  for (int i = 0; i < (int)task_threads.size(); i++) {
    Queue(&exiter[i]);
  }
  for (int i = 0; i < (int)task_threads.size(); i++) {
    pthread_join(task_threads[i]->pthread_control_handler, NULL);
    delete task_threads[i];
  }
  pthread_cond_destroy(&done_condition);
  pthread_cond_destroy(&todo_condition);
  pthread_mutex_destroy(&queue_lock);
  pthread_attr_destroy(&attr);
  delete[] exiter;
}
void NimbusThreadQueue::Queue(TASK* task) {
  pthread_mutex_lock(&queue_lock);
  ++queue_length_;
  queue.push_back(task);
  if (inactive_threads != 0) {
    pthread_cond_signal(&todo_condition);
  }
  pthread_mutex_unlock(&queue_lock);
}
void NimbusThreadQueue::Wait() {
  pthread_mutex_lock(&queue_lock);
  dbg(DBG_WORKER_BD,
      DBG_WORKER_BD_S"%d tasks before the sync.\n", queue_length_);
  queue_length_ = 0;
  --active_threads;
  ++blocked_threads;
  // In case the parent thread blocked the children thread.
  if (inactive_threads != 0) {
    pthread_cond_signal(&todo_condition);
  }
  while (!queue.empty() || active_threads != 0) {
    pthread_cond_wait(&done_condition, &queue_lock);
  }
  ++active_threads;
  --blocked_threads;
  pthread_mutex_unlock(&queue_lock);
}
int NimbusThreadQueue::Number_Of_Threads() {
  return task_threads.size();
}
int NimbusThreadQueue::get_active_threads() {
  int result;
  pthread_mutex_lock(&queue_lock);
  result = active_threads;
  pthread_mutex_unlock(&queue_lock);
  return result;
}

}  // namespace nimbus
