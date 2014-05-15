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

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_THREAD_QUEUE_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_THREAD_QUEUE_H_

#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <list>
#include <vector>
#include "shared/nimbus.h"
#include "worker/thread_queue_proto.h"

namespace nimbus {

class NimbusTaskThread;

class NimbusThreadQueue
    : public PhysBAM::THREAD_QUEUE, public ThreadQueueProto {
  friend class NimbusTaskThread;
 public:
  using PhysBAM::THREAD_QUEUE::TASK;
  NimbusThreadQueue(int thread_count);
  virtual ~NimbusThreadQueue();

  virtual void Queue(TASK* task);
  virtual void Wait();
  virtual int Number_Of_Threads();

  struct EXITER : public TASK {
    void Run(const int threadid) {
      pthread_exit(0);
    }
  };
  int get_active_threads();
 private:
  int queue_length_;
  std::vector<NimbusTaskThread*> task_threads;
  std::list<TASK*> queue;
  pthread_attr_t attr;
  pthread_mutex_t queue_lock;
  pthread_cond_t done_condition, todo_condition;
  int active_threads, inactive_threads;
  int blocked_threads, thread_quota;
};

class NimbusTaskThread {
 public:
  NimbusTaskThread(NimbusThreadQueue* nimbus_thread_queue, int thread_rank) {
    queue_ = nimbus_thread_queue;
    thread_rank_ = thread_rank;
  }
  ~NimbusTaskThread() {};
  static void* ThreadRoutine(void* parameter) {
    NimbusTaskThread* nimbus_task_thread =
        static_cast<NimbusTaskThread*>(parameter);
    nimbus_task_thread->Run();
    assert(false);
    return NULL;
  }
  void Run() {
    while (1) {
      pthread_mutex_lock(&queue_->queue_lock);
      while (queue_->queue.empty() ||
             queue_->active_threads > queue_->thread_quota){
        queue_->active_threads--;
        if (queue_->active_threads == 0) {
          pthread_cond_signal(&queue_->done_condition);
        }
        queue_->inactive_threads++;
        pthread_cond_wait(&queue_->todo_condition, &queue_->queue_lock);
        queue_->active_threads++;
        queue_->inactive_threads--;
      }
      NimbusThreadQueue::TASK* task = queue_->queue.front();
      queue_->queue.pop_front();
      pthread_mutex_unlock(&queue_->queue_lock);
      // PhysBAM rank starts from 1.
      task->Run(thread_rank_ + 1);
      delete task;
    }
  }
  pthread_t pthread_control_handler;
 private:
  int thread_rank_;
  NimbusThreadQueue* queue_;
};

}  // namespace nimbus
#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_THREAD_QUEUE_H_

