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

#include <list>

#include "src/worker/worker_manager.h"
#include "src/worker/worker_thread.h"

namespace nimbus {

WorkerThread::WorkerThread(WorkerManager* worker_manager) {
  worker_manager_ = worker_manager;
#ifndef __MACH__
  used_cpu_set_ = NULL;
#endif
}

WorkerThread::~WorkerThread() {
}

#ifndef __MACH__
void WorkerThread::SetThreadAffinity(const cpu_set_t* cpuset) {
  if (used_cpu_set_ == NULL) {
    used_cpu_set_ = new cpu_set_t;
  } else {
    if (CPU_EQUAL(used_cpu_set_, cpuset)) {
      return;
    } else {
      *used_cpu_set_ = *cpuset;
    }
  }
  pthread_setaffinity_np(thread_id, sizeof(cpu_set_t), cpuset);
  for (TaskThreadPool::TaskThreadList::iterator iter =
       allocated_threads.begin();
       iter != allocated_threads.end();
       ++iter) {
    pthread_setaffinity_np((*iter)->thread_handle(),
                           sizeof(cpu_set_t), cpuset);
  }
}
#endif

}  // namespace nimbus
