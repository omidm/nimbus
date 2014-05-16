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
 * Malloc profiler.
 *      
 * Author: Andrew Lim <alim16@stanford.edu>
 */

#ifndef NIMBUS_SHARED_PROFILER_MALLOC_H_
#define NIMBUS_SHARED_PROFILER_MALLOC_H_

#include <boost/thread/recursive_mutex.hpp>
#include <pthread.h>
#include <stdint.h>
#include <inttypes.h>
#include <map>
#include <vector>
#include <string>

namespace nimbus {

class ProfilerMalloc {
 public:
  static bool IsInit();
  static bool IsMapInclude();
  static void Initialize();
  static void IncreaseAlloc(size_t size);
  static void DecreaseAlloc(size_t size);
  static void InsertAllocPointer(void *ptr, size_t size);
  static void DeleteAllocPointer(void *ptr);
  static size_t AllocSize(void *ptr);
  static uint64_t CurrentAlloc();
  static size_t AllocMax();
  static size_t AllocMaxTid(pthread_t tid);
  static size_t AllocCurr();
  static void Exit();
  static void Enable();
  static void EnableTid(pthread_t tid);
  static void Disable();
  static void DisableTid(pthread_t tid);
  static bool IsEnabled();
  static bool IsEnabledTid(pthread_t tid);
  static void ResetThreadStatistics();
  static void ResetThreadStatisticsByTid(pthread_t tid);
  static void RegisterThreads(std::vector<pthread_t> tids);

 private:
  struct ThreadAllocState {
    ThreadAllocState() {
      max_alloc = 0;
      curr_alloc = 0;
      on = true;
    }
    uint64_t max_alloc;
    uint64_t curr_alloc;
    bool on;
  };

 private:
  typedef std::map<void *, size_t> MallocMap;
  typedef std::map<pthread_t, ThreadAllocState> ThreadMap;

 private:
  static uint64_t alloc_;
  static MallocMap *alloc_map_;
  static bool init_;
  static bool map_include_;
  static ThreadMap *thread_alloc_map_;
  static bool enabled_;
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_PROFILER_MALLOC_H_
