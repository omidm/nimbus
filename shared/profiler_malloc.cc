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
 * Author: Andrew Lim <alim16@stanford.edu>
 */

#include <shared/profiler_malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include <assert.h>
#include <signal.h>
#include <string.h>

namespace nimbus {
  size_t ProfilerMalloc::alloc_ = 0;
  ProfilerMalloc::MallocMap *ProfilerMalloc::alloc_map_ = NULL;
  ProfilerMalloc::ThreadMap *ProfilerMalloc::thread_alloc_map_ = NULL;
  bool ProfilerMalloc::init_ = false;
  bool ProfilerMalloc::map_include_ = false;
  bool ProfilerMalloc::enabled_ = false;

  static pthread_mutex_t mutex_ = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;

  extern "C"
  void *malloc(size_t size) {
    static void *(*mallocp)(size_t size);
    char *error;
    void *ptr;

    /* get address of libc malloc */
    if (!mallocp) {
      mallocp = reinterpret_cast<void *(*)(size_t size)>(dlsym(RTLD_NEXT, "malloc"));
      if ((error = dlerror()) != NULL) {
        fputs(error, stderr);
        exit(1);
      }
    }

    if (!ProfilerMalloc::IsEnabled() || !ProfilerMalloc::IsEnabledTid(pthread_self())) {
      ptr = mallocp(size);
      return ptr;
    }

    ptr = mallocp(size);

    pthread_mutex_lock(&mutex_);
    if (ProfilerMalloc::IsMapInclude()) {
      ProfilerMalloc::InsertAllocPointer(ptr, size);
      pthread_mutex_unlock(&mutex_);
      ProfilerMalloc::IncreaseAlloc(size);
    } else {
      pthread_mutex_unlock(&mutex_);
    }

    return ptr;
  }

  extern "C"
  void free(void *ptr) {
    static void (*freep)(void *ptr);
    char *error;

    /* get address of libc free */
    if (!freep) {
      freep = reinterpret_cast<void (*)(void *ptr)>(dlsym(RTLD_NEXT, "free"));
      if ((error = dlerror()) != NULL) {
        fputs(error, stderr);
        exit(1);
      }
    }

    if (!ProfilerMalloc::IsEnabled() || !ProfilerMalloc::IsEnabledTid(pthread_self())) {
      freep(ptr);
      return;
    }

    pthread_mutex_lock(&mutex_);
    size_t size = ProfilerMalloc::AllocSize(ptr);
    if (size > 0) {
      ProfilerMalloc::DeleteAllocPointer(ptr);
      pthread_mutex_unlock(&mutex_);
      ProfilerMalloc::DecreaseAlloc(size);
    } else {
      pthread_mutex_unlock(&mutex_);
    }
    freep(ptr);
  }

  bool ProfilerMalloc::IsInit() {
    return init_;
  }

  bool ProfilerMalloc::IsMapInclude() {
    return map_include_;
  }

  void ProfilerMalloc::Initialize() {
    init_ = true;
    alloc_map_ = new MallocMap;
    thread_alloc_map_ = new ThreadMap;
    map_include_ = true;
  }

  uint64_t ProfilerMalloc::CurrentAlloc() {
    return alloc_;
  }

  void ProfilerMalloc::IncreaseAlloc(size_t size) {
    pthread_mutex_lock(&mutex_);
    alloc_ += size;
    pthread_mutex_unlock(&mutex_);
    pthread_t tid = pthread_self();
    /* Statistics are only maintained for threads that have been registered with the profiler. */
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      state.curr_alloc +=size;
      if (state.curr_alloc > state.max_alloc) {
        state.max_alloc = state.curr_alloc;
      }
    }
  }

  void ProfilerMalloc::DecreaseAlloc(size_t size) {
    pthread_mutex_lock(&mutex_);
    alloc_ -= size;
    pthread_mutex_unlock(&mutex_);
    pthread_t tid = pthread_self();
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      state.curr_alloc -=size;
    }
  }

  void ProfilerMalloc::InsertAllocPointer(void *ptr, size_t size) {
    map_include_ = false;
    if (alloc_map_ != NULL) {
      (*alloc_map_)[ptr] = size;
    }
    map_include_ = true;
  }

  void ProfilerMalloc::DeleteAllocPointer(void *ptr) {
    if (alloc_map_ != NULL && alloc_map_->find(ptr) != alloc_map_->end()) {
      alloc_map_->erase(ptr);
    }
  }

  size_t ProfilerMalloc::AllocSize(void *ptr) {
    if (alloc_map_ != NULL && alloc_map_->find(ptr) != alloc_map_->end()) {
      size_t size = (*alloc_map_)[ptr];
      return size;
    } else {
      return 0;
    }
  }

  size_t ProfilerMalloc::AllocMax() {
    pthread_t tid = pthread_self();
    return AllocMaxTid(tid);
  }

  size_t ProfilerMalloc::AllocMaxTid(pthread_t tid) {
    pthread_mutex_lock(&mutex_);
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      pthread_mutex_unlock(&mutex_);
      return state.max_alloc;
    } else {
      pthread_mutex_unlock(&mutex_);
      return 0;
    }
  }

  size_t ProfilerMalloc::AllocCurr() {
    pthread_t tid = pthread_self();
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      return state.curr_alloc;
    } else {
      return 0;
    }
  }

  /* Only the main worker thread should call exit. All other threads should have
     terminated at this time. */
  void ProfilerMalloc::Exit() {
    pthread_mutex_lock(&mutex_);
    delete thread_alloc_map_;
    thread_alloc_map_ = NULL;
    delete alloc_map_;
    alloc_map_ = NULL;
    enabled_ = false;
    pthread_mutex_unlock(&mutex_);
  }

  /* Only the main worker thread should call enable. */
  void ProfilerMalloc::Enable() {
    pthread_mutex_lock(&mutex_);
    enabled_ = true;
    pthread_mutex_unlock(&mutex_);
  }

  void ProfilerMalloc::EnableTid(pthread_t tid) {
    // pthread_mutex_lock(&mutex_);
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      state.on = true;
    }
    // pthread_mutex_unlock(&mutex_);
  }

  /* Only the main worker thread should call disable. */
  void ProfilerMalloc::Disable() {
    pthread_mutex_lock(&mutex_);
    enabled_ = false;
    pthread_mutex_unlock(&mutex_);
  }

  void ProfilerMalloc::DisableTid(pthread_t tid) {
    // pthread_mutex_lock(&mutex_);
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      state.on = false;
    }
    // pthread_mutex_unlock(&mutex_);
  }

  bool ProfilerMalloc::IsEnabled() {
    return enabled_;
  }

  bool ProfilerMalloc::IsEnabledTid(pthread_t tid) {
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      return state.on;
    } else {
      return false;
    }
  }

  void ProfilerMalloc::ResetThreadStatistics() {
    // pthread_mutex_lock(&mutex_);
    pthread_t tid = pthread_self();
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      state.curr_alloc = 0;
      state.max_alloc = 0;
    }
    // pthread_mutex_unlock(&mutex_);
  }

  void ProfilerMalloc::ResetThreadStatisticsByTid(pthread_t tid) {
    // pthread_mutex_lock(&mutex_);
    if (thread_alloc_map_ != NULL && thread_alloc_map_->find(tid) != thread_alloc_map_->end()) {
      ThreadAllocState& state = (*thread_alloc_map_)[tid];
      state.curr_alloc = 0;
      state.max_alloc = 0;
    }
    // pthread_mutex_unlock(&mutex_);
  }

  /* Worker thread should register tids for threads in the thread pool. */
  void ProfilerMalloc::RegisterThreads(std::vector<pthread_t> tids) {
    pthread_mutex_lock(&mutex_);
    ProfilerMalloc::Initialize();
    assert(thread_alloc_map_ != NULL);
    for (std::vector<pthread_t>::iterator tid = tids.begin(); tid != tids.end(); ++tid) {
      ThreadAllocState& state = (*thread_alloc_map_)[*tid];
      state.on = true;
    }
    // enabled_ = true;
    pthread_mutex_unlock(&mutex_);
  }

}  // namespace nimbus
