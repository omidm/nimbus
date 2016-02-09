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
 * 
 * Routines to log memory allocations on a per thread basis. 
 */

#ifdef __MACH__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include <assert.h>
#include <signal.h>
#include <string.h>
#include "src/shared/profiler_malloc.h"

namespace nimbus {
  __thread size_t curr_alloc = 0;
  __thread size_t max_alloc = 0;
  __thread size_t base_alloc = 0;

  void ProfilerMalloc::AddPointerToThreadLocalData(void *ptr) {
#ifdef __MACH__
    size_t size = malloc_size(ptr);
#else
    size_t size = malloc_usable_size(ptr);
#endif
    curr_alloc += size;
    if (curr_alloc > max_alloc) {
      max_alloc = curr_alloc;
    }
  }

  void ProfilerMalloc::RemovePointerFromThreadLocalData(void *ptr) {
#ifdef __MACH__
    size_t size = malloc_size(ptr);
#else
    size_t size = malloc_usable_size(ptr);
#endif
    curr_alloc -= size;
  }

  void ProfilerMalloc::ResetBaseAlloc() {
    base_alloc = curr_alloc;
    max_alloc = curr_alloc;
  }

  size_t ProfilerMalloc::GetMaxAlloc() {
    return max_alloc - base_alloc;
  }

  size_t ProfilerMalloc::GetCurrAlloc() {
    return curr_alloc;
  }

  size_t ProfilerMalloc::GetBaseAlloc() {
    return base_alloc;
  }

  void *ProfilerMalloc::p_malloc(size_t size) {
    static void *(*mallocp)(size_t size) = NULL;
    char *error;
    void *ptr;

    if (!mallocp) {
      mallocp = reinterpret_cast<void *(*)(size_t size)>(dlsym(RTLD_NEXT, "malloc"));
      if ((error = dlerror()) != NULL) {
        fputs(error, stderr);
        exit(1);
      }
    }

    ptr = mallocp(size);
    ProfilerMalloc::AddPointerToThreadLocalData(ptr);
    return ptr;
  }

  void ProfilerMalloc::p_free(void *ptr) {
    static void (*freep)(void *ptr) = NULL;
    char *error;

    if (!freep) {
      freep = reinterpret_cast<void (*)(void *ptr)>(dlsym(RTLD_NEXT, "free"));
      if ((error = dlerror()) != NULL) {
        fputs(error, stderr);
        exit(1);
      }
    }
    ProfilerMalloc::RemovePointerFromThreadLocalData(ptr);
    freep(ptr);
  }

}  // namespace nimbus
