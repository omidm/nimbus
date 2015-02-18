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
  * Nimbus log interface.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */


#include "shared/fast_log.hh"

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif  // __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <sys/syscall.h>

#include <cstdlib>
#include <string>

namespace {
pthread_mutex_t map_lock;
}  // namespace
namespace nimbus {

namespace timer {

TimersMap timers_map;
CountersMap counters_map;
pthread_key_t timer_keys[kMaxTimer];
pthread_key_t counter_keys[kMaxCounter];

std::string TimerName(TimerType timer_type) {
  switch (timer_type) {
    case kTotal: return "kTotal";
    case kExecuteComputationJob: return "kExecuteComputationJob";
    case kExecuteCopyJob: return "kExecuteCopyJob";
    case kAssemblingCache: return "kAssemblingCache";
    case kSumCyclesTotal: return "kSumCyclesTotal";
    case kSumCyclesBlock: return "kSumCyclesBlock";
    case kSumCyclesRun: return "kSumCyclesRun";
    default: return "Unknown";
  };
}

std::string CounterName(CounterType counter_type) {
  switch (counter_type) {
    case kNumIteration: return "kNumIteration";
    default: return "Unknown";
  }
}

void InitializeKeys() {
  pthread_mutex_init(&map_lock, NULL);
  for (int i = 0; i < kMaxTimer; ++i) {
    pthread_key_create(&(timer_keys[i]), NULL);
  }
  for (int i = 0; i < kMaxCounter; ++i) {
    pthread_key_create(&(counter_keys[i]), NULL);
  }
}

void InitializeTimers() {
  pid_t pid = syscall(SYS_gettid);
  for (int i = 0; i < kMaxTimer; ++i) {
    TimerRecord* record = new(malloc(4096)) TimerRecord;
    pthread_setspecific(timer_keys[i], record);
    pthread_mutex_lock(&map_lock);
    timers_map[std::make_pair(pid, static_cast<TimerType>(i))] = record;
    pthread_mutex_unlock(&map_lock);
  }
  for (int i = 0; i < kMaxCounter; ++i) {
    CounterRecord* record = new(malloc(4096)) CounterRecord;
    pthread_setspecific(counter_keys[i], record);
    pthread_mutex_lock(&map_lock);
    counters_map[std::make_pair(pid, static_cast<CounterType>(i))] = record;
    pthread_mutex_unlock(&map_lock);
  }
}

void PrintTimerSummary(FILE* output) {
  pthread_mutex_lock(&map_lock);
  struct timespec now;
  clock_gettime(CLOCK_REALTIME, &now);
  for (TimersMap::iterator iter = timers_map.begin();
       iter != timers_map.end(); ++iter) {
    TimerRecord* record = iter->second;
    assert(record != NULL);
    if (record->sum == 0 && record->depth == 0) {
      continue;
    }
    int64_t result = record->sum +
        record->depth * (
            (now.tv_sec - record->old_timestamp.tv_sec) * 1e9
            + now.tv_nsec - record->old_timestamp.tv_nsec);
    fprintf(output, "tid %d name %s time %.9f\n",
           iter->first.first, TimerName(iter->first.second).c_str(),
           static_cast<double>(result) / 1e9);
  }
  for (CountersMap::iterator iter = counters_map.begin();
       iter != counters_map.end(); ++iter) {
    CounterRecord* record = iter->second;
    assert(record != NULL);
    if (record->sum == 0) {
      continue;
    }
    fprintf(output, "tid %d name %s counts %"PRId64"\n",
           iter->first.first, CounterName(iter->first.second).c_str(),
           iter->second->sum);
  }
  pthread_mutex_unlock(&map_lock);
}

}  // namespace timer

}  // namespace nimbus
