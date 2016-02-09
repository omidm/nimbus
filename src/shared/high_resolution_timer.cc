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


#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <cassert>
#include <fstream>  // NOLINT

#include "src/shared/high_resolution_timer.h"

namespace nimbus {

HighResolutionTimer::HighResolutionTimer() {
  timer_map_ = new HighResolutionTimerMap;
  pthread_mutex_init(&lock_, 0);
}

HighResolutionTimer::~HighResolutionTimer() {
  delete timer_map_;
  pthread_mutex_destroy(&lock_);
}

void HighResolutionTimer::Start(job_id_t id) {
  assert(timer_map_ != NULL);
  pthread_mutex_lock(&lock_);
  HighResolutionTimerState& state = (*timer_map_)[id];
  pthread_mutex_unlock(&lock_);
  if (!state.on) {
#ifdef __MACH__  // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    state.start_time.tv_sec = mts.tv_sec;
    state.start_time.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, &state.start_time);
#endif
    state.on = true;
  }
}

double HighResolutionTimer::Stop(job_id_t id) {
  assert(timer_map_ != NULL);
  struct timespec end_time;

#ifdef __MACH__  // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  end_time.tv_sec = mts.tv_sec;
  end_time.tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME, &end_time);
#endif

  pthread_mutex_lock(&lock_);
  HighResolutionTimerState& state = (*timer_map_)[id];
  pthread_mutex_unlock(&lock_);
  assert(state.on);
  state.on = false;
  return Diff(&state.start_time, &end_time);
}

double HighResolutionTimer::Diff(timespec *start_time, timespec *end_time) {
  timespec elapsed_time;
  if ((end_time->tv_nsec-start_time->tv_nsec) < 0) {
    elapsed_time.tv_sec = end_time->tv_sec-start_time->tv_sec-1;
    elapsed_time.tv_nsec = 1000000000+end_time->tv_nsec-start_time->tv_nsec;
  } else {
    elapsed_time.tv_sec = end_time->tv_sec-start_time->tv_sec;
    elapsed_time.tv_nsec = end_time->tv_nsec-start_time->tv_nsec;
  }

  double diff  = (static_cast<double>(elapsed_time.tv_sec)) +
    .000000001 * (static_cast<double>(elapsed_time.tv_nsec));

  return diff;
}

}  // namespace nimbus
