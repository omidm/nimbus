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

#include <shared/timer.h>
#include <cassert>
#include <fstream>  // NOLINT

namespace nimbus {

Timer::TimerMap* Timer::timer_map_ = NULL;
pthread_mutex_t* Timer::lock_ = NULL;

void Timer::Initialize() {
  if (timer_map_ == NULL) {
    timer_map_ = new TimerMap;
  }
  if (lock_ == NULL) {
    lock_ = new pthread_mutex_t();
    pthread_mutexattr_t attr;
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(lock_, &attr);
  }
}

void Timer::Print(const std::string& comment) {
  assert(timer_map_ != NULL);
  std::ofstream output("timing.out", std::ofstream::app);
  output << comment << std::endl;
  assert(lock_ != NULL);
  pthread_mutex_lock(lock_);
  for (TimerMap::iterator iter = timer_map_->begin();
       iter != timer_map_->end();
       ++iter) {
    Reset(iter->first);
    output << iter->first << " : " << iter->second.total_time << std::endl;
  }
  pthread_mutex_unlock(lock_);
  output.close();
}

Timer::Timer(const std::string& timer_name) {
  assert(timer_map_ != NULL);
  Start(timer_name);
  timer_name_ = timer_name;
}

Timer::~Timer() {
  Stop(timer_name_);
}

void Timer::Start(const std::string& timer_name) {
  assert(timer_map_ != NULL);
  assert(lock_ != NULL);
  pthread_mutex_lock(lock_);
  TimerState& state = (*timer_map_)[timer_name];
  pthread_mutex_unlock(lock_);
  if (!state.on) {
    gettimeofday(&state.start_time, NULL);
    state.on = true;
  }
}

void Timer::Stop(const std::string& timer_name) {
  assert(timer_map_ != NULL);
  struct timeval end_time;
  gettimeofday(&end_time, NULL);
  assert(lock_ != NULL);
  pthread_mutex_lock(lock_);
  TimerState& state = (*timer_map_)[timer_name];
  pthread_mutex_unlock(lock_);
  assert(state.on);
  double temp =
      (static_cast<double>(end_time.tv_sec - state.start_time.tv_sec)) +
      .000001 *
      (static_cast<double>(end_time.tv_usec - state.start_time.tv_usec));
  state.total_time += temp;
  state.on = false;
}

void Timer::Reset(const std::string& timer_name) {
  assert(timer_map_ != NULL);
  assert(lock_ != NULL);
  pthread_mutex_lock(lock_);
  TimerState& state = (*timer_map_)[timer_name];
  pthread_mutex_unlock(lock_);
  if (state.on) {
    Stop(timer_name);
    Start(timer_name);
  }
}

}  // namespace nimbus
