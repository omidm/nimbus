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

std::map<std::string, double>* Timer::timer_map_ = NULL;

void Timer::Initialize() {
  timer_map_ = new std::map<std::string, double>;
}

void Timer::Print() {
  assert(timer_map_ != NULL);
  std::ofstream output("timing.out");
  for (std::map<std::string, double>::iterator iter = timer_map_->begin();
       iter != timer_map_->end();
       ++iter) {
    output << iter->first << " : " << iter->second << std::endl;
  }
  output.close();
}

Timer::Timer(const std::string& timer_name) {
  assert(timer_map_ != NULL);
  if (timer_map_->find(timer_name) == timer_map_->end()) {
    (*timer_map_)[timer_name] = 0;
  }
  timer_name_ = timer_name;
  on_ = true;
  gettimeofday(&start_time_, NULL);
  total_time_ = 0;
}

Timer::~Timer() {
  Stop();
}

void Timer::Resume() {
  if (on_) {
    return;
  }
  gettimeofday(&start_time_, NULL);
  on_ = true;
}

void Timer::Pause() {
  if (!on_) {
    return;
  }
  struct timeval end_time;
  gettimeofday(&end_time, NULL);
  double temp =
      (static_cast<double>(end_time.tv_sec - start_time_.tv_sec)) +
      .000001 * (static_cast<double>(end_time.tv_usec - start_time_.tv_usec));
  total_time_ += temp;
  on_ = false;
}

void Timer::Stop() {
  Pause();
  (*timer_map_)[timer_name_] += total_time_;
  total_time_ = 0;
  on_ = false;
}

}  // namespace nimbus
