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
  * Nimbus chronometer. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/chronometer.h"

Chronometer::Chronometer() {
    Initialize();
}

Chronometer::~Chronometer() {
}

void Chronometer::Initialize() {
  gettimeofday(&start_time_, NULL);
  timer_ = 0;
  timer_is_on_ = false;
}

struct timeval* Chronometer::start_time() {
  return &start_time_;
}

void Chronometer::set_start_time(struct timeval* time) {
  start_time_.tv_sec = time->tv_sec;
  start_time_.tv_usec = time->tv_usec;
}

double Chronometer::timer() {
  if (timer_is_on_) {
    struct timeval t;
    gettimeofday(&t, NULL);
    double temp  = (static_cast<double>(t.tv_sec - timer_start_time_.tv_sec)) +
      .000001 * (static_cast<double>(t.tv_usec - timer_start_time_.tv_usec));
    return timer_ + temp;
  } else {
    return timer_;
  }
}

void Chronometer::Reset() {
  timer_ = 0;
  timer_is_on_ = false;
}

void Chronometer::Start() {
  timer_ = 0;
  gettimeofday(&timer_start_time_, NULL);
  timer_is_on_ = true;
}

void Chronometer::Resume() {
  if (timer_is_on_)
    return;
  gettimeofday(&timer_start_time_, NULL);
  timer_is_on_ = true;
}

void Chronometer::Stop() {
  if (!timer_is_on_)
    return;
  struct timeval t;
  gettimeofday(&t, NULL);
  timer_  += (static_cast<double>(t.tv_sec - timer_start_time_.tv_sec)) +
  .000001 * (static_cast<double>(t.tv_usec - timer_start_time_.tv_usec));
  timer_is_on_ = false;
}

double Chronometer::GetTime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  double time  = (static_cast<double>(t.tv_sec - start_time_.tv_sec)) +
  .000001 * (static_cast<double>(t.tv_usec - start_time_.tv_usec));
  return time;
}

double Chronometer::GetRawTime() {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  return time_sum;
}

