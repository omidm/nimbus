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
  * Utilities for timing measurement. Not thread-safe.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_TIMER_H_
#define NIMBUS_SRC_SHARED_TIMER_H_

#include <pthread.h>
#include <sys/time.h>
#include <string>
#include <map>

namespace nimbus {

class Timer {
 public:
  // Started in constructor.
  explicit Timer(const std::string& timer_name);
  ~Timer();
  static void Initialize();
  static void Print(const std::string& comment);
  static void Start(const std::string& timer_name);
  static void Stop(const std::string& timer_name);
  static void Reset(const std::string& timer_name);

 private:
  struct TimerState {
    TimerState() {
      on = false;
      total_time = 0;
    }
    struct timeval start_time;
    bool on;
    double total_time;
  };
  std::string timer_name_;

 public:
  typedef std::map<std::string, TimerState> TimerMap;

 private:
  static pthread_mutex_t* lock_;
  static TimerMap* timer_map_;
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_TIMER_H_
