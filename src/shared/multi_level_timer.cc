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
  * Multi level timer, to measure performance of multi-threaded code.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "src/shared/multi_level_timer.h"

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif  // __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <sys/syscall.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <cstdlib>
#include <string>

namespace nimbus {

MultiLevelTimer::MultiLevelTimer() {
  depth_ = 0;
  sum_ = 0;
}

MultiLevelTimer::MultiLevelTimer(std::string name) {
  name_ = name;
  depth_ = 0;
  sum_ = 0;
}

MultiLevelTimer::~MultiLevelTimer() {
}

std::string MultiLevelTimer::name() {
  return name_;
}

void MultiLevelTimer::set_name(std::string name) {
  name_ = name;
}

void MultiLevelTimer::Print(FILE* output) {
  struct timespec now;
#ifdef __MACH__  // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  now.tv_sec = mts.tv_sec;
  now.tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME, &now);
#endif
  int64_t result = sum_ + depth_ * (
            (now.tv_sec - old_timestamp_.tv_sec) * 1e9
            + now.tv_nsec - old_timestamp_.tv_nsec);
    fprintf(output, "tid 0000 name %s time %.9f\n",
        name_.c_str(), static_cast<double>(result) / 1e9);
}

}  // namespace nimbus
