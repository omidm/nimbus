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


#include "shared/multi_level_timer.h"

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif  // __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <sys/syscall.h>

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
  clock_gettime(CLOCK_REALTIME, &now);
  int64_t result = sum_ + depth_ * (
            (now.tv_sec - old_timestamp_.tv_sec) * 1e9
            + now.tv_nsec - old_timestamp_.tv_nsec);
    fprintf(output, "tid 0000 name %s time %.9f\n",
        name_.c_str(), static_cast<double>(result) / 1e9);
}

}  // namespace nimbus
