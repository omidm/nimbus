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

#ifndef NIMBUS_SRC_SHARED_MULTI_LEVEL_TIMER_H_
#define NIMBUS_SRC_SHARED_MULTI_LEVEL_TIMER_H_

#include <sys/types.h>

#include <cassert>
#include <cstdio>
#include <string>
#include <map>
#include <utility>

namespace nimbus {

class MultiLevelTimer {
  public:
    MultiLevelTimer();
    explicit MultiLevelTimer(std::string name);
    ~MultiLevelTimer();

    void Print(FILE* output = stdout);

    std::string name();
    void set_name(std::string name);

    inline void Start(size_t d) {
      assert(d > 0);
      clock_gettime(CLOCK_REALTIME, &(new_timestamp_));
      sum_ += depth_ * (
            (new_timestamp_.tv_sec - old_timestamp_.tv_sec) * 1e9
            + new_timestamp_.tv_nsec - old_timestamp_.tv_nsec);
      old_timestamp_ = new_timestamp_;
      depth_ += d;
    }

    inline void Stop(size_t d) {
      assert(d > 0);
      clock_gettime(CLOCK_REALTIME, &(new_timestamp_));
      sum_ += depth_ * (
            (new_timestamp_.tv_sec - old_timestamp_.tv_sec) * 1e9
            + new_timestamp_.tv_nsec - old_timestamp_.tv_nsec);
      old_timestamp_ = new_timestamp_;
      depth_ -= d;
      assert(depth_ >= 0);
    }

    inline int64_t Read() {
      struct timespec now;
      clock_gettime(CLOCK_REALTIME, &(now));
      return sum_ + depth_ * (
            (now.tv_sec - old_timestamp_.tv_sec) * 1e9
            + now.tv_nsec - old_timestamp_.tv_nsec);
    }

  private:
    std::string name_;
    struct timespec old_timestamp_;
    struct timespec new_timestamp_;
    int64_t depth_;
    int64_t sum_;
};


void InitializeTimers();

}  // namespace nimbus


#endif  // NIMBUS_SRC_SHARED_MULTI_LEVEL_TIMER_H_
