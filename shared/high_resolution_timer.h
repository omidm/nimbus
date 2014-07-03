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
 * High resolution timer for profiling job times.
 *                                                                                                     
 * Author: Andrew Lim <alim16@stanford.edu>
 */

#ifndef NIMBUS_SHARED_HIGH_RESOLUTION_TIMER_H_
#define NIMBUS_SHARED_HIGH_RESOLUTION_TIMER_H_

#include <sys/time.h>
#include <map>
#include "shared/id.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class HighResolutionTimer {
 public:
  explicit HighResolutionTimer();
  ~HighResolutionTimer();
  void Start(job_id_t id);
  double Stop(job_id_t id);

 private:
  struct HighResolutionTimerState {
    HighResolutionTimerState() {
      on = false;
    }
    bool on;
    struct timespec start_time;
  };

 private:
  typedef std::map<job_id_t, HighResolutionTimerState> HighResolutionTimerMap;
  static double Diff(timespec *start_time, timespec *end_time);

 private:
  HighResolutionTimerMap* timer_map_;
  pthread_mutex_t lock_;
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_HIGH_RESOLUTION_TIMER_H_
