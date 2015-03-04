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
  * Testing the correctness of unordered_map_view.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <sched.h>
#include <sys/time.h>
#include <unistd.h>
#include <boost/unordered_map.hpp>
#include <shared/unordered_map_view.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>

void set_aff(int proc_num) {
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(proc_num, &set);
  if (sched_setaffinity(getpid(), sizeof(cpu_set_t), &set)) {
    perror("sched_setaffinity");
  }
}

int main() {
  srand(1991);
  struct timespec t;

  set_aff(0);
  boost::unordered_map<int64_t, int64_t> map;
  boost::unordered_map_view<int64_t, int64_t> view;

  int64_t total = 64000;
  for (int64_t i = 0; i < total; ++i) {
    map[i] = i + 1;
  }
  view.snapshot(map);

  int64_t sample = 100;
  int64_t loop = 10;
  set_aff(2);
  clock_gettime(CLOCK_REALTIME, &t);

  double time_1 = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  for (int64_t i = 0; i < sample; ++i) {
    int64_t r = rand() % total;  // NOLINT
    for (int64_t j = 0; j < loop; ++j) {
      int64_t k = (r+j) % total;
      assert(view[k] == k + 1);
      // assert(map[k] == k + 1);
    }
  }
  clock_gettime(CLOCK_REALTIME, &t);
  double time_2 = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);

  set_aff(4);
  clock_gettime(CLOCK_REALTIME, &t);
  double time_3 = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  int64_t sum;
  for (int64_t i = 0; i < sample; ++i) {
    int64_t r = rand() % total;  // NOLINT
    for (int64_t j = 0; j < loop; ++j) {
      int64_t k = (r+j) % total;
      sum += k;
    }
  }
  clock_gettime(CLOCK_REALTIME, &t);
  double time_4 = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);

  printf("time per query: %.9f ns\n",
         ((time_2 - time_1)-(time_4 - time_3))/sample/loop*1e9);

  return 0;
}
