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
  * Some helper functions used for implementation.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "src/shared/helpers.h"

namespace nimbus {

std::string exec(const char* cmd) {
  FILE* pipe = popen(cmd, "r");
  if (!pipe) return "ERROR";
  char buffer[128];
  std::string result = "";
  while (!feof(pipe)) {
    if (fgets(buffer, 128, pipe) != NULL)
      result += buffer;
  }
  pclose(pipe);
  if (result.size() > 0) {
    if (*result.rbegin() == '\n') {
      result.erase(--result.end());
    }
  }
  return result;
}

std::string int2string(uint64_t num) {
  std::stringstream ss;
  ss << num;
  return ss.str();
}

void spin_wait(size_t spin_wait_us, size_t per_iter_mul) {
  if (spin_wait_us == 0) {
    return;
  }

  struct timeval start;
  gettimeofday(&start, NULL);
  size_t dummy = 1;
  size_t elapsed = 0;
  size_t loop_count = 0;
  while (elapsed < spin_wait_us) {
    for (size_t i = 0; i < per_iter_mul; ++i) {
      dummy *= i;
    }
    ++loop_count;
    struct timeval end;
    gettimeofday(&end, NULL);
    elapsed =
      (static_cast<size_t>(end.tv_sec - start.tv_sec)) * 1000000 +
      (static_cast<size_t>(end.tv_usec - start.tv_usec));
  }

  return;
}


}  // namespace nimbus

