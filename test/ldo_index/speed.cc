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
  * This file tests whether the Logical Data Object Index is working
  * properly.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#define DEBUG_MODE

#include <sys/time.h>
#include <iostream> // NOLINT


#include "shared/geometric_region.h"
#include "shared/logical_data_object.h"
#include "shared/ldo_index.h"

#define NUM_LDOS 1000000

using namespace nimbus;  // NOLINT

int main(int argc, char *argv[]) {
  LdoIndex index;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  for (int i = 0; i < NUM_LDOS; i++) {
    std::string str = "pressure";
    GeometricRegion* region = new GeometricRegion(lrand48(), lrand48(),
                                                 lrand48(), lrand48(),
                                                 lrand48(), lrand48());
    LogicalDataObject* obj = new LogicalDataObject(i, str, region);

    index.AddObject(obj);
  }
  gettimeofday(&end, NULL);
  int64_t time = (end.tv_sec - start.tv_sec) * 1000000;
  time += (end.tv_usec - start.tv_usec);
  double dtime = static_cast<double>(time);
  dtime /= 1000000.0;
  printf("Inserting %i LDOs took %f seconds (%f usec/LDO).\n",
         NUM_LDOS,
         dtime,
         ((dtime * 1000000.0) / static_cast<double>(NUM_LDOS)));

  gettimeofday(&start, NULL);
  for (int i = 0; i < NUM_LDOS; i++) {
    LogicalDataObject* obj = index.SpecificObject(i);
    if (obj == NULL) {
      fprintf(stderr, "Trying to delete object %i fails because index has no such object.\n", i);
    } else {
      index.RemoveObject(i);
    }
    delete obj;
  }
  gettimeofday(&end, NULL);
  time = (end.tv_sec - start.tv_sec) * 1000000;
  time += (end.tv_usec - start.tv_usec);
  dtime = static_cast<double>(time);
  dtime /= 1000000.0;
  printf("Deleting %i LDOs took %f seconds (%f usec/LDO).\n",
         NUM_LDOS,
         dtime,
         (dtime * 1000000.0) / static_cast<double>(NUM_LDOS));
}
