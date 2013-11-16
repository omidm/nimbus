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
  * This program tests the memory use (and potential memory leaks) of
  * the ClusterMap.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#define NUM_COMPUTERS 39
#define NUM_SWITCHES 17


#include "shared/cluster.h"
#include "shared/dbg.h"

int main(int argc, char *argv[]) {
  nimbus::SchedulerWorker* workers[NUM_COMPUTERS];
  nimbus::cluster_map_id_t computers[NUM_COMPUTERS];
  nimbus::cluster_map_id_t switches[NUM_SWITCHES];
  nimbus::cluster_map_id_t links[NUM_COMPUTERS * NUM_SWITCHES];
  nimbus::cluster_map_id_t swich;

  nimbus::ClusterMap* cm = new nimbus::ClusterMap();

  dbg_init();

  printf("Inserting computers.\n");
  for (int i = 0; i < NUM_COMPUTERS; i++) {
    nimbus::SchedulerWorker* sw = new nimbus::SchedulerWorker(i, NULL, NULL);
    workers[i] = sw;
    computers[i] = cm->CreateComputer(sw,
                                      (lrand48() % 64) + 1,  // mem
                                      (lrand48() % 20) + 1,  // cores
                                      ((lrand48() % 30) + 1) * 100,  // mhz
                                      lrand48() % 30,  // mbps
                                      i);  // ip
  }
  printf("Inserting switches and links.\n");
  for (int i = 0; i < NUM_SWITCHES; i++) {
    swich = cm->CreateSwitch(i,
                             NUM_COMPUTERS,
                             1000,
                             10000,
                             0x7f000001);
    for (int j = 0; j < NUM_COMPUTERS; j++) {
      if (j & 0x1) {
        cm->AddLink(swich, computers[j], 37);
      }
      if (j & 0x02) {
        cm->AddLink(computers[j], swich, 68);
      }
    }
  }

  printf("Deleting computers.\n");
  for (int i = 0; i < NUM_COMPUTERS; i++) {
    cm->Delete(computers[i]);
    delete(workers[i]);
  }
  printf("Leaving switches for destructor.\n");

  delete cm;
}
