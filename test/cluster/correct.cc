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
  * This program tests whether the cluster map is working properly.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#define NUM_COMPUTERS 50
#define SWITCH_ID 6505

#include "shared/cluster.h"
#include "shared/dbg.h"

int main(int argc, char *argv[]) {
  nimbus::cluster_map_id_t computers[NUM_COMPUTERS];
  nimbus::cluster_map_id_t links[NUM_COMPUTERS];
  nimbus::cluster_map_id_t swich;

  nimbus::ClusterMap* cm = new nimbus::ClusterMap();

  for (int i = 0; i < NUM_COMPUTERS; i++) {
    nimbus::SchedulerWorker* sw = new nimbus::SchedulerWorker(i, NULL, NULL);
    computers[i] = cm->CreateComputer(sw,
                                      (lrand48() % 64) + 1,  // mem
                                      (lrand48() % 20) + 1,  // cores
                                      ((lrand48() % 30) + 1) * 100,  // mhz
                                      lrand48() % 30,  // mbps
                                      i);  // ip
  }

  printf("Testing computer insertion.\n");
  for (int i = 0; i < NUM_COMPUTERS; i++) {
    nimbus::cluster_map_id_t cmid = cm->LookupWorkerId(i);
    if (cmid != computers[i] || true) {
      printf("Worker id %i has cluster map id %i: %s.\n",
             i, cmid, (cmid == computers[i])? "CORRECT":"WRONG");
    }
    nimbus::Computer* comp = cm->LookupComputer(cmid);
    nimbus::Node* comp2 = cm->LookupNode(cmid);
    if (comp != comp2 || true) {
      printf("Looking up node and computer by cmid %i: %s\n",
             cmid, (comp == comp2)? "MATCH":"FAIL");
    }
  }

  swich = cm->CreateSwitch(SWITCH_ID,
                           NUM_COMPUTERS,
                           1000,
                           10000,
                           0x7f000001);

  printf("Testing links.\n");
  for (int i = 0; i < NUM_COMPUTERS; i++) {
    links[i] = cm->AddLink(computers[i], swich, 100);
    links[i] = cm->AddLink(swich, computers[i], 100);
  }

  for (int i = 0; i < NUM_COMPUTERS; i++) {
    nimbus::Computer* comp = cm->LookupComputer(computers[i]);
    nimbus::LinkPtrSet* out = comp->out_links();
    nimbus::LinkPtrSet::iterator it = out->begin();
    for (; it != out->end(); ++it) {
      nimbus::Link* link = *it;
      nimbus::Node* dest = link->destination();
      if (dest->type() != nimbus::CLUSTER_SWITCH ||
          dest->id() != swich) {
        printf("Destination of edge from computer %i is incorrect. It should be to a switch but isn't.\n", i);  // NOLINT
      }
    }
    cm->Delete(computers[i]);
  }
}
