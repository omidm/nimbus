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
  * Testing command parsing and communication between a scheduler and
  * a worker.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include <pthread.h>
#include <iostream>  // NOLINT

#include "shared/scheduler_server.h"
#include "shared/scheduler_command.h"

using ::std::cout;
using ::std::endl;
using namespace nimbus;; // NOLINT
#define PORT 11714
#define NUM_COMMANDS 2

nimbus::SchedulerServer server(PORT);
nimbus::SchedulerCommand* commands[NUM_COMMANDS];

void* run(void* p) {
  server.Run();
  return NULL;
}

void create_commands();

int main(int argc, char *argv[]) {
  create_commands();

  cout << "Starting server." << std::endl;
  pthread_t thread;
  pthread_create(&thread, NULL, run, NULL);

  int i = 0;
  while (server.worker_num() < 1) {
    cout << i << ". Waiting for workers" << std::endl;
    sleep(1);
    i++;
  }

  for (int i = 0; i < NUM_COMMANDS; i++) {
    std::string val = commands[i]->toString();
    cout << "Broadcasting command of length " << val.length() << ": ";
    cout << commands[i]->toStringWTags() << std::endl;
    server.BroadcastCommand(commands[i]);
  }

  nimbus::SchedulerCommandList commandList;
  server.ReceiveCommands(&commandList, NUM_COMMANDS);

  nimbus::SchedulerCommandList::iterator it = commandList.begin();
  while (it != commandList.end()) {
    nimbus::SchedulerCommand* cmd = *it;
    cout << cmd->toStringWTags() << std::endl;
    ++it;
  }
}

IDSet<logical_data_id_t> ldo_set_a;
IDSet<logical_data_id_t> ldo_set_b;
IDSet<job_id_t> job_set_a;
IDSet<job_id_t> job_set_b;
ID<logical_data_id_t> ldo_a;
ID<logical_data_id_t> ldo_b;
ID<job_id_t> job_id;
ID<job_id_t> parent;
ID<job_id_t> future;
bool sterile;
Parameter params;
std::string name = "test-string";

void initalize_sets() {
  for (int i = 0; i < 125; i++) {
    ldo_set_a.insert(10000 + i);
  }
  for (int i = 0; i < 119; i++) {
    ldo_set_b.insert(11000 + i);
  }
  for (int i = 0; i < 21; i++) {
    job_set_a.insert(20000 + i);
  }
  for (int i = 0; i < 34; i++) {
    job_set_b.insert(21000 + i);
  }

  ldo_a.set_elem(10125);
  ldo_b.set_elem(11119);
  job_id.set_elem(2014);
  parent.set_elem(1940);
  future.set_elem(1977);
  sterile = false;

  SerializedData d("serialized data");
  params.set_ser_data(d);
}

void create_commands() {
  commands[0] = new SpawnComputeJobCommand((const std::string)name,
                                           (const ID<job_id_t>)job_id,
                                           (const IDSet<logical_data_id_t>)ldo_set_a,
                                           (const IDSet<logical_data_id_t>)ldo_set_b,
                                           (const IDSet<job_id_t>)job_set_a,
                                           (const IDSet<job_id_t>)job_set_b,
                                           (const ID<job_id_t>)parent,
                                           (const ID<job_id_t>)future,
                                           (const bool)sterile,
                                           (const Parameter)params);
  commands[1] = new SpawnCopyJobCommand((const ID<job_id_t>)job_id,
                                        (const ID<logical_data_id_t>)ldo_a,
                                        (const ID<logical_data_id_t>)ldo_b,
                                        (const IDSet<job_id_t>)job_set_a,
                                        (const IDSet<job_id_t>)job_set_b,
                                        (const ID<job_id_t>)parent);
}
