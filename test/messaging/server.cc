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

#include "./create_commands.h"

using ::std::cout;
using ::std::endl;
using namespace nimbus;; // NOLINT
#define PORT 11714

nimbus::SchedulerServer server(PORT);

void* run(void* p) {
  server.Run();
  return NULL;
}

int main(int argc, char *argv[]) {
  dbg_init();
  load_command_prototypes();
  create_commands();
  server.set_worker_command_table(&prototypes);

  cout << "Starting server." << std::endl;
  pthread_t thread;
  pthread_create(&thread, NULL, run, NULL);

  int i = 0;
  while (server.worker_num() < 1) {
    cout << i << ". Waiting for workers" << std::endl;
    sleep(1);
    i++;
  }

  cout << "Starting command sends." << std::endl;
  for (int i = 0; i < NUM_COMMANDS; i++) {
    std::string val = commands[i]->ToNetworkData();
    cout << "  " << i << " of " << NUM_COMMANDS << ": broadcasting command of length " << val.length() << ": ";  // NOLINT
    cout << commands[i]->ToString() << std::endl;
    server.BroadcastCommand(commands[i]);
  }

  cout << "Starting command receives." << std::endl;
  nimbus::SchedulerCommandList commandList;
  for (int i = 0; i < NUM_COMMANDS;) {
    server.ReceiveCommands(&commandList, NUM_COMMANDS);
    while (!commandList.empty()) {
      nimbus::SchedulerCommand* cmd = commandList.front();
      commandList.pop_front();
      std::cout << "  " << i << " of " << NUM_COMMANDS << ": received " << cmd->ToString() << std::endl; // NOLINT
      ++i;
      delete cmd;
    }
  }
}
