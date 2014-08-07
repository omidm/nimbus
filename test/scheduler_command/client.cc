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

#include "shared/scheduler_client.h"
#include "shared/scheduler_command.h"
#include "shared/protobuf_compiled/commands.pb.h"

using ::std::cout;
using ::std::endl;
using namespace nimbus; // NOLINT

#define PORT 11714
#define NUM_COMMANDS 2


nimbus::SchedulerClient client("127.0.0.1", 11714);
nimbus::SchedulerCommand::PrototypeTable prototypes;

void load_command_prototypes();

void* run(void* p) {
  client.Run();
  return NULL;
}

int main(int argc, char *argv[]) {
  dbg_init();
  load_command_prototypes();
  nimbus::SchedulerCommand* commands[NUM_COMMANDS];
  pthread_t thread;

  cout << "Starting client." << std::endl;
  client.set_scheduler_command_table(&prototypes);
  client.Run();

  // pthread_create(&thread, NULL, run, NULL);

  for (int i = 0; i < NUM_COMMANDS; i++) {
    nimbus::SchedulerCommand* cmd = client.ReceiveCommand();
    cout << "Received command " << cmd->toStringWTags() << std::endl;
  }

  while (1) {
    sleep(1);
  }
}

void load_command_prototypes() {
  prototypes[SchedulerPBuf_Type_HANDSHAKE] = new HandshakeCommand();
  prototypes[SchedulerPBuf_Type_JOB_DONE] = new JobDoneCommand();
  prototypes[SchedulerPBuf_Type_EXECUTE_COMPUTE] = new ComputeJobCommand();
  prototypes[SchedulerPBuf_Type_SPAWN_COMPUTE] = new SpawnComputeJobCommand();
  prototypes[SchedulerPBuf_Type_SPAWN_COPY] = new SpawnCopyJobCommand();
  prototypes[SchedulerPBuf_Type_CREATE_DATA] = new CreateDataCommand();
  prototypes[SchedulerPBuf_Type_REMOTE_SEND] = new RemoteCopySendCommand();
  prototypes[SchedulerPBuf_Type_REMOTE_RECEIVE] = new RemoteCopyReceiveCommand(); // NOLINT
  prototypes[SchedulerPBuf_Type_LOCAL_COPY] = new LocalCopyCommand();
  prototypes[SchedulerPBuf_Type_LDO_ADD] = new LdoAddCommand();
  prototypes[SchedulerPBuf_Type_LDO_REMOVE] = new LdoRemoveCommand();
  prototypes[SchedulerPBuf_Type_PARTITION_ADD] = new PartitionAddCommand();
  prototypes[SchedulerPBuf_Type_PARTITION_REMOVE] = new PartitionRemoveCommand();
  prototypes[SchedulerPBuf_Type_TERMINATE] = new TerminateCommand();
}
