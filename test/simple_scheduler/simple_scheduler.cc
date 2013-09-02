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
  * Simple Nimbus scheduler that is supposed to run the application over a
  * single worker. It is intended to check the command exchange interface, the
  * mapping logics and generally the system abstraction soundness.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./simple_scheduler.h"
#define WORKER_NUM 4

SimpleScheduler::SimpleScheduler(unsigned int p)
: Scheduler(p) {
}

void SimpleScheduler::SchedulerCoreProcessor() {
  Log::dbg_printLine("Simple Scheduler Core Processor");

  while (server_->workers()->begin() == server_->workers()->end()) {
    sleep(1);
    std::cout << "Waiting for the first worker to connect ..." << std::endl;
  }

  std::string str = "spawnjob main {0} {} {} {} {} COMP none";
  SchedulerCommand* cm;
  SchedulerCommand::GenerateSchedulerCommandChild(
      str, &worker_command_set_, cm);
  std::cout << "Sending command: " << cm->toString() << std::endl;
  server_->SendCommand(*(server_->workers()->begin()), cm);

  while (true) {
    // sleep(1);
    SchedulerCommandList* storage = new SchedulerCommandList();
    if (server_->ReceiveCommands(storage, (uint32_t)10)) {
      SchedulerCommandList::iterator iter = storage->begin();
      for (; iter != storage->end(); iter++) {
        SchedulerCommand* comm = *iter;
        dbg(DBG_NET, "Iterating across command (of %i) %s\n",
            storage->size(), comm->toString().c_str());
        ProcessSchedulerCommand(comm);
        delete comm;
      }
    }
    delete storage;
  }
}

void SimpleScheduler::ProcessSpawnJobCommand(SpawnJobCommand* cm) {
    server_->SendCommand(*(server_->workers()->begin()), cm);
}

void SimpleScheduler::ProcessDefineDataCommand(DefineDataCommand* cm) {
    server_->SendCommand(*(server_->workers()->begin()), cm);
}

