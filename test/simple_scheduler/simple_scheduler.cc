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

//  std::string test = "spawnjob main {0} {1,2} {4,5} ";
//  test += " {6,7,8} {10,20,30} COMP none";
//  SchedulerCommand* gen_com = NULL;
//  bool cond = SchedulerCommand::GenerateSchedulerCommandChild(
//      test, &worker_command_set_, gen_com);
//  if (cond) {
//    std::cout << "*****parsed properly" << std::endl;
//    std::cout << gen_com->toStringWTags() << std::endl;
//    assert(gen_com);
//  } else {
//    std::cout << "ERROR: *****did not generate anything" << std::endl;
//    assert(gen_com == NULL);
//  }

  SchedulerWorkerList::iterator iter = server_->workers()->begin();
//  while (server_->workers()->begin() == server_->workers()->end()) {
//    sleep(1);
//    std::cout << "Waiting for the first worker to connect ..." << std::endl;
//  }

  while (true) {
    sleep(1);
    int ready_num = 0;
    for (iter = server_->workers()->begin();
        iter != server_->workers()->end(); iter++) {
      if ((*iter)->handshake_done()) {
        ready_num++;
      } else {
        ID<worker_id_t> worker_id;
        worker_id.set_elem((*iter)->worker_id());
        std::string ip("you-know");
        ID<port_t> port(0);

        SchedulerCommand* cm = new HandshakeCommand(worker_id, ip, port);
        server_->SendCommand(*iter, cm);
        delete cm;
      }
    }
    if (ready_num >= WORKER_NUM)
      break;

    std::cout << ready_num << " workers are registered, waiting for " <<
      WORKER_NUM - ready_num << " more workers to join ..."  << std::endl;

    SchedulerCommandList* storage = new SchedulerCommandList();
    if (server_->ReceiveCommands(storage, (uint32_t)10)) {
      SchedulerCommandList::iterator iter = storage->begin();
      for (; iter != storage->end(); iter++) {
        SchedulerCommand* comm = *iter;
        ProcessSchedulerCommand(comm);
        delete comm;
      }
    }
    delete storage;
  }

//  std::string str = "spawnjob name:main id:{0} read:{} write:{}";
//  str += " before:{} after:{} type:operation param:none";
//  SchedulerCommand cm(str);
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
        // std::cout << "OMID Received command: " << comm->toString() << std::endl;
        // std::cout << "OMID Sending command: " << comm->toString() << std::endl;
        // server_->SendCommand(*(server_->workers()->begin()), comm);
        ProcessSchedulerCommand(comm);
        delete comm;
      }
    }
    delete storage;
  }
}

