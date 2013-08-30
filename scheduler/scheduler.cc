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
  * Nimbus scheduler. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/scheduler.h"

using namespace nimbus; // NOLINT

Scheduler::Scheduler(port_t p)
: listening_port_(p) {
  appId_ = 0;
}

void Scheduler::Run() {
  Log::dbg_printLine("Running the Scheduler");

  SetupWorkerInterface();
  SetupUserInterface();

  SchedulerCoreProcessor();
}

void Scheduler::ProcessSchedulerCommand(SchedulerCommand* cm) {
  std::string command_name = cm->name();

  if (command_name == "spawnjob") {
    server_->SendCommand(*(server_->workers()->begin()), cm);
  } else if (command_name == "definedata") {
    server_->SendCommand(*(server_->workers()->begin()), cm);
  } else if (command_name == "handshake") {
    HandshakeCommand* hsc = reinterpret_cast<HandshakeCommand*>(cm);

    SchedulerWorkerList::iterator iter = server_->workers()->begin();
    for (; iter != server_->workers()->end(); iter++) {
      if ((*iter)->worker_id() == hsc->worker_id().elem()) {
        // std::string ip =
        //   (*iter)->connection()->socket()->remote_endpoint().address().to_string();
        (*iter)->set_ip(hsc->ip());
        (*iter)->set_port(hsc->port().elem());
        (*iter)->set_handshake_done(true);
        std::cout << "Registered worker, id: " << (*iter)->worker_id() <<
          " IP: " << (*iter)->ip() << " port: " << (*iter)->port() << std::endl;
        break;
      }
    }
  } else {
    std::cout << "ERROR: " << cm->toString() <<
      " have not been implemented in ProcessSchedulerCommand yet." <<
      std::endl;
  }
}

void Scheduler::SetupWorkerInterface() {
  LoadWorkerCommands();
  server_ = new SchedulerServer(listening_port_);
  server_->set_worker_command_set(&worker_command_set_);
  worker_interface_thread_ = new boost::thread(boost::bind(&SchedulerServer::Run, server_));
}

void Scheduler::SetupUserInterface() {
  LoadUserCommands();
  user_interface_thread_ = new boost::thread(
      boost::bind(&Scheduler::GetUserCommand, this));
}

void Scheduler::GetUserCommand() {
  while (true) {
    std::cout << "command: " << std::endl;
    std::string token("runapp");
    std::string str, cm;
    std::vector<int> args;
    getline(std::cin, str);
    parseCommand(str, user_command_set_, cm, args);
    std::cout << "you typed: " << cm << std::endl;
  }
}

void Scheduler::LoadUserCommands() {
  std::stringstream cms("loadapp runapp killapp haltapp resumeapp quit");
  while (true) {
    std::string word;
    cms >> word;
    if (cms.fail()) {
      break;
    }
    user_command_set_.insert(word);
  }
}

void Scheduler::LoadWorkerCommands() {
  // std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  worker_command_set_.insert(
      std::make_pair(std::string("spawnjob"), COMMAND_SPAWN_JOB));
  worker_command_set_.insert(
      std::make_pair(std::string("definedata"), COMMAND_DEFINE_DATA));
  worker_command_set_.insert(
      std::make_pair(std::string("handshake"), COMMAND_HANDSHAKE));
  worker_command_set_.insert(
      std::make_pair(std::string("jobdone"), COMMAND_JOB_DONE));
}


