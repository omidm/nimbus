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

#include "lib/scheduler.h"

Scheduler::Scheduler(unsigned int port) {
  appId = 0;
  this->port = port;
}

void Scheduler::run() {
  cout << "Running the Scheduler" << endl;

  loadUserCommands();
  loadWorkerCommands();

  user_interface_thread = new boost::thread(
      boost::bind(&Scheduler::setupUI, this));

  worker_interface_thread = new boost::thread(
      boost::bind(&Scheduler::setupWI, this));

  user_interface_thread->join();
  worker_interface_thread->join();
}

void Scheduler::setupWI() {
  server = new SchedulerServer(port, this);
  server->run();
}

void Scheduler::setupUI() {
  while (true) {
    cout << "command: ";
    string token("runapp");
    string str, cm;
    vector<int> args;
    getline(cin, str);
    parseCommand(str, userCmSet, cm, args);
    cout << "you typed: " << cm << endl;
  }
}

void Scheduler::loadUserCommands() {
  stringstream cms("loadapp runapp killapp haltapp resumeapp quit");
  while (true) {
    string word;
    cms >> word;
    if (cms.fail()) {
      break;
    }
    userCmSet.insert(word);
  }
}

void Scheduler::loadWorkerCommands() {
  stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  while (true) {
    string word;
    cms >> word;
    if (cms.fail()) {
      break;
    }
    workerCmSet.insert(word);
  }
}


