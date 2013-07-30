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
  * A Nimbus worker. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./worker.h"


Worker::Worker(Application* a)
: app(a)
{}

void Worker::run(SchedulerClient* c) {
  std::cout << "Running the Worker" << std::endl;

  client = c;

  // I think the main run loop, which pulls commands off the queue
  // and dispatches them, should be here. Spawn a separate thread that
  // reads commands, or do so whenever a job completes. I suspect
  // the way to do this might be to have a selective thread join.
  // I.e., "spawn these three threads, join on any one completing".
  loadSchedulerCommands();

  // Start the app. The scheduler client is up, so the app can
  // now send commands to the scheduler to dispatch. Start()
  // is the call where the application will initialize data and seed
  // the first set of jobs.
  app->start(c);

  while (true) {
    std::cout << "Worker running core loop." << std::endl;
    // pull commands from SchedulerClient
    // dispatch commands in queue
    //
  }
}

void Worker::loadSchedulerCommands() {
  std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  while (true) {
    std::string word;
    cms >> word;
    if (cms.fail()) {
      break;
    }
    schedulerCmSet.insert(word);
  }
}
