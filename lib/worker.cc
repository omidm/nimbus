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

#include "lib/worker.h"


Worker::Worker(unsigned int p, Application* a)
: port(p),
  app(a) {
}

void Worker::run() {
  std::cout << "Running the Worker" << std::endl;

  setupSchedulerInterface();
  app->start(client);

  workerCoreProcessor();
}

void Worker::processSchedulerCommand(SchedulerCommand cm) {
  std::string command_name = cm.getName();

  if (command_name == "spawnjob") {
      std::string job_name = cm.getParameters()["name"].getArg();
      Job * j = app->cloneJob(job_name);

      std::vector<Data*> da;
      IDSet::IDSetIter iter;

      IDSet read = cm.getParameters()["read"].getIDSet();
      for (iter = read.set.begin(); iter != read.set.end(); iter++)
        da.push_back(dataMap[*iter]);

      IDSet write = cm.getParameters()["write"].getIDSet();
      for (iter = write.set.begin(); iter != write.set.end(); iter++)
        da.push_back(dataMap[*iter]);

      std::string param = cm.getParameters()["param"].getArg();

      j->Execute(param, da);
  } else if (command_name == "definedata") {
      std::string data_name = cm.getParameters()["name"].getArg();
      Data * d = app->cloneData(data_name);
      d->Create();
      IDSet id = cm.getParameters()["id"].getIDSet();
      addData(*(id.set.begin()), d);
  }
}

void Worker::setupSchedulerInterface() {
  loadSchedulerCommands();
  client = new SchedulerClient(port);
  client->run();
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

void Worker::addJob(int id, Job* job) {
  jobMap[id] = job;
}

void Worker::delJob(int id) {
}

void Worker::addData(int id, Data* data) {
  dataMap[id] = data;
}

void Worker::delData(int id) {
}

