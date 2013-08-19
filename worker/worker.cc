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

#include "worker/worker.h"

using namespace nimbus; // NOLINT

Worker::Worker(unsigned int p, Application* a)
: port_(p),
  application_(a) {
    log.InitTime();
}

void Worker::run() {
  std::cout << "Running the Worker" << std::endl;

  setupSchedulerInterface();
  application_->start(client_);

  workerCoreProcessor();
}

void Worker::processSchedulerCommand(SchedulerCommand* cm) {
  std::string command_name = cm->name();

  if (command_name == "spawnjob") {
      std::string job_name = (*(cm->parameters()))["name"].value();
      Job * j = application_->cloneJob(job_name);

      std::vector<Data*> da;
      IDSet<data_id_t>::IDSetIter iter;

      IDSet<data_id_t>* read = (*(cm->parameters()))["read"].identifier_set();
      for (iter = read->begin(); iter != read->end(); iter++)
        da.push_back(data_map_[*iter]);

      IDSet<data_id_t>* write = (*(cm->parameters()))["write"].identifier_set();
      for (iter = write->begin(); iter != write->end(); iter++)
        da.push_back(data_map_[*iter]);

      std::string param = (*(cm->parameters()))["param"].value();

      IDSet<data_id_t>* id_set = (*(cm->parameters()))["id"].identifier_set();
      job_id_t id = *(id_set->begin());

      log.StartTimer();
      j->execute(param, da);
      log.StopTimer();

      char buff[MAX_BUFF_SIZE];
      snprintf(buff, sizeof(buff),
          "Execute Job, name: %25s  id: %4ld  length(ms): %6.3lf  time(s): %6.3lf",
          job_name.c_str(), id, 1000 * log.timer(), log.GetTime());

      log.writeToFile(std::string(buff), LOG_INFO);
  } else if (command_name == "definedata") {
      std::string data_name = (*(cm->parameters()))["name"].value();
      Data * d = application_->cloneData(data_name);

      IDSet<data_id_t>* id_set = (*(cm->parameters()))["id"].identifier_set();
      data_id_t id = *(id_set->begin());

      log.StartTimer();
      d->create();
      d->set_id(id);
      addData(d);
      log.StopTimer();

      char buff[MAX_BUFF_SIZE];
      snprintf(buff, sizeof(buff),
          "Create Data, name: %25s  id: %4ld  length(ms): %6.3lf  time(s): %6.3lf",
          data_name.c_str(), id, 1000 * log.timer(), log.GetTime());

      log.writeToFile(std::string(buff), LOG_INFO);
  }
}

void Worker::setupSchedulerInterface() {
  loadSchedulerCommands();
  client_ = new SchedulerClient(port_);
  client_->run();
}

void Worker::loadSchedulerCommands() {
  std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  while (true) {
    std::string word;
    cms >> word;
    if (cms.fail()) {
      break;
    }
    scheduler_command_set_.insert(word);
  }
}

void Worker::addJob(Job* job) {
  job_map_[job->id()] = job;
}

void Worker::deleteJob(int id) {
}

void Worker::addData(Data* data) {
  data_map_[data->id()] = data;
}

void Worker::deleteData(int id) {
}

