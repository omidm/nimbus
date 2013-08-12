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
  * Nimbus abstraction of an application. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "lib/application.h"

Application::Application() {
  job_id_ = 1;
  data_id_ = 1;
}

void Application::load() {
  std::cout << "Loaded Nimbus base application." << std::endl;
}

void Application::start(SchedulerClient* sc) {
  std::cout << "Running Nimbus application: " << id_ << std::endl;
  client_ = sc;
  load();
}

void Application::registerJob(std::string name, Job* j) {
  job_table_[name] = j;
}

void Application::registerData(std::string name, Data* d) {
  data_table_[name] = d;
}

void Application::spawnJob(std::string name, int id, IDSet before, IDSet after,
        IDSet read, IDSet write, std::string params) {
  IDSet idset;
  idset.insert(id);
  std::string str = "spawnjob";
  str = str + " name:" + name;
  str = str + " id:" + idset.toString();
  str = str + " before:" + before.toString();
  str = str + " after:" + after.toString();
  str = str + " read:" + read.toString();
  str = str + " write:" + write.toString();
  str = str + " param:" + params;
  Command cm(str);
  client_->sendCommand(&cm);
}

void Application::defineData(std::string name, int id) {
  IDSet idset;
  idset.insert(id);
  std::string str = "definedata";
  str = str + " name:" + name;
  str = str + " id:" + idset.toString();
  Command cm(str);
  client_->sendCommand(&cm);
}

Job* Application::cloneJob(std::string name) {
  return job_table_[name]->clone();
}

Data* Application::cloneData(std::string name) {
  return data_table_[name]->clone();
}

void Application::getNewJobID(int req_num, std::vector<int>* result) {
  for (int i = 0; i < req_num; i++) {
    result->push_back(job_id_);
    job_id_++;
  }
}

void Application::getNewDataID(int req_num, std::vector<int>* result) {
  for (int i = 0; i < req_num; i++) {
    result->push_back(data_id_);
    data_id_++;
  }
}




