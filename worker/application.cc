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

#include "worker/application.h"

Application::Application() {
  job_id_ = 1;
  data_id_ = 1;
  app_data_ = NULL;
}

void Application::Load() {
  std::cout << "Loaded Nimbus base application." << std::endl;
}

void Application::Start(SchedulerClient* client, IDMaker* id_maker) {
  std::cout << "Running Nimbus application: " << id_ << std::endl;
  client_ = client;
  id_maker_ = id_maker;
  Load();
}

void Application::RegisterJob(std::string name, Job* j) {
  job_table_[name] = j;
}

void Application::RegisterData(std::string name, Data* d) {
  data_table_[name] = d;
}

void Application::SpawnJob(const std::string& name,
    const job_id_t& id,
    const IDSet<data_id_t>& read,
    const IDSet<data_id_t>& write,
    const IDSet<job_id_t>& before,
    const IDSet<job_id_t>& after,
    const JobType& type,
    const Parameter& params) {
  IDSet<job_id_t> id_set;
  id_set.insert(id);

  SpawnJobCommand cm(name, id_set, read, write, before, after, type, params);
  client_->sendCommand(&cm);
}

void Application::SpawnComputeJob(const std::string& name,
    const job_id_t& id,
    const IDSet<data_id_t>& read,
    const IDSet<data_id_t>& write,
    const IDSet<job_id_t>& before,
    const IDSet<job_id_t>& after,
    const Parameter& params) {

  SpawnComputeJobCommand cm(name, ID<job_id_t>(id), read, write, before, after, params);
  client_->sendCommand(&cm);
}

void Application::SpawnCopyJob(const job_id_t& id,
    const data_id_t& from_id,
    const data_id_t& to_id,
    const IDSet<job_id_t>& before,
    const IDSet<job_id_t>& after,
    const Parameter& params) {

  SpawnCopyJobCommand cm(ID<job_id_t>(id), ID<data_id_t>(from_id),
      ID<data_id_t>(to_id), before, after, params);
  client_->sendCommand(&cm);
}

void Application::DefineData(const std::string& name,
    const data_id_t& id,
    const partition_t& partition_id,
    const IDSet<partition_t>& neighbor_partitions,
    const Parameter& params) {
  ID<data_id_t> id_made(id);
  ID<partition_t> partition_id_made(partition_id);

  DefineDataCommand cm(name, id_made, partition_id_made, neighbor_partitions, params);
  client_->sendCommand(&cm);
}

Job* Application::CloneJob(std::string name) {
  return job_table_[name]->Clone();
}

Data* Application::CloneData(std::string name) {
  return data_table_[name]->Clone();
}

bool Application::GetNewJobID(std::vector<job_id_t>* result, size_t req_num) {
  return id_maker_->GetNewJobID(result, req_num);
}

bool Application::GetNewDataID(std::vector<data_id_t>* result, size_t req_num) {
  return id_maker_->GetNewDataID(result, req_num);
}


// TODO(omidm) remove both of them.
void Application::GetNewJobID(int req_num, std::vector<int>* result) {
  for (int i = 0; i < req_num; i++) {
    result->push_back(job_id_);
    job_id_++;
  }
}

void Application::GetNewDataID(int req_num, std::vector<int>* result) {
  for (int i = 0; i < req_num; i++) {
    result->push_back(data_id_);
    data_id_++;
  }
}



void* Application::app_data() {
    return app_data_;
}

void Application::set_app_data(void* data) {
    app_data_ = data;
}
