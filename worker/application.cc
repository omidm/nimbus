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
  app_data_ = NULL;
}

void Application::Load() {
  std::cout << "Loaded Nimbus base application." << std::endl;
}

void Application::Start(SchedulerClient* client,
    IDMaker* id_maker,
    WorkerLdoMap* ldo_map) {
  std::cout << "Running Nimbus application: " << id_ << std::endl;
  client_ = client;
  id_maker_ = id_maker;
  ldo_map_ = ldo_map;
  Load();
}

void Application::RegisterJob(std::string name, Job* j) {
  job_table_[name] = j;
}

void Application::RegisterData(std::string name, Data* d) {
  data_table_[name] = d;
}

void Application::SpawnComputeJob(const std::string& name,
    const job_id_t& id,
    const IDSet<logical_data_id_t>& read,
    const IDSet<logical_data_id_t>& write,
    const IDSet<job_id_t>& before,
    const IDSet<job_id_t>& after,
    const job_id_t& parent_id,
    const Parameter& params) {

  SpawnComputeJobCommand cm(name, ID<job_id_t>(id), read, write, before, after,
      ID<job_id_t>(parent_id), params);
  client_->sendCommand(&cm);
}

void Application::SpawnCopyJob(const job_id_t& id,
    const logical_data_id_t& from_logical_id,
    const logical_data_id_t& to_logical_id,
    const IDSet<job_id_t>& before,
    const IDSet<job_id_t>& after,
    const job_id_t& parent_id,
    const Parameter& params) {

  SpawnCopyJobCommand cm(ID<job_id_t>(id), ID<logical_data_id_t>(from_logical_id),
      ID<logical_data_id_t>(to_logical_id), before, after, ID<job_id_t>(parent_id), params);
  client_->sendCommand(&cm);
}

void Application::DefineData(const std::string& name,
    const logical_data_id_t& logical_data_id,
    const partition_id_t& partition_id,
    const IDSet<partition_id_t>& neighbor_partitions,
    const job_id_t& parent_id,
    const Parameter& params) {
  ID<logical_data_id_t> logical_id_made(logical_data_id);
  ID<partition_id_t> partition_id_made(partition_id);
  ID<job_id_t> parent_id_made(parent_id);

  DefineDataCommand cm(name, logical_id_made, partition_id_made,
      neighbor_partitions, parent_id_made, params);
  client_->sendCommand(&cm);
}

void Application::DefinePartition(const ID<partition_id_t>& partition_id,
     const GeometricRegion& r,
     const Parameter& params) {
  DefinePartitionCommand cm(partition_id, r, params);
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

bool Application::GetNewLogicalDataID(std::vector<logical_data_id_t>* result, size_t req_num) {
  return id_maker_->GetNewLogicalDataID(result, req_num);
}

const LogicalDataObject* Application::GetLogicalObject(logical_data_id_t id) {
  if (ldo_map_ == NULL) {
    std::cout << "Error: GetLogicalObject, ldo_map_ has not been set." << std::endl;
    exit(-1);
  } else {
    return ldo_map_->FindLogicalObject(id);
  }
}

int Application::GetCoveredLogicalObjects(CLdoVector* result,
     std::string& variable,
     GeometricRegion* r) {
  if (ldo_map_ == NULL) {
    return false;
  } else {
    return ldo_map_->FindCoveredLogicalObjects(variable, r, result);
  }
}

int Application::GetAdjacentLogicalObjects(CLdoVector* result,
                                           std::string& variable,
                                           GeometricRegion* r) {
  if (ldo_map_ == NULL) {
    return false;
  } else {
    return ldo_map_->FindAdjacentLogicalObjects(variable, r, result);
  }
}

int Application::GetIntersectingLogicalObjects(CLdoVector* result,
                                               std::string& variable,
                                               GeometricRegion* r) {
  if (ldo_map_ == NULL) {
    return false;
  } else {
    return ldo_map_->FindIntersectingLogicalObjects(variable, r, result);
  }
}

void* Application::app_data() {
    return app_data_;
}

void Application::set_app_data(void* data) {
    app_data_ = data;
}
