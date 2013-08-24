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
  * Nimbus abstraction of an application. Programmers use this base class to
  * write various application served by Nimbus.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#ifndef NIMBUS_WORKER_APPLICATION_H_
#define NIMBUS_WORKER_APPLICATION_H_

#include <map>
#include <string>
#include <vector>
#include "worker/job.h"
#include "worker/data.h"
#include "shared/nimbus_types.h"
#include "shared/scheduler_client.h"
#include "shared/scheduler_command.h"

using namespace nimbus; // NOLINT

namespace nimbus {

class Application;
typedef std::map<int, Application*> AppMap;

class Application {
 public:
  Application();
  ~Application() {}

  virtual void load();
  virtual void start(SchedulerClient* scheduler);

  void registerJob(std::string name, Job* job);
  void registerData(std::string name, Data* data);
  void SpawnJob(const std::string& name,
      const job_id_t& id,
      const IDSet<data_id_t>& read,
      const IDSet<data_id_t>& write,
      const IDSet<job_id_t>& before,
      const IDSet<job_id_t>& after,
      const JobType& type,
      std::string params);

  void DefineData(const std::string& name,
      const data_id_t& id,
      const partition_t& partition_id,
      const IDSet<partition_t>& neighbor_partition,
      std::string params);

  Job* cloneJob(std::string name);
  Data* cloneData(std::string name);
  void getNewJobID(int req_num, std::vector<int>* result);
  void getNewDataID(int req_num, std::vector<int>* result);
  void* app_data();
  void set_app_data(void* data);

 private:
  uint64_t id_;
  uint64_t priority_;
  uint64_t job_id_;
  uint64_t data_id_;

  IDSet<job_id_t> job_ID_pool_;
  JobTable job_table_;
  DataTable data_table_;
  SchedulerClient* client_;
  void* app_data_;
};

}  //  namespace nimbus

#endif  // NIMBUS_WORKER_APPLICATION_H_
