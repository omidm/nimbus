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

#ifndef NIMBUS_WORKER_WORKER_H_
#define NIMBUS_WORKER_WORKER_H_

#include <boost/thread.hpp>
#include <string>
#include <vector>
#include <map>
#include "worker/data.h"
#include "worker/job.h"
#include "worker/application.h"
#include "shared/nimbus_types.h"
#include "shared/id_maker.h"
#include "shared/scheduler_client.h"
#include "shared/scheduler_command_include.h"
#include "shared/worker_data_exchanger.h"
#include "shared/cluster.h"
#include "shared/parser.h"
#include "shared/log.h"

namespace nimbus {

class Worker;
typedef std::map<int, Worker*> WorkerMap;

class Worker {
 public:
  Worker(std::string scheuler_ip, port_t scheduler_port,
      port_t listening_port_, Application* application);

  virtual void Run();
  virtual void WorkerCoreProcessor();
  virtual void ScanBlockedJobs();
  virtual void ScanPendingTransferJobs();
  virtual void GetJobsToRun(JobList* list, size_t max_num);
  virtual void ExecuteJob(Job* job);
  virtual void ProcessSchedulerCommand(SchedulerCommand* command);
  virtual void ProcessComputeJobCommand(ComputeJobCommand* command);
  virtual void ProcessCreateDataCommand(CreateDataCommand* command);
  virtual void ProcessRemoteCopySendCommand(RemoteCopySendCommand* command);
  virtual void ProcessRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* command);
  virtual void ProcessLocalCopyCommand(LocalCopyCommand* command);
  virtual void ProcessHandshakeCommand(HandshakeCommand* command);
  virtual void ProcessJobDoneCommand(JobDoneCommand* command);
  virtual void ProcessLdoAddCommand(LdoAddCommand* command);
  virtual void ProcessLdoRemoveCommand(LdoRemoveCommand* command);
  virtual void ProcessPartitionAddCommand(PartitionAddCommand* command);
  virtual void ProcessPartitionRemoveCommand(PartitionRemoveCommand* command);
  virtual void ProcessTerminateCommand(TerminateCommand* command);

  worker_id_t id();
  void set_id(worker_id_t id);
  virtual PhysicalDataMap* data_map();

 protected:
  SchedulerClient* client_;
  WorkerDataExchanger* data_exchanger_;
  WorkerLdoMap* ldo_map_;
  IDMaker id_maker_;
  SchedulerCommand::PrototypeTable scheduler_command_table_;
  worker_id_t id_;
  std::string scheduler_ip_;
  port_t scheduler_port_;
  port_t listening_port_;

 private:
  Log log_;
  Computer host_;
  boost::thread* client_thread_;
  boost::thread* data_exchanger_thread_;
  PhysicalDataMap data_map_;
  JobList ready_jobs_;
  JobList blocked_jobs_;
  JobList pending_transfer_jobs_;
  // std::map<job_id_t, IDSet<job_id_t> > done_jobs_;
  Application* application_;

  virtual void SetupSchedulerInterface();

  virtual void SetupDataExchangerInterface();

  virtual void AddData(Data* data);
  virtual void DeleteData(physical_data_id_t physical_data_id);
  virtual void LoadSchedulerCommands();
};

}  // namespace nimbus

#endif  // NIMBUS_WORKER_WORKER_H_
