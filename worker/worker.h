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

#ifndef MUTE_LOG
#define MUTE_LOG
#endif  // MUTE_LOG

#include <boost/thread.hpp>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include "worker/data.h"
#include "worker/job.h"
#include "worker/application.h"
#include "worker/physical_data_map.h"
#include "worker/worker_job_graph/worker_job_graph.h"
#include "shared/nimbus_types.h"
#include "shared/id_maker.h"
#include "shared/scheduler_client.h"
#include "shared/scheduler_command_include.h"
#include "shared/worker_data_exchanger.h"
#include "shared/cluster.h"
#include "shared/parser.h"
#include "shared/log.h"
#include "shared/high_resolution_timer.h"
#include "shared/profiler.h"

namespace nimbus {

#ifndef DBG_WORKER_FD_S
#define DBG_WORKER_FD_S "\n[WORKER_FD] "
#endif  // DBG_WORKER_FD_S

#ifndef DBG_WORKER_BD_S
#define DBG_WORKER_BD_S "\n[WORKER_BD] "
#endif  // DBG_WORKER_BD_S

class Worker;
class WorkerManager;
typedef std::map<int, Worker*> WorkerMap;
class WorkerThreadMonitor;

class Worker {
  friend class WorkerThreadMonitor;
 public:
  Worker(std::string scheuler_ip, port_t scheduler_port,
      port_t listening_port_, Application* application);
  virtual ~Worker();

  virtual void Run();
  virtual void WorkerCoreProcessor();

  // virtual void ExecuteJob(Job* job);
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
  void set_ip_address(std::string ip);
  virtual PhysicalDataMap* data_map();

  // TODO(quhang) maybe not a good interface.
  void SendCommand(SchedulerCommand* command) {
    client_->SendCommand(command);
  }

  Log *cache_log;

 protected:
  SchedulerClient* client_;
  WorkerDataExchanger* data_exchanger_;
  WorkerLdoMap* ldo_map_;
  IDMaker id_maker_;
  SchedulerCommand::PrototypeTable scheduler_command_table_;
  worker_id_t id_;
  std::string ip_address_;
  std::string scheduler_ip_;
  port_t scheduler_port_;
  port_t listening_port_;

 private:
  FILE* event_log;
  WorkerJobGraph worker_job_graph_;
  Log log_;
  Log version_log_;
  Log data_hash_log_;
  Computer host_;
  boost::thread* client_thread_;
  boost::thread* data_exchanger_thread_;
  boost::thread* profiler_thread_;
  // TODO(quhang) a strong assumption is made that the data map is never changed
  // during the runtime. Indeed, for now, it is only changed at the very
  // beginning of the simulation, so it keeps the same during the simulation,
  // which might break in the future.

  PhysicalDataMap data_map_;
  JobList ready_jobs_;
  Application* application_;
  job_id_t DUMB_JOB_ID;
  WorkerManager* worker_manager_;
  HighResolutionTimer timer_;
  Profiler profiler_;

  virtual void SetupSchedulerInterface();

  virtual void SetupDataExchangerInterface();

  virtual void LoadSchedulerCommands();

  virtual void AddJobToGraph(Job* job);
  virtual void NotifyLocalJobDone(Job* job);
  virtual void NotifyJobDone(job_id_t job_id);
  virtual void ClearAfterSet(WorkerJobVertex* vertex);
  virtual void NotifyTransmissionDone(job_id_t job_id);

  void PrintTimeStamp(const char* format, ...);

 public:
  void ResolveDataArray(Job* job);
};
}  // namespace nimbus

#endif  // NIMBUS_WORKER_WORKER_H_
