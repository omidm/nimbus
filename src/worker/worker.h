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

#ifndef NIMBUS_SRC_WORKER_WORKER_H_
#define NIMBUS_SRC_WORKER_WORKER_H_

#ifndef MUTE_LOG
#define MUTE_LOG
#endif  // MUTE_LOG

#include <boost/thread.hpp>
#include <boost/unordered_set.hpp>
#include <cstdio>
#include <list>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "src/worker/data.h"
#include "src/worker/job.h"
#include "src/worker/application.h"
#include "src/worker/physical_data_map.h"
#include "src/worker/worker_job_graph/worker_job_graph.h"
#include "src/shared/nimbus_types.h"
#include "src/shared/id_maker.h"
#include "src/shared/scheduler_client.h"
#include "src/shared/scheduler_command_include.h"
#include "src/shared/worker_data_exchanger.h"
#include "src/shared/distributed_db.h"
#include "src/shared/cluster.h"
#include "src/shared/parser.h"
#include "src/shared/high_resolution_timer.h"
#include "src/shared/multi_level_timer.h"
#include "src/shared/profiler.h"
#include "src/shared/execution_template.h"
#include "src/shared/helpers.h"

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
  virtual void CreateModules();
  virtual void SetupTimers();
  virtual void SetupApplication();
  virtual void SetupWorkerManager();
  virtual void SetupSchedulerClient();
  virtual void SetupWorkerDataExchanger();

  virtual void SetupCommandProcessor();
  virtual void RunCommandProcessor();

  virtual void SetupReceiveEventProcessor();
  virtual void RunReceiveEventProcessor();

  virtual void WorkerCoreProcessor();

  // virtual void ExecuteJob(Job* job);
  virtual void ProcessSchedulerCommand(SchedulerCommand* command);
  virtual void ProcessComputeJobCommand(ComputeJobCommand* command);
  virtual void ProcessCreateDataCommand(CreateDataCommand* command);
  virtual void ProcessRemoteCopySendCommand(RemoteCopySendCommand* command);
  virtual void ProcessRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* command);
  virtual void ProcessMegaRCRCommand(MegaRCRCommand* command);
  virtual void ProcessLocalCopyCommand(LocalCopyCommand* command);
  virtual void ProcessHandshakeCommand(HandshakeCommand* command);
  virtual void ProcessJobDoneCommand(JobDoneCommand* command);
  virtual void ProcessLdoAddCommand(LdoAddCommand* command);
  virtual void ProcessLdoRemoveCommand(LdoRemoveCommand* command);
  virtual void ProcessPartitionAddCommand(PartitionAddCommand* command);
  virtual void ProcessPartitionRemoveCommand(PartitionRemoveCommand* command);
  virtual void ProcessTerminateCommand(TerminateCommand* command);
  virtual void ProcessDefinedTemplateCommand(DefinedTemplateCommand* command);
  virtual void ProcessSaveDataCommand(SaveDataCommand* command);
  virtual void ProcessLoadDataCommand(LoadDataCommand* command);
  virtual void ProcessPrepareRewindCommand(PrepareRewindCommand* command);
  virtual void ProcessRequestStatCommand(RequestStatCommand *command);
  virtual void ProcessPrintStatCommand(PrintStatCommand *command);
  virtual void ProcessStartCommandTemplateCommand(StartCommandTemplateCommand* command);
  virtual void ProcessEndCommandTemplateCommand(EndCommandTemplateCommand* command);
  virtual void ProcessSpawnCommandTemplateCommand(SpawnCommandTemplateCommand* command);

  virtual void NotifyLocalJobDone(Job* job);

  worker_id_t id();
  void set_id(worker_id_t id);
  void set_ip_address(std::string ip);
  void set_execution_template_active(bool flag);
  virtual PhysicalDataMap* data_map();

  // TODO(quhang) maybe not a good interface.
  void SendCommand(SchedulerCommand* command) {
    client_->SendCommand(command);
  }

 protected:
  SchedulerClient* client_;
  WorkerDataExchanger* data_exchanger_;
  DistributedDB *ddb_;
  WorkerLdoMap* ldo_map_;
  IDMaker *id_maker_;
  SchedulerCommand::PrototypeTable scheduler_command_table_;
  worker_id_t id_;
  std::string ip_address_;
  std::string scheduler_ip_;
  port_t scheduler_port_;
  port_t listening_port_;
  bool execution_template_active_;

 private:
  static const uint64_t max_hint_size_ = 10000;
  boost::unordered_set<job_id_t> hint_set_;
  std::list<job_id_t> hint_queue_;
  void AddFinishHintSet(const job_id_t job_id);
  bool InFinishHintSet(const job_id_t job_id);
  WorkerJobGraph worker_job_graph_;
  boost::recursive_mutex job_graph_mutex_;
  boost::condition_variable_any job_graph_cond_;

  boost::recursive_mutex *receive_event_mutex_;
  boost::condition_variable_any *receive_event_cond_;

  boost::recursive_mutex *command_processor_mutex_;
  boost::condition_variable_any *command_processor_cond_;

  Computer host_;
  boost::thread* client_thread_;
  boost::thread* profiler_thread_;
  boost::thread* command_processor_thread_;
  boost::thread* receive_event_thread_;
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

  typedef std::map<template_id_t, WorkerDataExchanger::EventList> EventMap;

  bool filling_execution_template_;
  std::string execution_template_in_progress_;
  std::map<std::string, ExecutionTemplate*> execution_templates_;
  std::map<template_id_t, ExecutionTemplate*> active_execution_templates_;
  EventMap pending_events_;

  bool prepare_rewind_phase_;

  virtual void LoadSchedulerCommands();

  virtual void AddJobToGraph(Job* job);
  virtual void NotifyJobDone(job_id_t job_id, bool final);
  virtual void ClearAfterSet(WorkerJobVertex* vertex);
  virtual void ProcessReceiveEvents(const WorkerDataExchanger::EventList& events);
  virtual void ProcessRCREvent(const WorkerDataExchanger::Event& event);
  virtual void ProcessMegaRCREvent(const WorkerDataExchanger::Event& event);

  virtual void SendJobDoneAndDeleteJob(Job* job);

  virtual void ClearBlockedJobs();
  virtual bool AllReadyJobsAreDone();

 public:
  void StatAddJob(size_t num);
  void StatEndJob(size_t num);
  void StatDispatchJob(size_t num);
  void ResolveDataArray(Job* job);
  void GetTimerStat(int64_t* idle, int64_t* block, int64_t* run);
  void PrintTimerStat();

 private:
  size_t stat_blocked_job_num_, stat_ready_job_num_;
  size_t stat_busy_cores_, stat_blocked_cores_, stat_idle_cores_;
  MultiLevelTimer run_timer_;
  MultiLevelTimer block_timer_;
  MultiLevelTimer total_timer_;
  boost::recursive_mutex stat_mutex_;
};
}  // namespace nimbus

#endif  // NIMBUS_SRC_WORKER_WORKER_H_
