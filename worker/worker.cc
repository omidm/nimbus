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

#include <unistd.h>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <string>
#include <ctime>
#include <list>
#include <limits>
#include "shared/fast_log.hh"
#include "shared/profiler_malloc.h"
#include "worker/worker.h"
#include "worker/worker_ldo_map.h"
#include "worker/worker_manager.h"
#include "worker/util_dumping.h"
#include "data/physbam/physbam_data.h"

#define SCHEDULER_COMMAND_GROUP_QUOTA 10

using boost::hash;

namespace nimbus {

Worker::Worker(std::string scheduler_ip, port_t scheduler_port,
    port_t listening_port, Application* a)
: scheduler_ip_(scheduler_ip),
  scheduler_port_(scheduler_port),
  listening_port_(listening_port),
  application_(a) {
    {
      struct sched_param param;
      param.sched_priority = 0;
      int st = pthread_setschedparam(pthread_self(), SCHED_OTHER, &param);
      if (st != 0) {
        // Scheduling setting goes wrong.
        exit(1);
      }
    }
    id_ = -1;
    ip_address_ = NIMBUS_RECEIVER_KNOWN_IP;
    worker_manager_ = new WorkerManager();
    ddb_ = new DistributedDB();
    DUMB_JOB_ID = std::numeric_limits<job_id_t>::max();
    worker_job_graph_.AddVertex(
        DUMB_JOB_ID,
        new WorkerJobEntry(
            DUMB_JOB_ID, NULL, WorkerJobEntry::CONTROL));
}

Worker::~Worker() {
  WorkerJobVertex* vertex = NULL;
  worker_job_graph_.GetVertex(DUMB_JOB_ID, &vertex);
  delete vertex->entry();
  worker_job_graph_.RemoveVertex(DUMB_JOB_ID);
  delete worker_manager_;
  delete ddb_;
}

void Worker::Run() {
  std::cout << "Running the Worker" << std::endl;

  SetupSchedulerInterface();

  ldo_map_ = new WorkerLdoMap();
  application_->Start(client_, &id_maker_, ldo_map_);

  SetupDataExchangerInterface();

  WorkerCoreProcessor();
}

void Worker::WorkerCoreProcessor() {
  timer::InitializeKeys();
  timer::InitializeTimers();
  // Since this timer is used for idle time computation, it should start after
  // the initial data mep generation trigered by main job is over. -omidm
  // timer::StartTimer(timer::kSumCyclesTotal,
  //                   WorkerManager::across_job_parallism);
  stat_blocked_job_num_ = 0;
  stat_ready_job_num_ = 0;
  stat_busy_cores_ = 0;
  stat_blocked_cores_ = 0;
  stat_idle_cores_ = WorkerManager::across_job_parallism;
  run_timer_.set_name("kSumCyclesRun");
  block_timer_.set_name("kSumCyclesBlock");
  total_timer_.set_name("kSumCyclesTotal");
  std::cout << "Base Worker Core Processor" << std::endl;
  worker_manager_->worker_ = this;
  dbg(DBG_WORKER_FD, DBG_WORKER_FD_S"Launching worker threads.\n");
  worker_manager_->StartWorkerThreads();
  dbg(DBG_WORKER_FD, DBG_WORKER_FD_S"Finishes launching worker threads.\n");
  worker_manager_->TriggerScheduling();

  JobList local_job_done_list;
  while (true) {
    bool process_jobs = false;
    // Process command.
    SchedulerCommandList storage;
    client_->ReceiveCommands(&storage, SCHEDULER_COMMAND_GROUP_QUOTA);
    SchedulerCommandList::iterator iter = storage.begin();
    for (; iter != storage.end(); ++iter) {
      timer::StartTimer(timer::kCoreCommand);
      SchedulerCommand *comm = *iter;
      dbg(DBG_WORKER, "Receives command from scheduler: %s\n",
          comm->ToString().c_str());
      dbg(DBG_WORKER_FD,
          DBG_WORKER_FD_S"Scheduler command arrives(%s).\n",
          comm->name().c_str());
      process_jobs = true;
      ProcessSchedulerCommand(comm);
      delete comm;
      timer::StopTimer(timer::kCoreCommand);
    }
    // Poll jobs that finish receiving.
    job_id_t receive_job_id;
    int quota = 10;
    while (data_exchanger_->GetReceiveEvent(&receive_job_id)) {
      timer::StartTimer(timer::kCoreTransmission);
      dbg(DBG_WORKER_FD,
          DBG_WORKER_FD_S"Receive-job transmission is done(job #%d)\n",
          receive_job_id);
      process_jobs = true;
      NotifyTransmissionDone(receive_job_id);
      timer::StopTimer(timer::kCoreTransmission);
      if (--quota <= 0) {
        break;
      }
    }
    // Job done processing.
    // worker_manager_->GetLocalJobDoneList(&local_job_done_list);
    // quota = 10;
    // while (!local_job_done_list.empty()) {
    //   timer::StartTimer(timer::kCoreJobDone);
    //   Job* job = local_job_done_list.front();
    //   local_job_done_list.pop_front();
    //   process_jobs = true;
    //   NotifyLocalJobDone(job);
    //   timer::StopTimer(timer::kCoreJobDone);
    //   if (--quota <= 0) {
    //     break;
    //   }
    // }

    if (!process_jobs) {
//      typename WorkerJobVertex::Iter iter = worker_job_graph_.begin();
//      for (; iter != worker_job_graph_.end(); ++iter) {
//        if (iter->second->incoming_edges()->size() != 0) {
//          Job* job = iter->second->entry()->get_job();
//          std::string name = job->name();
//          if (name.find("Copy") == std::string::npos &&
//              name.find("extrapolate_phi") != std::string::npos) {
//              std::cout << "OMID: " << job->id().elem() << " "
//                << iter->second->incoming_edges()->size() << " ";
//            WorkerJobEdge::Map::iterator it = iter->second->incoming_edges()->begin();
//            for (; it != iter->second->incoming_edges()->end(); ++it) {
//              WorkerJobEdge* edge = it->second;
//              std::cout << edge->start_vertex()->entry()->get_job_id() << " ";
//            }
//            std::cout << job->before_set().ToNetworkData() << std::endl;
//          }
//        }
//      }

      usleep(10);
    }
  }
}

// Extracts data objects from the read/write set to data array.
void Worker::ResolveDataArray(Job* job) {
  dbg(DBG_WORKER_FD, DBG_WORKER_FD_S"Job(name %s, #%d) ready to run.\n",
      job->name().c_str(), job->id().elem());
  if ((dynamic_cast<CreateDataJob*>(job) != NULL)) {  // NOLINT
    assert(job->read_set().size() == 0);
    assert(job->write_set().size() == 1);
    job->data_array.clear();
    job->data_array.push_back(
        data_map_.AcquireAccess(*job->write_set().begin(), job->id().elem(),
                                PhysicalDataMap::INIT));
    return;
  }
  job->data_array.clear();
  IDSet<physical_data_id_t>::IDSetIter iter;

  const IDSet<physical_data_id_t>& read = job->get_read_set();
  for (iter = read.begin(); iter != read.end(); iter++) {
    job->data_array.push_back(
        data_map_.AcquireAccess(*iter, job->id().elem(),
                                PhysicalDataMap::READ));
  }
  const IDSet<physical_data_id_t>& write = job->get_write_set();
  for (iter = write.begin(); iter != write.end(); iter++) {
    job->data_array.push_back(
        data_map_.AcquireAccess(*iter, job->id().elem(),
                                PhysicalDataMap::WRITE));
  }
}

void Worker::ProcessSchedulerCommand(SchedulerCommand* cm) {
  switch (cm->type()) {
    case SchedulerCommand::HANDSHAKE:
      ProcessHandshakeCommand(reinterpret_cast<HandshakeCommand*>(cm));
      break;
    case SchedulerCommand::JOB_DONE:
      ProcessJobDoneCommand(reinterpret_cast<JobDoneCommand*>(cm));
      break;
    case SchedulerCommand::EXECUTE_COMPUTE:
      ProcessComputeJobCommand(reinterpret_cast<ComputeJobCommand*>(cm));
      break;
    case SchedulerCommand::CREATE_DATA:
      ProcessCreateDataCommand(reinterpret_cast<CreateDataCommand*>(cm));
      break;
    case SchedulerCommand::REMOTE_SEND:
      ProcessRemoteCopySendCommand(reinterpret_cast<RemoteCopySendCommand*>(cm));
      break;
    case SchedulerCommand::REMOTE_RECEIVE:
      ProcessRemoteCopyReceiveCommand(reinterpret_cast<RemoteCopyReceiveCommand*>(cm));
      break;
    case SchedulerCommand::LOCAL_COPY:
      ProcessLocalCopyCommand(reinterpret_cast<LocalCopyCommand*>(cm));
      break;
    case SchedulerCommand::LDO_ADD:
      ProcessLdoAddCommand(reinterpret_cast<LdoAddCommand*>(cm));
      break;
    case SchedulerCommand::LDO_REMOVE:
      ProcessLdoRemoveCommand(reinterpret_cast<LdoRemoveCommand*>(cm));
      break;
    case SchedulerCommand::PARTITION_ADD:
      ProcessPartitionAddCommand(reinterpret_cast<PartitionAddCommand*>(cm));
      break;
    case SchedulerCommand::PARTITION_REMOVE:
      ProcessPartitionRemoveCommand(reinterpret_cast<PartitionRemoveCommand*>(cm));
      break;
    case SchedulerCommand::TERMINATE:
      ProcessTerminateCommand(reinterpret_cast<TerminateCommand*>(cm));
      break;
    case SchedulerCommand::DEFINED_TEMPLATE:
      ProcessDefinedTemplateCommand(reinterpret_cast<DefinedTemplateCommand*>(cm));
      break;
    case SchedulerCommand::SAVE_DATA:
      ProcessSaveDataCommand(reinterpret_cast<SaveDataCommand*>(cm));
      break;
    case SchedulerCommand::LOAD_DATA:
      ProcessLoadDataCommand(reinterpret_cast<LoadDataCommand*>(cm));
      break;
    case SchedulerCommand::PREPARE_REWIND:
      ProcessPrepareRewindCommand(reinterpret_cast<PrepareRewindCommand*>(cm));
      break;
    case SchedulerCommand::REQUEST_STAT:
      ProcessRequestStatCommand(reinterpret_cast<RequestStatCommand*>(cm));
      break;
    default:
      std::cout << "ERROR: " << cm->ToNetworkData() <<
        " have not been implemented in ProcessSchedulerCommand yet." <<
        std::endl;
      exit(-1);
  }
}

// Processes handshake command. Configures the worker based on the handshake
// command and responds by sending another handshake command.
void Worker::ProcessHandshakeCommand(HandshakeCommand* cm) {
  double time = Log::GetRawTime();
  ID<port_t> port(listening_port_);
  HandshakeCommand new_cm = HandshakeCommand(cm->worker_id(),
      // boost::asio::ip::host_name(), port);
      // "127.0.0.1", port);
      ip_address_, port, time);
  client_->SendCommand(&new_cm);

  if (ip_address_ == NIMBUS_RECEIVER_KNOWN_IP) {
    ip_address_ = cm->ip();
  }
  id_ = cm->worker_id().elem();
  id_maker_.Initialize(id_);
  ddb_->Initialize(ip_address_, id_);

  std::string wstr = int2string(id_);
}

// Processes jobdone command. Moves a job from blocked queue to ready queue if
// its before set is satisfied.
void Worker::ProcessJobDoneCommand(JobDoneCommand* cm) {
  NotifyJobDone(cm->job_id().elem(), cm->final());
  // If job main is over the data map creation phase is over, now start the
  // timer that goes toward the idle time computation.  -omidm
  if (cm->job_id().elem() == NIMBUS_KERNEL_JOB_ID + 1) {
    // timer::StartTimer(timer::kSumCyclesTotal, WorkerManager::across_job_parallism);
    total_timer_.Start(WorkerManager::across_job_parallism);
  }
}

// Processes computejob command. Generates the corresponding job and pushes it
// to the blocking queue.
void Worker::ProcessComputeJobCommand(ComputeJobCommand* cm) {
  Job* job = application_->CloneJob(cm->job_name());
  job->set_name("Compute:" + cm->job_name());
  job->set_id(cm->job_id());
  job->set_read_set(cm->read_set());
  job->set_write_set(cm->write_set());
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  job->set_future_job_id(cm->future_job_id());
  job->set_sterile(cm->sterile());
  job->set_region(cm->region());
  job->set_parameters(cm->params());
  AddJobToGraph(job);
}

// Processes createdata command. Generates the corresponding data and pushes a
// data creation job to the blocking queue.
void Worker::ProcessCreateDataCommand(CreateDataCommand* cm) {
  Data * data = application_->CloneData(cm->data_name());
  data->set_logical_id(cm->logical_data_id().elem());
  data->set_physical_id(cm->physical_data_id().elem());
  // data->set_name(cm->data_name());
  const LogicalDataObject* ldo;
  ldo = ldo_map_->FindLogicalObject(cm->logical_data_id().elem());
  data->set_region(*(ldo->region()));

  boost::unique_lock<boost::recursive_mutex> lock(job_graph_mutex_);
  data_map_.AddMapping(data->physical_id(), data);

  Job * job = new CreateDataJob();
  job->set_name("CreateData:" + cm->data_name());
  job->set_id(cm->job_id());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  AddJobToGraph(job);
}

void Worker::ProcessRemoteCopySendCommand(RemoteCopySendCommand* cm) {
  RemoteCopySendJob * job = new RemoteCopySendJob(data_exchanger_, application_);
  data_exchanger_->AddContactInfo(cm->to_worker_id().elem(),
      cm->to_ip(), cm->to_port().elem());
  job->set_name("RemoteCopySend");
  job->set_id(cm->job_id());
  job->set_receive_job_id(cm->receive_job_id());
  job->set_to_worker_id(cm->to_worker_id());
  job->set_to_ip(cm->to_ip());
  job->set_to_port(cm->to_port());
  IDSet<physical_data_id_t> read_set;
  read_set.insert(cm->from_physical_data_id().elem());
  job->set_read_set(read_set);
  job->set_before_set(cm->before_set());
  AddJobToGraph(job);
}

void Worker::ProcessRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* cm) {
  Job * job = new RemoteCopyReceiveJob(application_);
  job->set_name("RemoteCopyReceive");
  job->set_id(cm->job_id());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->to_physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  AddJobToGraph(job);
}

void Worker::ProcessSaveDataCommand(SaveDataCommand* cm) {
  SaveDataJob * job = new SaveDataJob(ddb_, application_);
  job->set_name("SaveData");
  job->set_id(cm->job_id());
  job->set_checkpoint_id(cm->checkpoint_id().elem());
  IDSet<physical_data_id_t> read_set;
  read_set.insert(cm->from_physical_data_id().elem());
  job->set_read_set(read_set);
  job->set_before_set(cm->before_set());
  AddJobToGraph(job);
}

void Worker::ProcessLoadDataCommand(LoadDataCommand* cm) {
  LoadDataJob * job = new LoadDataJob(ddb_, application_);
  job->set_name("LoadData");
  job->set_id(cm->job_id());
  job->set_handle(cm->handle());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->to_physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  AddJobToGraph(job);
}

void Worker::ProcessPrepareRewindCommand(PrepareRewindCommand* cm) {
  boost::unique_lock<boost::recursive_mutex> lock(job_graph_mutex_);

  // First remove all blocked jobs.
  ClearBlockedJobs(&worker_job_graph_);

  // Wait untill all runing jobs finish.
  while (!IsEmptyGraph(&worker_job_graph_)) {
    JobList local_job_done_list;
    worker_manager_->GetLocalJobDoneList(&local_job_done_list);
    bool new_done = false;
    while (!local_job_done_list.empty()) {
      Job* job = local_job_done_list.front();
      local_job_done_list.pop_front();
      new_done = true;
      NotifyLocalJobDone(job);
    }
    if (!new_done) {
      usleep(10);
    }
  }

  PrepareRewindCommand command(ID<worker_id_t>(id_), cm->checkpoint_id());
  client_->SendCommand(&command);
}

void Worker::ProcessRequestStatCommand(RequestStatCommand* cm) {
  int64_t idle, block, run;
  GetTimerStat(&idle, &block, &run);
  RespondStatCommand command(cm->query_id(), id_, run, block, idle);
  client_->SendCommand(&command);
}


void Worker::ProcessLocalCopyCommand(LocalCopyCommand* cm) {
  Job * job = new LocalCopyJob(application_);
  job->set_name("LocalCopy");
  job->set_id(cm->job_id());
  IDSet<physical_data_id_t> read_set;
  read_set.insert(cm->from_physical_data_id().elem());
  job->set_read_set(read_set);
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->to_physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  AddJobToGraph(job);
}

void Worker::ProcessLdoAddCommand(LdoAddCommand* cm) {
    const LogicalDataObject* ldo = cm->object();
    if (!ldo_map_->AddLogicalObject(ldo->id(), ldo->variable(), *(ldo->region()) ) )
        dbg(DBG_ERROR, "Worker could not add logical object %i to ldo map\n", (ldo->id()));
}

void Worker::ProcessLdoRemoveCommand(LdoRemoveCommand* cm) {
    const LogicalDataObject* ldo = cm->object();
    if (!ldo_map_->RemoveLogicalObject(ldo->id()))
        dbg(DBG_ERROR, "Worker could not remove logical object %i to ldo map\n", (ldo->id()));
}

void Worker::ProcessPartitionAddCommand(PartitionAddCommand* cm) {
    GeometricRegion r = *(cm->region());
    if (!ldo_map_->AddPartition(cm->id().elem(), r))
        dbg(DBG_ERROR, "Worker could not add partition %i to ldo map\n", cm->id().elem());
}

void Worker::ProcessPartitionRemoveCommand(PartitionRemoveCommand* cm) {
  if (!ldo_map_->RemovePartition(cm->id().elem()))
    dbg(DBG_ERROR, "Worker could not remove partition %i from ldo map\n", cm->id().elem());
}

void Worker::ProcessTerminateCommand(TerminateCommand* cm) {
  // profiler_thread_->interrupt();
  // profiler_thread_->join();
  std::string file_name = int2string(id_) + "_time_per_thread.txt";
  FILE* temp = fopen(file_name.c_str(), "w");
  total_timer_.Print(temp);
  block_timer_.Print(temp);
  run_timer_.Print(temp);
  timer::PrintTimerSummary(temp);
  fclose(temp);
  exit(cm->exit_status().elem());
}

void Worker::ProcessDefinedTemplateCommand(DefinedTemplateCommand* cm) {
  application_->DefinedTemplate(cm->job_graph_name());
}


void Worker::SetupDataExchangerInterface() {
  data_exchanger_ = new WorkerDataExchanger(listening_port_);
  data_exchanger_thread_ = new boost::thread(
      boost::bind(&WorkerDataExchanger::Run, data_exchanger_));
}

void Worker::SetupSchedulerInterface() {
  LoadSchedulerCommands();
  client_ = new SchedulerClient(scheduler_ip_, scheduler_port_);
  client_->set_scheduler_command_table(&scheduler_command_table_);
  // client_->Run();
  client_thread_ = new boost::thread(
      boost::bind(&SchedulerClient::Run, client_));

  // profiler_thread_ = new boost::thread(
  //     boost::bind(&Profiler::Run, &profiler_, client_, id_));
}

void Worker::LoadSchedulerCommands() {
  // std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  scheduler_command_table_[SchedulerCommand::HANDSHAKE] = new HandshakeCommand();
  scheduler_command_table_[SchedulerCommand::JOB_DONE] = new JobDoneCommand();
  scheduler_command_table_[SchedulerCommand::EXECUTE_COMPUTE] = new ComputeJobCommand();
  scheduler_command_table_[SchedulerCommand::CREATE_DATA] = new CreateDataCommand();
  scheduler_command_table_[SchedulerCommand::REMOTE_SEND] = new RemoteCopySendCommand();
  scheduler_command_table_[SchedulerCommand::REMOTE_RECEIVE] = new RemoteCopyReceiveCommand(); // NOLINT
  scheduler_command_table_[SchedulerCommand::LOCAL_COPY] = new LocalCopyCommand();
  scheduler_command_table_[SchedulerCommand::LDO_ADD] = new LdoAddCommand();
  scheduler_command_table_[SchedulerCommand::LDO_REMOVE] = new LdoRemoveCommand();
  scheduler_command_table_[SchedulerCommand::PARTITION_ADD] = new PartitionAddCommand();
  scheduler_command_table_[SchedulerCommand::PARTITION_REMOVE] = new PartitionRemoveCommand();
  scheduler_command_table_[SchedulerCommand::TERMINATE] = new TerminateCommand();
  scheduler_command_table_[SchedulerCommand::DEFINED_TEMPLATE] = new DefinedTemplateCommand();
  scheduler_command_table_[SchedulerCommand::SAVE_DATA] = new SaveDataCommand();
  scheduler_command_table_[SchedulerCommand::LOAD_DATA] = new LoadDataCommand();
  scheduler_command_table_[SchedulerCommand::PREPARE_REWIND] = new PrepareRewindCommand();
  scheduler_command_table_[SchedulerCommand::START_COMMAND_TEMPLATE] = new StartCommandTemplateCommand(); // NOLINT
  scheduler_command_table_[SchedulerCommand::END_COMMAND_TEMPLATE] = new EndCommandTemplateCommand(); //NOLINT
  scheduler_command_table_[SchedulerCommand::SPAWN_COMMAND_TEMPLATE] = new SpawnCommandTemplateCommand(); // NOLINT
  scheduler_command_table_[SchedulerCommand::REQUEST_STAT] = new RequestStatCommand();
}

worker_id_t Worker::id() {
  return id_;
}

void Worker::set_id(worker_id_t id) {
  id_ = id;
}

void Worker::set_ip_address(std::string ip) {
  ip_address_ = ip;
}

PhysicalDataMap* Worker::data_map() {
  return &data_map_;
}

void Worker::AddJobToGraph(Job* job) {
  boost::unique_lock<boost::recursive_mutex> lock(job_graph_mutex_);

  // Print periteration timer stats - omidm
  static uint64_t calculate_dt_counter = 0;
  if (job->name() == "Compute:calculate_dt") {
    if ((++calculate_dt_counter) % 8 == 1) {
      PrintTimerStat();
    }
  }

  // TODO(quhang): when a job is received.
  StatAddJob();
  assert(job != NULL);
  job_id_t job_id = job->id().elem();
  dbg(DBG_WORKER_FD,
      DBG_WORKER_FD_S"Job(%s, #%d) is added to the local job graph.\n",
      job->name().c_str(), job_id);
  assert(job_id != DUMB_JOB_ID);

  // Add vertex for the new job.
  WorkerJobVertex* vertex = NULL;
  if (worker_job_graph_.HasVertex(job_id)) {
    // The job is in the graph but not received.
    worker_job_graph_.GetVertex(job_id, &vertex);
    assert(vertex->entry()->get_job() == NULL);
    switch (vertex->entry()->get_state()) {
      case WorkerJobEntry::PENDING: {
        // If the job is a receive job, add a dumb edge.
        if (dynamic_cast<RemoteCopyReceiveJob*>(job)) {  // NOLINT
          worker_job_graph_.AddEdge(DUMB_JOB_ID, job_id);
        }
        break;
      }
      case WorkerJobEntry::PENDING_DATA_RECEIVED: {
        // Flag shows that the data is already received.
        SerializedData* ser_data = NULL;
        data_version_t version;
        if (data_exchanger_->ReceiveSerializedData(
                job_id, &ser_data, version)) {
          RemoteCopyReceiveJob* receive_job
              = dynamic_cast<RemoteCopyReceiveJob*>(job);  // NOLINT
          assert(receive_job != NULL);
          receive_job->set_serialized_data(ser_data);
          receive_job->set_data_version(version);
        } else {
          assert(false);
        }
        break;
      }
      default:
        assert(false);
    }
  } else {
    // The job is new.
    worker_job_graph_.AddVertex(job_id, new WorkerJobEntry());
    worker_job_graph_.GetVertex(job_id, &vertex);
    if (dynamic_cast<RemoteCopyReceiveJob*>(job)) {  // NOLINT
      worker_job_graph_.AddEdge(DUMB_JOB_ID, job_id);
    }
  }
  vertex->entry()->set_job_id(job_id);
  vertex->entry()->set_job(job);
  vertex->entry()->set_state(WorkerJobEntry::BLOCKED);
  // Add edges for the new job.
  IDSet<job_id_t> before_set = job->before_set();
  for (IDSet<job_id_t>::IDSetIter iter = before_set.begin();
       iter != before_set.end();
       ++iter) {
    job_id_t before_job_id = *iter;
    if (InFinishHintSet(before_job_id)) {
      continue;
    }
    WorkerJobVertex* before_job_vertex = NULL;
    if (worker_job_graph_.HasVertex(before_job_id)) {
      // The job is already known.
      worker_job_graph_.GetVertex(before_job_id, &before_job_vertex);
    } else {
      if (IDMaker::SchedulerProducedJobID(before_job_id)) {
        // Local job is acknowledged locally.
        continue;
      }
      // The job is unknown.
      worker_job_graph_.AddVertex(before_job_id, new WorkerJobEntry());
      worker_job_graph_.GetVertex(before_job_id, &before_job_vertex);
      before_job_vertex->entry()->set_job_id(before_job_id);
      before_job_vertex->entry()->set_job(NULL);
      before_job_vertex->entry()->set_state(WorkerJobEntry::PENDING);
    }
    // If that job is not finished.
    if (before_job_vertex->entry()->get_state() != WorkerJobEntry::FINISH) {
      worker_job_graph_.AddEdge(before_job_vertex, vertex);
    }
  }
  // If the job has no dependency, it is ready.
  if (vertex->incoming_edges()->empty()) {
    vertex->entry()->set_state(WorkerJobEntry::READY);
    ResolveDataArray(job);
    int success_flag = worker_manager_->PushJob(job);
    vertex->entry()->set_job(NULL);
    assert(success_flag);
  }
}

void Worker::ClearAfterSet(WorkerJobVertex* vertex) {
  boost::unique_lock<boost::recursive_mutex> lock(job_graph_mutex_);

  WorkerJobEdge::Map* outgoing_edges = vertex->outgoing_edges();
  // Deletion inside loop is dangerous.
  std::list<WorkerJobVertex*> deletion_list;
  for (WorkerJobEdge::Iter iter = outgoing_edges->begin();
       iter != outgoing_edges->end();
       ++iter) {
    WorkerJobVertex* after_job_vertex = (iter->second)->end_vertex();
    assert(after_job_vertex != NULL);
    deletion_list.push_back(after_job_vertex);
  }
  std::list<Job*> job_list;
  for (std::list<WorkerJobVertex*>::iterator iter = deletion_list.begin();
       iter != deletion_list.end();
       ++iter) {
    WorkerJobVertex* after_job_vertex = *iter;
    worker_job_graph_.RemoveEdge(vertex, after_job_vertex);
    if (after_job_vertex->incoming_edges()->empty()) {
      after_job_vertex->entry()->set_state(WorkerJobEntry::READY);
      assert(after_job_vertex->entry()->get_job() != NULL);
      ResolveDataArray(after_job_vertex->entry()->get_job());
      job_list.push_back(after_job_vertex->entry()->get_job());
      after_job_vertex->entry()->set_job(NULL);
    }
  }
  bool success_flag = worker_manager_->PushJobList(&job_list);
  assert(success_flag);
}

void Worker::NotifyLocalJobDone(Job* job) {
  {
    boost::unique_lock<boost::recursive_mutex> lock(job_graph_mutex_);
    StatEndJob(1);

    job_id_t job_id = job->id().elem();
    data_map_.ReleaseAccess(job_id);
    // Job done for unknown job is not handled.
    if (!worker_job_graph_.HasVertex(job_id)) {
      // The job must be in the local job graph.
      assert(false);
      return;
    }
    WorkerJobVertex* vertex = NULL;
    worker_job_graph_.GetVertex(job_id, &vertex);
    assert(vertex->incoming_edges()->empty());
    ClearAfterSet(vertex);
    delete vertex->entry();
    worker_job_graph_.RemoveVertex(job_id);
    if (!IDMaker::SchedulerProducedJobID(job_id)) {
      AddFinishHintSet(job_id);
    }
    // vertex->entry()->set_state(WorkerJobEntry::FINISH);
    // vertex->entry()->set_job(NULL);
  }

  Parameter params;
  SaveDataJob *j = dynamic_cast<SaveDataJob*>(job); // NOLINT
  if (j != NULL) {
    SaveDataJobDoneCommand cm(j->id(), j->run_time(), j->wait_time(), j->max_alloc(),
                              ID<checkpoint_id_t>(j->checkpoint_id()), j->handle());
    client_->SendCommand(&cm);
  } else if ((!IDMaker::SchedulerProducedJobID(job->id().elem())) || (!job->sterile())) {
    JobDoneCommand cm(job->id(), job->run_time(), job->wait_time(), job->max_alloc(), false);
    client_->SendCommand(&cm);
  }

  delete job;
}

void Worker::NotifyJobDone(job_id_t job_id, bool final) {
  dbg(DBG_WORKER_FD,
      DBG_WORKER_FD_S"Job(#%d) is removed in the local job graph.\n", job_id);
  if (IDMaker::SchedulerProducedJobID(job_id)) {
    // Jobdone command for local job is not handled.
    return;
  }

  boost::unique_lock<boost::recursive_mutex> lock(job_graph_mutex_);

  if (final) {
    // Job done for unknown job is not handled.
    if (!worker_job_graph_.HasVertex(job_id)) {
      return;
    }
    WorkerJobVertex* vertex = NULL;
    worker_job_graph_.GetVertex(job_id, &vertex);
    assert(vertex->incoming_edges()->empty());
    assert(vertex->entry()->get_job() == NULL);
    if (vertex->entry()->get_state() != WorkerJobEntry::FINISH) {
      ClearAfterSet(vertex);
    }
    delete vertex->entry();
    worker_job_graph_.RemoveVertex(job_id);
  } else {
    if (worker_job_graph_.HasVertex(job_id)) {
      WorkerJobVertex* vertex = NULL;
      worker_job_graph_.GetVertex(job_id, &vertex);
      assert(vertex->incoming_edges()->empty());
      assert(vertex->entry()->get_job() == NULL);
      vertex->entry()->set_state(WorkerJobEntry::FINISH);
      ClearAfterSet(vertex);
    } else {
      AddFinishHintSet(job_id);
    }
  }
}

void Worker::NotifyTransmissionDone(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(job_graph_mutex_);

  SerializedData* ser_data = NULL;
  data_version_t version;
  WorkerJobVertex* vertex = NULL;
  // Can a receive-job receives multiple data?
  if (worker_job_graph_.HasVertex(job_id)) {
    worker_job_graph_.GetVertex(job_id, &vertex);
    switch (vertex->entry()->get_state()) {
      case WorkerJobEntry::PENDING: {
        // The job is already in the graph and not received.
        vertex->entry()->set_state(WorkerJobEntry::PENDING_DATA_RECEIVED);
        break;
      }
      case WorkerJobEntry::BLOCKED: {
        // The job is already in the graph and received.
        assert(vertex->entry()->get_job() != NULL);
        RemoteCopyReceiveJob* receive_job =
            dynamic_cast<RemoteCopyReceiveJob*>(vertex->entry()->get_job());  // NOLINT
        assert(receive_job != NULL);
        if (data_exchanger_->ReceiveSerializedData(
                job_id, &ser_data, version)) {
          receive_job->set_serialized_data(ser_data);
          receive_job->set_data_version(version);
        } else {
          assert(false);
        }
        // Remove edge and should be ready now.
        worker_job_graph_.RemoveEdge(DUMB_JOB_ID, job_id);
        if (vertex->incoming_edges()->empty()) {
          vertex->entry()->set_state(WorkerJobEntry::READY);
          ResolveDataArray(receive_job);
          int success_flag = worker_manager_->PushJob(receive_job);
          vertex->entry()->set_job(NULL);
          assert(success_flag);
        }
        break;
      }
      default:
        assert(false);
    }  // End switch.
  } else {
    // The job is not in the graph and not received.
    worker_job_graph_.AddVertex(job_id, new WorkerJobEntry());
    worker_job_graph_.GetVertex(job_id, &vertex);
    vertex->entry()->set_job_id(job_id);
    vertex->entry()->set_job(NULL);
    vertex->entry()->set_state(WorkerJobEntry::PENDING_DATA_RECEIVED);
  }
}


void Worker::AddFinishHintSet(const job_id_t job_id) {
  if (hint_set_.find(job_id) != hint_set_.end()) {
    return;
  }
  if (hint_set_.size() < max_hint_size_) {
    hint_set_.insert(job_id);
    hint_queue_.push_back(job_id);
  } else {
    hint_set_.erase(hint_queue_.front());
    hint_queue_.pop_front();
    hint_set_.insert(job_id);
    hint_queue_.push_back(job_id);
  }
}

bool Worker::InFinishHintSet(const job_id_t job_id) {
  return hint_set_.find(job_id) != hint_set_.end();
}


void Worker::ClearBlockedJobs(WorkerJobGraph* job_graph) {
  std::list<job_id_t> list_to_remove;
  typename WorkerJobVertex::Iter iter = job_graph->begin();
  for (; iter != job_graph->end(); ++iter) {
    WorkerJobEntry* job_entry = iter->second->entry();
    if (job_entry->get_state() != WorkerJobEntry::CONTROL &&
        job_entry->get_state() != WorkerJobEntry::READY) {
      list_to_remove.push_back(job_entry->get_job_id());
      if (job_entry->get_job()) {
        delete job_entry->get_job();
      }
      job_entry->set_job(NULL);
      delete job_entry;
    }
  }

  std::list<job_id_t>::iterator it = list_to_remove.begin();
  for (; it != list_to_remove.end(); ++it) {
    job_graph->RemoveVertex(*it);
  }
}

bool Worker::IsEmptyGraph(WorkerJobGraph* job_graph) {
  typename WorkerJobVertex::Iter iter = job_graph->begin();
  for (; iter != job_graph->end(); ++iter) {
    WorkerJobEntry* job_entry = iter->second->entry();
    if (job_entry->get_state() != WorkerJobEntry::CONTROL) {
      return false;
    }
  }
  return true;
}

void Worker::StatAddJob() {
  boost::unique_lock<boost::recursive_mutex> lock(stat_mutex_);
  // printf("add %d %d %d %d %d\n",
  //        stat_busy_cores_, stat_blocked_cores_, stat_idle_cores_,
  //        stat_blocked_job_num_, stat_ready_job_num_);
  ++stat_blocked_job_num_;
  if (stat_idle_cores_ != 0) {
    --stat_idle_cores_;
    ++stat_blocked_cores_;
    // timer::StartTimer(timer::kSumCyclesBlock);
    block_timer_.Start(1);
  }
  // printf("#add %d %d %d %d %d\n",
  //        stat_busy_cores_, stat_blocked_cores_, stat_idle_cores_,
  //        stat_blocked_job_num_, stat_ready_job_num_);
}
void Worker::StatDispatchJob(int len) {
  boost::unique_lock<boost::recursive_mutex> lock(stat_mutex_);
  // printf("%d dis %d %d %d %d %d\n", len,
  //        stat_busy_cores_, stat_blocked_cores_, stat_idle_cores_,
  //        stat_blocked_job_num_, stat_ready_job_num_);
  stat_blocked_job_num_ -= len;
  stat_ready_job_num_ += len;
  if (stat_blocked_cores_ > 0) {
    int release_cores = std::min(stat_blocked_cores_, len);
    if (release_cores > 0) {
      stat_blocked_cores_ -= release_cores;
      // timer::StopTimer(timer::kSumCyclesBlock, release_cores);
      block_timer_.Stop(release_cores);
      stat_busy_cores_ += release_cores;
      // timer::StartTimer(timer::kSumCyclesRun, release_cores);
      run_timer_.Start(release_cores);
    }
  }
  // printf("#dis %d %d %d %d %d\n",
  //        stat_busy_cores_, stat_blocked_cores_, stat_idle_cores_,
  //        stat_blocked_job_num_, stat_ready_job_num_);
}
void Worker::StatEndJob(int len) {
  boost::unique_lock<boost::recursive_mutex> lock(stat_mutex_);
  // printf("%d end %d %d %d %d %d\n", len,
  //        stat_busy_cores_, stat_blocked_cores_, stat_idle_cores_,
  //        stat_blocked_job_num_, stat_ready_job_num_);
  using std::min;
  stat_ready_job_num_ -= len;
  int busy_cores = min(stat_ready_job_num_,
                       WorkerManager::across_job_parallism);
  int blocked_cores = min(stat_blocked_job_num_,
                          WorkerManager::across_job_parallism - busy_cores);
  int idle_cores =
      WorkerManager::across_job_parallism - busy_cores - blocked_cores;
  if (busy_cores != stat_busy_cores_) {
    // timer::StopTimer(timer::kSumCyclesRun, stat_busy_cores_ - busy_cores);
    run_timer_.Stop(stat_busy_cores_ - busy_cores);
  }
  if (blocked_cores != stat_blocked_cores_) {
    // timer::StartTimer(timer::kSumCyclesBlock, blocked_cores - stat_blocked_cores_);
    block_timer_.Start(blocked_cores - stat_blocked_cores_);
  }
  stat_busy_cores_ = busy_cores;
  stat_blocked_cores_ = blocked_cores;
  stat_idle_cores_ = idle_cores;
  // printf("#end %d %d %d %d %d\n",
  //        stat_busy_cores_, stat_blocked_cores_, stat_idle_cores_,
  //        stat_blocked_job_num_, stat_ready_job_num_);
}

// The unit is in nano-second.
void Worker::GetTimerStat(int64_t* idle, int64_t* block, int64_t* run) {
  boost::unique_lock<boost::recursive_mutex> lock(stat_mutex_);
  static int64_t l_idle = 0, l_block = 0, l_run = 0;
  // int64_t c_block = timer::ReadTimer(timer::kSumCyclesBlock);
  int64_t c_block = block_timer_.Read();
  // int64_t c_run = timer::ReadTimer(timer::kSumCyclesRun);
  int64_t c_run = run_timer_.Read();
  // int64_t c_idle = timer::ReadTimer(timer::kSumCyclesTotal) - c_block - c_run;
  int64_t c_idle = total_timer_.Read() - c_block - c_run;
  *idle = c_idle - l_idle;
  *block = c_block - l_block;
  *run = c_run - l_run;
  l_idle = c_idle;
  l_block = c_block;
  l_run = c_run;
}

// The unit is in nano-second.
void Worker::PrintTimerStat() {
  boost::unique_lock<boost::recursive_mutex> lock(stat_mutex_);
  std::string file_name = int2string(id_) + "_main_timers.txt";
  static FILE* temp = fopen(file_name.c_str(), "w");
  static int64_t l_idle = 0, l_block = 0, l_run = 0;
  int64_t c_block = block_timer_.Read();
  int64_t c_run = run_timer_.Read();
  int64_t c_idle = total_timer_.Read() - c_block - c_run;
  int64_t idle = c_idle - l_idle;
  int64_t block = c_block - l_block;
  int64_t run = c_run - l_run;
  l_idle = c_idle;
  l_block = c_block;
  l_run = c_run;
  fprintf(temp, "run_time: %.9f block_time: %.9f idle_time: %.9f \n",
      static_cast<double>(run) / 1e9,
      static_cast<double>(block) / 1e9,
      static_cast<double>(idle) / 1e9);
  fflush(temp);
}


}  // namespace nimbus
