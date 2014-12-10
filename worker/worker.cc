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
#include <cstdio>
#include <sstream>
#include <string>
#include <ctime>
#include <list>
#include <limits>
#include "shared/profiler_malloc.h"
#include "worker/worker.h"
#include "worker/worker_ldo_map.h"
#include "worker/worker_manager.h"
#include "worker/util_dumping.h"
#include "data/physbam/physbam_data.h"

#define SCHEDULER_COMMAND_GROUP_QUOTA 10

using boost::hash;

namespace nimbus {

/*
// Comment(quhang): I moved these three functions to a seperate file:
// worker/util_dumping.cc.
// So that these utilities can be called outside "Worker" class.
void DumpVersionInformation(Job *job, const DataArray& da, Log *log,
                            std::string tag);

void DumpDataHashInformation(Job *job, const DataArray& da, Log *log,
                             std::string tag);

void DumpDataOrderInformation(Job *job, const DataArray& da, Log *log,
                              std::string tag);
 */

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
    log_.InitTime();
    id_ = -1;
    ip_address_ = NIMBUS_RECEIVER_KNOWN_IP;
    cache_log = NULL;
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

void Worker::PrintTimeStamp(const char* format, ...) {
  va_list args;
  va_start(args, format);
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(event_log, "%f ", time_sum);
  vfprintf(event_log, format, args);
  va_end(args);
}

void Worker::Run() {
  if (cache_log) {
    std::stringstream msg;
    msg << "~~~ Worker starts : " << cache_log->GetTime();
    cache_log->WriteToFile(msg.str());
  }

  std::cout << "Running the Worker" << std::endl;

  SetupSchedulerInterface();

  ldo_map_ = new WorkerLdoMap();
  application_->Start(client_, &id_maker_, ldo_map_);

  SetupDataExchangerInterface();

  WorkerCoreProcessor();
}

void Worker::WorkerCoreProcessor() {
  std::cout << "Base Worker Core Processor" << std::endl;
  worker_manager_->worker_ = this;
  worker_manager_->SetLoggingInterface(&log_, &version_log_, &data_hash_log_,
                                       cache_log,
                                       &timer_);
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
      SchedulerCommand *comm = *iter;
      dbg(DBG_WORKER, "Receives command from scheduler: %s\n",
          comm->ToString().c_str());
      dbg(DBG_WORKER_FD,
          DBG_WORKER_FD_S"Scheduler command arrives(%s).\n",
          comm->name().c_str());
      process_jobs = true;
      ProcessSchedulerCommand(comm);
      delete comm;
    }
    // Poll jobs that finish receiving.
    job_id_t receive_job_id;
    int quota = 10;
    while (data_exchanger_->GetReceiveEvent(&receive_job_id)) {
      dbg(DBG_WORKER_FD,
          DBG_WORKER_FD_S"Receive-job transmission is done(job #%d)\n",
          receive_job_id);
      PrintTimeStamp("io_done %lu\n", receive_job_id);
      process_jobs = true;
      NotifyTransmissionDone(receive_job_id);
      if (--quota <= 0) {
        break;
      }
    }
    // Job done processing.
    worker_manager_->GetLocalJobDoneList(&local_job_done_list);
    quota = 10;
    while (!local_job_done_list.empty()) {
      Job* job = local_job_done_list.front();
      local_job_done_list.pop_front();
      PrintTimeStamp("local_done %lu\n", job->id().elem());
      process_jobs = true;
      NotifyLocalJobDone(job);
      if (--quota <= 0) {
        break;
      }
    }
    if (!process_jobs) {
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

  IDSet<physical_data_id_t> read = job->read_set();
  for (iter = read.begin(); iter != read.end(); iter++) {
    job->data_array.push_back(
        data_map_.AcquireAccess(*iter, job->id().elem(),
                                PhysicalDataMap::READ));
  }
  IDSet<physical_data_id_t> write = job->write_set();
  for (iter = write.begin(); iter != write.end(); iter++) {
    job->data_array.push_back(
        data_map_.AcquireAccess(*iter, job->id().elem(),
                                PhysicalDataMap::WRITE));
  }
  // DumpVersionInformation(job, da, &version_log_, "version_in");
  // DumpDataHashInformation(job, da, &data_hash_log_, "hash_in");
  // DumpDataOrderInformation(job, job->data_array, &data_hash_log_, "data_order");
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
  // TODO(quhang) thread-safety(log).
  event_log = fopen((wstr + "_event_fe.txt").c_str(), "w");
  alloc_log = fopen((wstr + "_data_objects.txt").c_str(), "w");
  worker_manager_->SetEventLog(wstr);
  version_log_.set_file_name(wstr + "_version_log.txt");
  data_hash_log_.set_file_name(wstr + "_data_hash_log.txt");
  data_exchanger_->WriteTimeDriftToLog(time - cm->time());
}

// Processes jobdone command. Moves a job from blocked queue to ready queue if
// its before set is satisfied.
void Worker::ProcessJobDoneCommand(JobDoneCommand* cm) {
  PrintTimeStamp("recv_job %s %lu\n", "JOB_DONE", cm->job_id().elem());
  NotifyJobDone(cm->job_id().elem(), cm->final());
}

// Processes computejob command. Generates the corresponding job and pushes it
// to the blocking queue.
void Worker::ProcessComputeJobCommand(ComputeJobCommand* cm) {
  Job* job = application_->CloneJob(cm->job_name());
  job->set_name("Compute:" + cm->job_name());
  job->set_id(cm->job_id());
  // TODO(print_log): Receive a compute job.
  PrintTimeStamp("recv_job %s %lu\n", job->name().c_str(), job->id().elem());
  job->set_read_set(cm->read_set());
  job->set_write_set(cm->write_set());
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  job->set_future_job_id(cm->future_job_id());
  job->set_sterile(cm->sterile());
  job->set_region(cm->region());
  job->set_parameters(cm->params());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
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

  data_map_.AddMapping(data->physical_id(), data);
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(alloc_log, "%f : %lu\n", time_sum, sizeof(*data));

  Job * job = new CreateDataJob();
  job->set_name("CreateData:" + cm->data_name());
  job->set_id(cm->job_id());
  PrintTimeStamp("recv_job %s %lu\n", job->name().c_str(), job->id().elem());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
  AddJobToGraph(job);
}

void Worker::ProcessRemoteCopySendCommand(RemoteCopySendCommand* cm) {
  RemoteCopySendJob * job = new RemoteCopySendJob(data_exchanger_, application_);
  data_exchanger_->AddContactInfo(cm->to_worker_id().elem(),
      cm->to_ip(), cm->to_port().elem());
  job->set_name("RemoteCopySend");
  job->set_id(cm->job_id());
  PrintTimeStamp("recv_job %s %lu\n", job->name().c_str(), job->id().elem());
  job->set_receive_job_id(cm->receive_job_id());
  job->set_to_worker_id(cm->to_worker_id());
  job->set_to_ip(cm->to_ip());
  job->set_to_port(cm->to_port());
  IDSet<physical_data_id_t> read_set;
  read_set.insert(cm->from_physical_data_id().elem());
  job->set_read_set(read_set);
  job->set_before_set(cm->before_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
  AddJobToGraph(job);
}

void Worker::ProcessRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* cm) {
  Job * job = new RemoteCopyReceiveJob(application_);
  job->set_name("RemoteCopyReceive");
  job->set_id(cm->job_id());
  PrintTimeStamp("recv_job %s %lu\n", job->name().c_str(), job->id().elem());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->to_physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
  AddJobToGraph(job);
}

void Worker::ProcessSaveDataCommand(SaveDataCommand* cm) {
  SaveDataJob * job = new SaveDataJob(ddb_, application_);
  job->set_name("SaveData");
  job->set_id(cm->job_id());
  PrintTimeStamp("recv_job %s %lu\n", job->name().c_str(), job->id().elem());
  job->set_checkpoint_id(cm->checkpoint_id().elem());
  IDSet<physical_data_id_t> read_set;
  read_set.insert(cm->from_physical_data_id().elem());
  job->set_read_set(read_set);
  job->set_before_set(cm->before_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
  AddJobToGraph(job);
}

void Worker::ProcessLoadDataCommand(LoadDataCommand* cm) {
  LoadDataJob * job = new LoadDataJob(ddb_, application_);
  job->set_name("LoadData");
  job->set_id(cm->job_id());
  PrintTimeStamp("load_job %s %lu\n", job->name().c_str(), job->id().elem());
  job->set_handle(cm->handle());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->to_physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
  AddJobToGraph(job);
}

void Worker::ProcessPrepareRewindCommand(PrepareRewindCommand* cm) {
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
      PrintTimeStamp("local_done %lu\n", job->id().elem());
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

void Worker::ProcessLocalCopyCommand(LocalCopyCommand* cm) {
  Job * job = new LocalCopyJob(application_);
  job->set_name("LocalCopy");
  job->set_id(cm->job_id());
  PrintTimeStamp("recv_job %s %lu\n", job->name().c_str(), job->id().elem());
  IDSet<physical_data_id_t> read_set;
  read_set.insert(cm->from_physical_data_id().elem());
  job->set_read_set(read_set);
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->to_physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
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
  if (cache_log) {
    std::stringstream msg;
    msg << "~~~ Completed application : " << cache_log->GetTime();
    cache_log->WriteToFile(msg.str());
  }
  // profiler_thread_->interrupt();
  // profiler_thread_->join();
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
    // TODO(print_log): dispatch a job, newly received with empty before set.
    // if (!(dynamic_cast<CreateDataJob*>(job) || // NOLINT
    //       dynamic_cast<LocalCopyJob*>(job) || // NOLINT
    //       dynamic_cast<RemoteCopySendJob*>(job) || // NOLINT
    //       dynamic_cast<RemoteCopyReceiveJob*>(job))) { // NOLINT
     PrintTimeStamp("dispatch_job(new) %s %lu\n",
                    job->name().c_str(), job->id().elem());
    // }
    int success_flag = worker_manager_->PushJob(job);
#ifndef MUTE_LOG
    double wait_time = timer_.Stop(job->id().elem());
    job->set_wait_time(wait_time);
#endif  // MUTE_LOG
    vertex->entry()->set_job(NULL);
    assert(success_flag);
  }
}

void Worker::ClearAfterSet(WorkerJobVertex* vertex) {
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
      Job* job = after_job_vertex->entry()->get_job();
      PrintTimeStamp("dispatch_job(job_done) %s %lu %lu\n",
                     job->name().c_str(),
                     job->id().elem(),
                     vertex->entry()->get_job_id());
      job_list.push_back(after_job_vertex->entry()->get_job());
#ifndef MUTE_LOG
      double wait_time = timer_.Stop(after_job_vertex->entry()->get_job()->id().elem());
      after_job_vertex->entry()->get_job()->set_wait_time(wait_time);
#endif  // MUTE_LOG
      after_job_vertex->entry()->set_job(NULL);
    }
  }
  bool success_flag = worker_manager_->PushJobList(&job_list);
  assert(success_flag);
}

void Worker::NotifyLocalJobDone(Job* job) {
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
  delete job;
}

void Worker::NotifyJobDone(job_id_t job_id, bool final) {
  dbg(DBG_WORKER_FD,
      DBG_WORKER_FD_S"Job(#%d) is removed in the local job graph.\n", job_id);
  if (IDMaker::SchedulerProducedJobID(job_id)) {
    // Jobdone command for local job is not handled.
    return;
  }
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
          PrintTimeStamp("dispatch_job(new) %s %lu\n",
                         receive_job->name().c_str(), receive_job->id().elem());
          int success_flag = worker_manager_->PushJob(receive_job);
#ifndef MUTE_LOG
          double wait_time = timer_.Stop(receive_job->id().elem());
          receive_job->set_wait_time(wait_time);
#endif  // MUTE_LOG
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


}  // namespace nimbus
