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

#include <boost/functional/hash.hpp>
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
    log_.InitTime();
    id_ = -1;
    ip_address_ = NIMBUS_RECEIVER_KNOWN_IP;
    cache_log = NULL;
    worker_manager_ = new WorkerManager();
    DUMB_JOB_ID = std::numeric_limits<job_id_t>::max();
    worker_job_graph_.AddVertex(
        DUMB_JOB_ID,
        new WorkerJobEntry(
            DUMB_JOB_ID, NULL, WorkerJobEntry::CONTROL));
}

Worker::~Worker() {
  worker_job_graph_.RemoveVertex(DUMB_JOB_ID);
  delete worker_manager_;
}

void Worker::Run() {
#ifdef CACHE_LOG
  std::stringstream msg;
  msg << "~~~ Worker starts : " << cache_log->GetTime();
  cache_log->WriteToFile(msg.str());
#endif

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

  JobList local_job_done_list;
  while (true) {
    // Process command.
    SchedulerCommand* comm = client_->receiveCommand();
    int quota = SCHEDULER_COMMAND_GROUP_QUOTA;
    while (comm != NULL) {
      dbg(DBG_WORKER, "Receives command from scheduler: %s\n",
          comm->toStringWTags().c_str());
      dbg(DBG_WORKER_FD,
          DBG_WORKER_FD_S"Scheduler command arrives(%s).\n",
          comm->name().c_str());
      ProcessSchedulerCommand(comm);
      delete comm;
      if (--quota == 0) {
        break;
      }
      comm = client_->receiveCommand();
    }
    // Poll jobs that finish receiving.
    job_id_t receive_job_id;
    while (data_exchanger_->GetReceiveEvent(&receive_job_id)) {
      dbg(DBG_WORKER_FD,
          DBG_WORKER_FD_S"Receive-job transmission is done(job #%d)\n",
          receive_job_id);
      NotifyTransmissionDone(receive_job_id);
    }
    // Job done processing.
    worker_manager_->GetLocalJobDoneList(&local_job_done_list);
    for (JobList::iterator index = local_job_done_list.begin();
         index != local_job_done_list.end();
         ++index) {
      dbg(DBG_WORKER_FD,
          DBG_WORKER_FD_S"Local job-done notification arrives(job #%d).\n",
          *index);
      NotifyLocalJobDone(*index);
    }
    local_job_done_list.clear();
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
  // DumpVersionInformation(job, da, &version_log_, "version_in");
  // DumpDataHashInformation(job, da, &data_hash_log_, "hash_in");

  IDSet<physical_data_id_t> write = job->write_set();
  for (iter = write.begin(); iter != write.end(); iter++) {
    job->data_array.push_back(
        data_map_.AcquireAccess(*iter, job->id().elem(),
                                PhysicalDataMap::WRITE));
  }
  DumpDataOrderInformation(job, job->data_array, &data_hash_log_, "data_order");
}

void Worker::ProcessSchedulerCommand(SchedulerCommand* cm) {
  switch (cm->type()) {
    case SchedulerCommand::HANDSHAKE:
      ProcessHandshakeCommand(reinterpret_cast<HandshakeCommand*>(cm));
      break;
    case SchedulerCommand::JOB_DONE:
      ProcessJobDoneCommand(reinterpret_cast<JobDoneCommand*>(cm));
      break;
    case SchedulerCommand::COMPUTE_JOB:
      ProcessComputeJobCommand(reinterpret_cast<ComputeJobCommand*>(cm));
      break;
    case SchedulerCommand::CREATE_DATA:
      ProcessCreateDataCommand(reinterpret_cast<CreateDataCommand*>(cm));
      break;
    case SchedulerCommand::REMOTE_COPY_SEND:
      ProcessRemoteCopySendCommand(reinterpret_cast<RemoteCopySendCommand*>(cm));
      break;
    case SchedulerCommand::REMOTE_COPY_RECEIVE:
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
    default:
      std::cout << "ERROR: " << cm->toString() <<
        " have not been implemented in ProcessSchedulerCommand yet." <<
        std::endl;
  }
}

// Processes handshake command. Configures the worker based on the handshake
// command and responds by sending another handshake command.
void Worker::ProcessHandshakeCommand(HandshakeCommand* cm) {
  ID<port_t> port(listening_port_);
  HandshakeCommand new_cm = HandshakeCommand(cm->worker_id(),
      // boost::asio::ip::host_name(), port);
      // "127.0.0.1", port);
      ip_address_, port);
  client_->sendCommand(&new_cm);

  id_ = cm->worker_id().elem();
  id_maker_.Initialize(id_);

  std::ostringstream ss;
  ss << id_;
  // TODO(quhang) thread-safety(log).
  version_log_.set_file_name(ss.str() + "_version_log.txt");
  data_hash_log_.set_file_name(ss.str() + "_data_hash_log.txt");
}

// Processes jobdone command. Moves a job from blocked queue to ready queue if
// its before set is satisfied.
void Worker::ProcessJobDoneCommand(JobDoneCommand* cm) {
  NotifyJobDone(cm->job_id().elem());
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
  job->set_parameters(cm->params());
  job->set_sterile(cm->sterile());
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

  Job * job = new CreateDataJob();
  job->set_name("CreateData:" + cm->data_name());
  job->set_id(cm->job_id());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
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
  job->set_receive_job_id(cm->receive_job_id());
  job->set_to_worker_id(cm->to_worker_id());
  job->set_to_ip(cm->to_ip());
  job->set_to_port(cm->to_port());
  IDSet<physical_data_id_t> read_set;
  read_set.insert(cm->from_physical_data_id().elem());
  job->set_read_set(read_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
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
  job->set_after_set(cm->after_set());
#ifndef MUTE_LOG
  timer_.Start(job->id().elem());
#endif  // MUTE_LOG
  AddJobToGraph(job);
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
  job->set_after_set(cm->after_set());
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
    GeometricRegion r = *(cm->region());
    if (!ldo_map_->RemovePartition(cm->id().elem()))
        dbg(DBG_ERROR, "Worker could not remove partition %i from ldo map\n", cm->id().elem());
}

void Worker::ProcessTerminateCommand(TerminateCommand* cm) {
#ifdef CACHE_LOG
  std::stringstream msg;
  msg << "~~~ Completed application : " << cache_log->GetTime();
  cache_log->WriteToFile(msg.str());
#endif
  // profiler_thread_->interrupt();
  // profiler_thread_->join();
  // ProfilerMalloc::Exit();
  exit(cm->exit_status().elem());
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
  client_->run();
  // client_thread_ = new boost::thread(
  //     boost::bind(&SchedulerClient::run, client_));

  // profiler_thread_ = new boost::thread(
  //     boost::bind(&Profiler::Run, &profiler_, client_, id_));
}

void Worker::LoadSchedulerCommands() {
  // std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  scheduler_command_table_.push_back(new HandshakeCommand());
  scheduler_command_table_.push_back(new JobDoneCommand());
  scheduler_command_table_.push_back(new ComputeJobCommand());
  scheduler_command_table_.push_back(new CreateDataCommand);
  scheduler_command_table_.push_back(new RemoteCopySendCommand());
  scheduler_command_table_.push_back(new RemoteCopyReceiveCommand());
  scheduler_command_table_.push_back(new LocalCopyCommand());
  scheduler_command_table_.push_back(new LdoAddCommand());
  scheduler_command_table_.push_back(new LdoRemoveCommand());
  scheduler_command_table_.push_back(new PartitionAddCommand());
  scheduler_command_table_.push_back(new PartitionRemoveCommand());
  scheduler_command_table_.push_back(new TerminateCommand());
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
    WorkerJobVertex* before_job_vertex = NULL;
    if (worker_job_graph_.HasVertex(before_job_id)) {
      // The job is already known.
      worker_job_graph_.GetVertex(before_job_id, &before_job_vertex);
    } else {
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
  for (std::list<WorkerJobVertex*>::iterator iter = deletion_list.begin();
       iter != deletion_list.end();
       ++iter) {
    WorkerJobVertex* after_job_vertex = *iter;
    worker_job_graph_.RemoveEdge(vertex, after_job_vertex);
    if (after_job_vertex->incoming_edges()->empty()) {
      after_job_vertex->entry()->set_state(WorkerJobEntry::READY);
      assert(after_job_vertex->entry()->get_job() != NULL);
      ResolveDataArray(after_job_vertex->entry()->get_job());
      int success_flag =
          worker_manager_->PushJob(after_job_vertex->entry()->get_job());
#ifndef MUTE_LOG
      double wait_time = timer_.Stop(after_job_vertex->entry()->get_job()->id().elem());
      after_job_vertex->entry()->get_job()->set_wait_time(wait_time);
#endif  // MUTE_LOG
      after_job_vertex->entry()->set_job(NULL);
      assert(success_flag);
    }
  }
}

void Worker::NotifyLocalJobDone(Job* job) {
  Parameter params;
  JobDoneCommand cm(job->id(), job->after_set(), params, job->run_time(), job->wait_time(),
                    job->max_alloc());
  client_->sendCommand(&cm);
  job_id_t job_id = job->id().elem();
  // Job done for unknown job is not handled.
  if (!worker_job_graph_.HasVertex(job_id)) {
    return;
  }
  WorkerJobVertex* vertex = NULL;
  worker_job_graph_.GetVertex(job_id, &vertex);
  assert(vertex->incoming_edges()->empty());
  ClearAfterSet(vertex);
  vertex->entry()->set_state(WorkerJobEntry::FINISH);
  delete job;
}

void Worker::NotifyJobDone(job_id_t job_id) {
  dbg(DBG_WORKER_FD,
      DBG_WORKER_FD_S"Job(#%d) is removed in the local job graph.\n", job_id);
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
  worker_job_graph_.RemoveVertex(job_id);
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

}  // namespace nimbus
