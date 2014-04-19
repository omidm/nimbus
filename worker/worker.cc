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
#include <ctime>
#include "worker/worker.h"
#include "worker/worker_ldo_map.h"
#include "worker/worker_manager.h"
#include "worker/util_dumping.h"
#include "data/physbam/physbam_data.h"

#define MAX_PARALLEL_JOB 10
// Configures the number of threads used to run jobs.
#define CORE_NUMBER 1

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
  std::cout << "Base Worker Core Processor" << std::endl;
  WorkerManager worker_manager;
  worker_manager.worker_ = this;
  worker_manager.SetLoggingInterface(&log_, &version_log_, &data_hash_log_,
                                     &timer_);
  worker_manager.StartWorkerThreads(CORE_NUMBER);

  while (true) {
    SchedulerCommand* comm = client_->receiveCommand();
    if (comm != NULL) {
      std::cout << "Received command: " << comm->toStringWTags()
        << std::endl;
      ProcessSchedulerCommand(comm);
      delete comm;
    }

    ScanBlockedJobs();

    ScanPendingTransferJobs();

    GetJobsToRun(&worker_manager, (size_t)(MAX_PARALLEL_JOB));
  }
}

void Worker::ScanBlockedJobs() {
  JobList::iterator iter;
  for (iter = blocked_jobs_.begin(); iter != blocked_jobs_.end();) {
    if ((*iter)->before_set().size() == 0) {
      if (dynamic_cast<RemoteCopyReceiveJob*>(*iter) == NULL) { // NOLINT
        ready_jobs_.push_back(*iter);
      } else {
        pending_transfer_jobs_.push_back(*iter);
      }

      double wait_time =  timer_.Stop((*iter)->id().elem());
      (*iter)->set_wait_time(wait_time);

      blocked_jobs_.erase(iter++);
    } else {
      ++iter;
    }
  }
}

void Worker::ScanPendingTransferJobs() {
  JobList::iterator iter;
  for (iter = pending_transfer_jobs_.begin(); iter != pending_transfer_jobs_.end();) {
    SerializedData* ser_data;
    data_version_t version;
    if (data_exchanger_->ReceiveSerializedData((*iter)->id().elem(), &ser_data, version)) {
      static_cast<RemoteCopyReceiveJob*>(*iter)->set_serialized_data(ser_data);
      static_cast<RemoteCopyReceiveJob*>(*iter)->set_data_version(version);
      ready_jobs_.push_back(*iter);
      pending_transfer_jobs_.erase(iter++);
    } else {
      ++iter;
    }
  }
}


// Extracts at most "max_num" jobs from queue "ready_jobs", resolves their data
// array, and push them to the worker manager.
void Worker::GetJobsToRun(WorkerManager* worker_manager, size_t max_num) {
  size_t ready_num = ready_jobs_.size();
  for (size_t i = 0; (i < max_num) && (i < ready_num); i++) {
    Job* job = ready_jobs_.front();
    ready_jobs_.pop_front();
    ResolveDataArray(job);
    int success_flag = worker_manager->PushComputationJob(job);
    assert(success_flag);
  }
}

// Extracts data objects from the read/write set to data array.
void Worker::ResolveDataArray(Job* job) {
  job->data_array.clear();
  IDSet<physical_data_id_t>::IDSetIter iter;

  IDSet<physical_data_id_t> read = job->read_set();
  for (iter = read.begin(); iter != read.end(); iter++) {
    job->data_array.push_back(data_map_[*iter]);
  }
  // DumpVersionInformation(job, da, &version_log_, "version_in");
  // DumpDataHashInformation(job, da, &data_hash_log_, "hash_in");

  IDSet<physical_data_id_t> write = job->write_set();
  for (iter = write.begin(); iter != write.end(); iter++) {
    job->data_array.push_back(data_map_[*iter]);
  }
  DumpDataOrderInformation(job, job->data_array, &data_hash_log_, "data_order");
}

/*
// Comment(quhang) the function is moved to:
// Line 70 (worker_thread_finish.cc).
// So that the worker core process doesn't need to clean up a job.
// TODO(quhang) operation on shared data structure "data_map" should be
// synchronized.
void Worker::UpdateDataVersion(Job* job) {
  DataArray daw;
  IDSet<physical_data_id_t> write = job->write_set();
  IDSet<physical_data_id_t>::IDSetIter iter;
  for (iter = write.begin(); iter != write.end(); iter++) {
    daw.push_back(data_map_[*iter]);
  }
  DumpDataHashInformation(job, daw, &data_hash_log_, "hash_out");

  if ((dynamic_cast<CreateDataJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<LocalCopyJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<RemoteCopySendJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<RemoteCopyReceiveJob*>(job) == NULL)) { // NOLINT
    for (iter = write.begin(); iter != write.end(); iter++) {
      Data *d = data_map_[*iter];
      data_version_t version = d->version();
      ++version;
      d->set_version(version);
    }
  }

  Parameter params;
  JobDoneCommand cm(job->id(), job->after_set(), params, job->run_time(), job->wait_time());
  client_->sendCommand(&cm);
  // ProcessJobDoneCommand(&cm);
  delete job;
}
*/

// Comment(quhang): This funciton is moved to
// Line 66, worker_thread_computation.cc.
/*
void Worker::ExecuteJob(Job* job) {
  log_.StartTimer();
  timer_.Start(job->id().elem());
  job->Execute(job->parameters(), job->data_array);
  double run_time = timer_.Stop(job->id().elem());
  log_.StopTimer();

  job->set_run_time(run_time);

  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff),
      "Execute Job, name: %35s  id: %6llu  length(s): %2.3lf  time(s): %6.3lf",
           job->name().c_str(), job->id().elem(), log_.timer(), log_.GetTime());
  log_.WriteToOutputStream(std::string(buff), LOG_INFO);

  char time_buff[LOG_MAX_BUFF_SIZE];
  snprintf(time_buff, sizeof(time_buff),
      "Queue Time: %2.9lf, Run Time: %2.9lf",
      job->wait_time(), job->run_time());
  log_.WriteToOutputStream(std::string(time_buff), LOG_INFO);

}
*/

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

void Worker::ProcessHandshakeCommand(HandshakeCommand* cm) {
  ID<port_t> port(listening_port_);
  HandshakeCommand new_cm = HandshakeCommand(cm->worker_id(),
      // boost::asio::ip::host_name(), port);
      "127.0.0.1", port);
  client_->sendCommand(&new_cm);

  id_ = cm->worker_id().elem();
  id_maker_.Initialize(id_);

  std::ostringstream ss;
  ss << id_;
  version_log_.set_file_name(ss.str() + "_version_log.txt");
  data_hash_log_.set_file_name(ss.str() + "_data_hash_log.txt");
}

void Worker::ProcessJobDoneCommand(JobDoneCommand* cm) {
  JobList::iterator iter;
  for (iter = blocked_jobs_.begin(); iter != blocked_jobs_.end(); iter++) {
    IDSet<job_id_t> req = (*iter)->before_set();
    req.remove(cm->job_id().elem());
    (*iter)->set_before_set(req);
  }
}

void Worker::ProcessComputeJobCommand(ComputeJobCommand* cm) {
  Job * job = application_->CloneJob(cm->job_name());
  job->set_name("Compute:" + cm->job_name());
  job->set_id(cm->job_id());
  job->set_read_set(cm->read_set());
  job->set_write_set(cm->write_set());
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  job->set_parameters(cm->params());
  job->set_sterile(cm->sterile());
  timer_.Start(job->id().elem());
  blocked_jobs_.push_back(job);
}

void Worker::ProcessCreateDataCommand(CreateDataCommand* cm) {
  Data * data = application_->CloneData(cm->data_name());
  data->set_logical_id(cm->logical_data_id().elem());
  data->set_physical_id(cm->physical_data_id().elem());
  // data->set_name(cm->data_name());
  const LogicalDataObject* ldo;
  ldo = ldo_map_->FindLogicalObject(cm->logical_data_id().elem());
  data->set_region(*(ldo->region()));
  AddData(data);

  Job * job = new CreateDataJob();
  job->set_name("CreateData:" + cm->data_name());
  job->set_id(cm->job_id());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  timer_.Start(job->id().elem());
  blocked_jobs_.push_back(job);
}

void Worker::ProcessRemoteCopySendCommand(RemoteCopySendCommand* cm) {
  RemoteCopySendJob * job = new RemoteCopySendJob(data_exchanger_);
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
  timer_.Start(job->id().elem());
  blocked_jobs_.push_back(job);
}

void Worker::ProcessRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* cm) {
  Job * job = new RemoteCopyReceiveJob();
  job->set_name("RemoteCopyReceive");
  job->set_id(cm->job_id());
  IDSet<physical_data_id_t> write_set;
  write_set.insert(cm->to_physical_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  timer_.Start(job->id().elem());
  blocked_jobs_.push_back(job);
}

void Worker::ProcessLocalCopyCommand(LocalCopyCommand* cm) {
  Job * job = new LocalCopyJob();
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
  timer_.Start(job->id().elem());
  blocked_jobs_.push_back(job);
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

void Worker::AddData(Data* data) {
  data_map_[data->physical_id()] = data;
}

void Worker::DeleteData(physical_data_id_t physical_data_id) {
  data_map_.erase(physical_data_id);
}

worker_id_t Worker::id() {
  return id_;
}

void Worker::set_id(worker_id_t id) {
  id_ = id;
}

PhysicalDataMap* Worker::data_map() {
  return &data_map_;
}

}  // namespace nimbus
