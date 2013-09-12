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

#include "worker/worker.h"

#define MAX_PARALLEL_JOB 10

using namespace nimbus; // NOLINT

Worker::Worker(std::string scheduler_ip, port_t scheduler_port,
    port_t listening_port, Application* a)
: scheduler_ip_(scheduler_ip),
  scheduler_port_(scheduler_port),
  listening_port_(listening_port),
  application_(a) {
    log.InitTime();
    id_ = -1;
}

void Worker::Run() {
  std::cout << "Running the Worker" << std::endl;

  SetupSchedulerInterface();
  application_->Start(client_, &id_maker_);

  SetupDataExchangerInterface();

  WorkerCoreProcessor();
}

void Worker::WorkerCoreProcessor() {
  std::cout << "Base Worker Core Processor" << std::endl;

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

    JobList jobs_to_run;
    GetJobsToRun(&jobs_to_run, (size_t)(MAX_PARALLEL_JOB));
    JobList::iterator iter = jobs_to_run.begin();
    for (; iter != jobs_to_run.end(); iter++) {
      ExecuteJob(*iter);
    }
  }
}

void Worker::ScanBlockedJobs() {
  JobList::iterator iter;
  for (iter = blocked_jobs_.begin(); iter != blocked_jobs_.end();) {
    IDSet<job_id_t> req = (*iter)->before_set();
    IDSet<job_id_t>::IDSetIter it;
    for (it = req.begin(); it != req.end();) {
      if (done_jobs_.count(*it) != 0)
        req.remove(it++);
      else
        ++it;
    }
    (*iter)->set_before_set(req);
    if ((*iter)->before_set().size() == 0) {
      if (dynamic_cast<RemoteCopyReceiveJob*>(*iter) == NULL) // NOLINT
        ready_jobs_.push_back(*iter);
      else
        pending_transfer_jobs_.push_back(*iter);
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
    if (data_exchanger_->ReceiveSerializedData((*iter)->id().elem(), &ser_data)) {
      static_cast<RemoteCopyReceiveJob*>(*iter)->set_serialized_data(ser_data);
      ready_jobs_.push_back(*iter);
      pending_transfer_jobs_.erase(iter++);
    } else {
      ++iter;
    }
  }
}

void Worker::GetJobsToRun(JobList* list, size_t max_num) {
  list->clear();
  size_t ready_num = ready_jobs_.size();
  for (size_t i = 0; (i < max_num) && (i < ready_num); i++) {
    Job* job = ready_jobs_.front();
    list->push_back(job);
    ready_jobs_.pop_front();
  }
}

void Worker::ExecuteJob(Job* job) {
  DataArray da;
  IDSet<data_id_t>::IDSetIter iter;

  IDSet<data_id_t> read = job->read_set();
  for (iter = read.begin(); iter != read.end(); iter++)
    da.push_back(data_map_[*iter]);

  IDSet<data_id_t> write = job->write_set();
  for (iter = write.begin(); iter != write.end(); iter++)
    da.push_back(data_map_[*iter]);

  log.StartTimer();
  job->Execute(job->parameters(), da);
  log.StopTimer();

  char buff[MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff),
      "Execute Job, name: %25s  id: %4ld  length(ms): %6.3lf  time(s): %6.3lf",
      job->name().c_str(), job->id().elem(), 1000 * log.timer(), log.GetTime());

  log.writeToFile(std::string(buff), LOG_INFO);

  JobDoneCommand cm(job->id(), job->after_set(), "nothing_yet");
  client_->sendCommand(&cm);
}




void Worker::ProcessSchedulerCommand(SchedulerCommand* cm) {
  std::string command_name = cm->name();

  if (command_name == "handshake") {
    ProcessHandshakeCommand(reinterpret_cast<HandshakeCommand*>(cm));
  } else if (command_name == "jobdone") {
    ProcessJobDoneCommand(reinterpret_cast<JobDoneCommand*>(cm));
  } else if (command_name == "computejob") {
    ProcessComputeJobCommand(reinterpret_cast<ComputeJobCommand*>(cm));
  } else if (command_name == "createdata") {
    ProcessCreateDataCommand(reinterpret_cast<CreateDataCommand*>(cm));
  } else if (command_name == "remotecopysend") {
    ProcessRemoteCopySendCommand(reinterpret_cast<RemoteCopySendCommand*>(cm));
  } else if (command_name == "remotecopyreceive") {
    ProcessRemoteCopyReceiveCommand(reinterpret_cast<RemoteCopyReceiveCommand*>(cm));
  } else if (command_name == "localcopy") {
    ProcessLocalCopyCommand(reinterpret_cast<LocalCopyCommand*>(cm));
  } else if (command_name == "spawnjob") {
    ProcessSpawnJobCommand(reinterpret_cast<SpawnJobCommand*>(cm));
  } else if (command_name == "definedata") {
    ProcessDefineDataCommand(reinterpret_cast<DefineDataCommand*>(cm));
  } else {
    std::cout << "ERROR: " << cm->toString() <<
      " have not been implemented in ProcessSchedulerCommand yet." <<
      std::endl;
  }
}

void Worker::ProcessHandshakeCommand(HandshakeCommand* cm) {
  ID<port_t> port(listening_port_);
  HandshakeCommand new_cm = HandshakeCommand(cm->worker_id(),
      boost::asio::ip::host_name(), port);
  client_->sendCommand(&new_cm);

  id_ = cm->worker_id().elem();
  id_maker_.Initialize(id_);
}

void Worker::ProcessJobDoneCommand(JobDoneCommand* cm) {
  std::map<job_id_t, IDSet<job_id_t> >::iterator iter;
  for (iter = done_jobs_.begin(); iter != done_jobs_.end();) {
    iter->second.remove(cm->job_id().elem());
    if (iter->second.size() == 0)
      done_jobs_.erase(iter++);
    else
      ++iter;
  }
  if (cm->after_set().size() != 0)
    done_jobs_[cm->job_id().elem()] = cm->after_set();
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
  blocked_jobs_.push_back(job);
}

void Worker::ProcessCreateDataCommand(CreateDataCommand* cm) {
  Data * data = application_->CloneData(cm->data_name());
  data->set_id(cm->data_id().elem());
  AddData(data);

  Job * job = new CreateDataJob();
  job->set_name("CreateData:" + cm->data_name());
  job->set_id(cm->job_id());
  IDSet<data_id_t> write_set;
  write_set.insert(cm->data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
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
  IDSet<data_id_t> read_set;
  read_set.insert(cm->from_data_id().elem());
  job->set_read_set(read_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  blocked_jobs_.push_back(job);
}

void Worker::ProcessRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* cm) {
  Job * job = new RemoteCopyReceiveJob();
  job->set_name("RemoteCopyReceive");
  job->set_id(cm->job_id());
  IDSet<data_id_t> write_set;
  write_set.insert(cm->to_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  blocked_jobs_.push_back(job);
}

void Worker::ProcessLocalCopyCommand(LocalCopyCommand* cm) {
  Job * job = new LocalCopyJob();
  job->set_name("LocalCopy");
  job->set_id(cm->job_id());
  IDSet<data_id_t> read_set;
  read_set.insert(cm->from_data_id().elem());
  job->set_read_set(read_set);
  IDSet<data_id_t> write_set;
  write_set.insert(cm->to_data_id().elem());
  job->set_write_set(write_set);
  job->set_before_set(cm->before_set());
  job->set_after_set(cm->after_set());
  blocked_jobs_.push_back(job);
}

void Worker::ProcessSpawnJobCommand(SpawnJobCommand* cm) {
  Job * j = application_->CloneJob(cm->job_name());

  std::vector<Data*> da;
  IDSet<data_id_t>::IDSetIter iter;

  IDSet<data_id_t> read = cm->read_set();
  for (iter = read.begin(); iter != read.end(); iter++)
    da.push_back(data_map_[*iter]);

  IDSet<data_id_t> write = cm->write_set();
  for (iter = write.begin(); iter != write.end(); iter++)
    da.push_back(data_map_[*iter]);

  job_id_t id = *(cm->job_id().begin());

  log.StartTimer();
  j->Execute(cm->params(), da);
  log.StopTimer();

  char buff[MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff),
      "Execute Job, name: %25s  id: %4ld  length(ms): %6.3lf  time(s): %6.3lf",
      cm->job_name().c_str(), id, 1000 * log.timer(), log.GetTime());

  log.writeToFile(std::string(buff), LOG_INFO);
}

void Worker::ProcessDefineDataCommand(DefineDataCommand* cm) {
  Data * d = application_->CloneData(cm->data_name());

  data_id_t id = cm->data_id().elem();

  log.StartTimer();
  d->Create();
  d->set_id(id);
  AddData(d);
  log.StopTimer();

  char buff[MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff),
      "Create Data, name: %25s  id: %4ld  length(ms): %6.3lf  time(s): %6.3lf",
      cm->data_name().c_str(), id, 1000 * log.timer(), log.GetTime());

  log.writeToFile(std::string(buff), LOG_INFO);
}

void Worker::SetupDataExchangerInterface() {
  data_exchanger_ = new WorkerDataExchanger(listening_port_);
  data_exchanger_thread_ = new boost::thread(
      boost::bind(&WorkerDataExchanger::Run, data_exchanger_));
}

void Worker::SetupSchedulerInterface() {
  LoadSchedulerCommands();
  client_ = new SchedulerClient(scheduler_ip_, scheduler_port_);
  client_->set_scheduler_command_set(&scheduler_command_set_);
  client_->run();
  // client_thread_ = new boost::thread(
  //     boost::bind(&SchedulerClient::run, client_));
}

void Worker::LoadSchedulerCommands() {
  // std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  scheduler_command_set_.insert(
      std::make_pair(std::string("spawnjob"), COMMAND_SPAWN_JOB));
  scheduler_command_set_.insert(
      std::make_pair(std::string("definedata"), COMMAND_DEFINE_DATA));
  scheduler_command_set_.insert(
      std::make_pair(std::string("handshake"), COMMAND_HANDSHAKE));
  scheduler_command_set_.insert(
      std::make_pair(std::string("jobdone"), COMMAND_JOB_DONE));
  scheduler_command_set_.insert(
      std::make_pair(std::string("computejob"), COMMAND_COMPUTE_JOB));
  scheduler_command_set_.insert(
      std::make_pair(std::string("createdata"), COMMAND_CREATE_DATA));
  scheduler_command_set_.insert(
      std::make_pair(std::string("remotecopysend"), COMMAND_REMOTE_COPY_SEND));
  scheduler_command_set_.insert(
      std::make_pair(std::string("remotecopyreceive"), COMMAND_REMOTE_COPY_RECEIVE));
  scheduler_command_set_.insert(
      std::make_pair(std::string("localcopy"), COMMAND_LOCAL_COPY));
}

void Worker::AddData(Data* data) {
  data_map_[data->id()] = data;
}

void Worker::DeleteData(data_id_t data_id) {
  data_map_.erase(data_id);
}

worker_id_t Worker::id() {
  return id_;
}

void Worker::set_id(worker_id_t id) {
  id_ = id;
}


