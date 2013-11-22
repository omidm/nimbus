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
  * Nimbus scheduler.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/scheduler.h"

namespace nimbus {

Scheduler::Scheduler(port_t p)
: listening_port_(p) {
  appId_ = 0;
  registered_worker_num_ = 0;
  data_manager_ = NULL;
  job_manager_ = NULL;
}

Scheduler::~Scheduler() {
  if (data_manager_ != NULL) {
    delete data_manager_;
  }
}

void Scheduler::Run() {
  Log::dbg_printLine("Running the Scheduler");

  SetupWorkerInterface();
  SetupUserInterface();
  SetupDataManager();
  SetupJobManager();
  id_maker_.Initialize(0);

  SchedulerCoreProcessor();
}

void Scheduler::SchedulerCoreProcessor() {
  // Worker registration phase before starting the main job.
  RegisterInitialWorkers(MIN_WORKERS_TO_JOIN);

  // Adding main job to the job manager.
  AddMainJob();

  // Main Loop of the scheduler.
  while (true) {
    RegisterPendingWorkers();
    ProcessQueuedSchedulerCommands((size_t)MAX_BATCH_COMMAND_NUM);
    AssignReadyJobs();
  }
}

void Scheduler::ProcessQueuedSchedulerCommands(size_t max_num) {
  SchedulerCommandList storage;
  if (server_->ReceiveCommands(&storage, max_num)) {
    SchedulerCommandList::iterator iter = storage.begin();
    for (; iter != storage.end(); iter++) {
      SchedulerCommand* comm = *iter;
      dbg(DBG_SCHED, "Processing command: %s.\n", comm->toStringWTags().c_str());
      ProcessSchedulerCommand(comm);
      delete comm;
    }
  }
}

void Scheduler::ProcessSchedulerCommand(SchedulerCommand* cm) {
  switch (cm->type()) {
    case SchedulerCommand::SPAWN_COMPUTE_JOB:
      ProcessSpawnComputeJobCommand(reinterpret_cast<SpawnComputeJobCommand*>(cm));
      break;
    case SchedulerCommand::SPAWN_COPY_JOB:
      ProcessSpawnCopyJobCommand(reinterpret_cast<SpawnCopyJobCommand*>(cm));
      break;
    case SchedulerCommand::DEFINE_DATA:
      ProcessDefineDataCommand(reinterpret_cast<DefineDataCommand*>(cm));
      break;
    case SchedulerCommand::HANDSHAKE:
      ProcessHandshakeCommand(reinterpret_cast<HandshakeCommand*>(cm));
      break;
    case SchedulerCommand::JOB_DONE:
      ProcessJobDoneCommand(reinterpret_cast<JobDoneCommand*>(cm));
      break;
    case SchedulerCommand::DEFINE_PARTITION:
      ProcessDefinePartitionCommand(reinterpret_cast<DefinePartitionCommand*>(cm));
      break;
    default:
      dbg(DBG_ERROR, "ERROR: %s have not been implemented in ProcessSchedulerCommand yet.\n",
          cm->toString().c_str());
  }
}

void Scheduler::ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm) {
  job_manager_->AddJobEntry(JOB_COMP,
        cm->job_name(), cm->job_id().elem(),
        cm->read_set(), cm->write_set(),
        cm->before_set(), cm->after_set(),
        cm->parent_job_id().elem(), cm->params());
}

void Scheduler::ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm) {
  std::string job_name = "copyjob";
  IDSet<logical_data_id_t> read_set;
  read_set.insert(cm->from_logical_id().elem());
  IDSet<logical_data_id_t> write_set;
  write_set.insert(cm->to_logical_id().elem());

  job_manager_->AddJobEntry(JOB_COPY,
        job_name, cm->job_id().elem(),
        read_set, write_set,
        cm->before_set(), cm->after_set(),
        cm->parent_job_id().elem(), cm->params());
}

void Scheduler::ProcessDefineDataCommand(DefineDataCommand* cm) {
  data_manager_->AddLogicalObject(cm->logical_data_id().elem(),
                                 cm->data_name(),
                                 cm->partition_id().elem());
  job_manager_->DefineData(cm->parent_job_id().elem(),
                          cm->logical_data_id().elem());
}

void Scheduler::ProcessDefinePartitionCommand(DefinePartitionCommand* cm) {
  GeometricRegion r = *(cm->region());
  data_manager_->AddPartition(cm->partition_id().elem(), r);
}

void Scheduler::ProcessHandshakeCommand(HandshakeCommand* cm) {
  SchedulerWorkerList::iterator iter = server_->workers()->begin();
  for (; iter != server_->workers()->end(); iter++) {
    if ((*iter)->worker_id() == cm->worker_id().elem()) {
      if ((*iter)->handshake_done()) {
        dbg(DBG_SCHED, "Worker already registered, id: %lu IP: %s port: %lu.\n",
            (*iter)->worker_id(), (*iter)->ip().c_str(), (*iter)->port());
      } else {
        // std::string ip =
        //   (*iter)->connection()->socket()->remote_endpoint().address().to_string();
        (*iter)->set_ip(cm->ip());
        (*iter)->set_port(cm->port().elem());
        (*iter)->set_handshake_done(true);
        ++registered_worker_num_;
        dbg(DBG_SCHED, "Registered new worker, id: %lu IP: %s port: %lu.\n",
            (*iter)->worker_id(), (*iter)->ip().c_str(), (*iter)->port());
      }
      break;
    }
  }
}

void Scheduler::ProcessJobDoneCommand(JobDoneCommand* cm) {
  job_manager_->JobDone(cm->job_id().elem());

  SchedulerWorkerList::iterator iter = server_->workers()->begin();
  for (; iter != server_->workers()->end(); iter++) {
    server_->SendCommand(*iter, cm);
  }
}

void Scheduler::AddMainJob() {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  job_manager_->AddJobEntry(JOB_COMP, "main", j[0], (job_id_t)(0));
}

bool Scheduler::GetWorkerToAssignJob(JobEntry* job, SchedulerWorker*& worker) {
  // Assumption is that partition Ids start from 0, and incrementally go up.
  size_t worker_num = server_->worker_num();
  size_t chunk = (data_manager_->max_defined_partition() + 1) / worker_num;
  std::vector<int> workers_rank(worker_num, 0);

  IDSet<logical_data_id_t> union_set = job->union_set();
  IDSet<logical_data_id_t>::IDSetIter iter;
  for (iter = union_set.begin(); iter != union_set.end(); ++iter) {
    const LogicalDataObject* ldo;
    ldo = data_manager_->FindLogicalObject(*iter);
    size_t poll = std::min((size_t)(ldo->partition()) / chunk, worker_num - 1);
    workers_rank[poll] = workers_rank[poll] + 1;
  }

  // find the worker that wins the poll.
  worker_id_t w_id = 1;
  int count = workers_rank[0];
  for (size_t i = 1; i < worker_num; ++i) {
    if (count < workers_rank[i]) {
      count = workers_rank[i];
      w_id = i + 1;
    }
  }

  return server_->GetSchedulerWorkerById(worker, w_id);
}

bool Scheduler::PrepareDataForJobAtWorker(JobEntry* job,
    SchedulerWorker* worker, logical_data_id_t l_id) {
  LogicalDataObject* ldo =
    const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(l_id));
  JobEntry::VersionTable vt = job->version_table();
  JobEntry::PhysicalTable pt;
  PhysicalDataVector pv;
  data_manager_->InstancesByWorkerAndVersion(ldo, worker->worker_id(), vt[l_id], &pv);
  if (pv.size() > 1) {
    PhysicalData p = pv[0];
    data_manager_->RemovePhysicalInstance(ldo, p);
    p.set_version(p.version() + 1);
    if (job->read_set().contains(l_id))
      p.set_last_job_read(job->job_id());
    if (job->write_set().contains(l_id))
      p.set_last_job_write(job->job_id());
    data_manager_->AddPhysicalInstance(ldo, pv[0]);
    pt[l_id] = p.id();
  } else if (pv.size() == 1) {
    // TODO(omidm): under progress.
    // ...
    // ...
  } else {
    // TODO(omidm): under progress.
    // ...
    // ...
  }

  // TODO(omidm): under progress.
  return false;
}

void Scheduler::SendJobToWorker(JobEntry* job, SchedulerWorker* worker) {
  if (job->job_type() == JOB_COMP) {
    ID<job_id_t> id(job->job_id());
    IDSet<physical_data_id_t> read_set, write_set;
    // TODO(omidm): check the return value of the following methods.
    job->GetPhysicalReadSet(&read_set);
    job->GetPhysicalWriteSet(&write_set);
    ComputeJobCommand cm(job->job_name(), id,
        read_set, write_set, job->before_set(), job->after_set(), job->params());
    dbg(DBG_SCHED, "Sending job %lu to worker %lu: ", job->job_id(), worker->worker_id());
    server_->SendCommand(*(server_->workers()->begin()), &cm);
  } else {
    // TODO(omidm): under progress.
    // ...
    // ...
  }
}

bool Scheduler::AssignJob(JobEntry* job) {
  SchedulerWorker* worker;
  GetWorkerToAssignJob(job, worker);

  IDSet<logical_data_id_t> union_set = job->union_set();
  IDSet<logical_data_id_t>::IDSetIter it;
  for (it = union_set.begin(); it != union_set.end(); ++it) {
    logical_data_id_t id = *it;
    PrepareDataForJobAtWorker(job, worker, id);
    SendJobToWorker(job, worker);
  }
  return false;
}

size_t Scheduler::AssignReadyJobs() {
  size_t count = 0;
  JobEntryList list;
  job_manager_->GetJobsReadyToAssign(&list, (size_t)(MAX_JOB_TO_ASSIGN));
  JobEntryList::iterator iter;
  for (iter = list.begin(); iter != list.end(); ++iter) {
    JobEntry* job = *iter;
    if (AssignJob(job))
      ++count;
  }
  return 0;
}

void Scheduler::RegisterInitialWorkers(size_t min_to_join) {
  while (registered_worker_num_ < min_to_join) {
    dbg(DBG_SCHED, "%d worker(s) are registered, waiting for %d more workers to join ...",
        registered_worker_num_, min_to_join - registered_worker_num_);
    ProcessQueuedSchedulerCommands((size_t)MAX_BATCH_COMMAND_NUM);
    RegisterPendingWorkers();
    sleep(1);
  }
}

size_t Scheduler::RegisterPendingWorkers() {
  size_t registered_num = 0;
  SchedulerWorkerList::iterator iter;
  for (iter = server_->workers()->begin();
      iter != server_->workers()->end(); iter++) {
    if (!(*iter)->handshake_done()) {
      ++registered_num;
      ID<worker_id_t> worker_id((*iter)->worker_id());
      std::string ip("you-know");
      ID<port_t> port(0);
      HandshakeCommand cm(worker_id, ip, port);
      dbg(DBG_SCHED, "Sending command: %s.\n", cm.toStringWTags().c_str());
      server_->SendCommand(*iter, &cm);
    }
  }
  return registered_num;
}

void Scheduler::SetupWorkerInterface() {
  LoadWorkerCommands();
  server_ = new SchedulerServer(listening_port_);
  server_->set_worker_command_table(&worker_command_table_);
  worker_interface_thread_ = new boost::thread(boost::bind(&SchedulerServer::Run, server_));
}

void Scheduler::SetupUserInterface() {
  LoadUserCommands();
  user_interface_thread_ = new boost::thread(
      boost::bind(&Scheduler::GetUserCommand, this));
}

void Scheduler::SetupDataManager() {
  data_manager_ = new DataManager(server_);
}

void Scheduler::SetupJobManager() {
  job_manager_ = new JobManager();
}

void Scheduler::LoadWorkerCommands() {
  // std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  worker_command_table_.push_back(new SpawnComputeJobCommand());
  worker_command_table_.push_back(new SpawnCopyJobCommand());
  worker_command_table_.push_back(new DefineDataCommand());
  worker_command_table_.push_back(new HandshakeCommand());
  worker_command_table_.push_back(new JobDoneCommand());
  worker_command_table_.push_back(new DefinePartitionCommand());
}

void Scheduler::LoadUserCommands() {
  std::stringstream cms("loadapp runapp killapp haltapp resumeapp quit");
  while (true) {
    std::string word;
    cms >> word;
    if (cms.fail()) {
      break;
    }
    user_command_set_.insert(word);
  }
}

void Scheduler::GetUserCommand() {
  while (true) {
    std::cout << "command: " << std::endl;
    std::string token("runapp");
    std::string str, cm;
    std::vector<int> args;
    getline(std::cin, str);
    parseCommand(str, user_command_set_, cm, args);
    std::cout << "you typed: " << cm << std::endl;
  }
}

}  // namespace nimbus
