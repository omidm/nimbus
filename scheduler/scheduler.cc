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

#define MAX_BATCH_COMMAND_NUM 10000
#define DEFAULT_MIN_WORKER_TO_JOIN 2
#define JOB_ASSIGNER_THREAD_NUM 0
#define MAX_JOB_TO_ASSIGN 1

Scheduler::Scheduler(port_t p)
: listening_port_(p) {
  appId_ = 0;
  registered_worker_num_ = 0;
  data_manager_ = NULL;
  job_manager_ = NULL;
  min_worker_to_join_ = DEFAULT_MIN_WORKER_TO_JOIN;
  terminate_application_flag_ = false;
  stamp_state_ = -1;
  compute_count_ = 0;
  copy_count_ = 0;
}

Scheduler::~Scheduler() {
  if (data_manager_ != NULL) {
    delete data_manager_;
  }
  if (job_manager_ != NULL) {
    delete job_manager_;
  }
}

void Scheduler::Run() {
  Log::log_PrintLine("Running the Scheduler");

  /*
   * First create the modules, then set them up.
   */

  CreateIDMaker();
  CreateSchedulerServer();
  CreateDataManager();
  CreateJobManager();
  CreateLoadBalancer();
  CreateJobAssigner();

  SetupIDMaker();
  SetupSchedulerServer();
  SetupDataManager();
  SetupJobManager();
  SetupLoadBalancer();
  SetupJobAssigner();
  // SetupUserInterface();

  SchedulerCoreProcessor();
}


void Scheduler::set_min_worker_to_join(size_t num) {
  min_worker_to_join_ = num;
}

void Scheduler::SchedulerCoreProcessor() {
  // Worker registration phase before starting the main job.
  RegisterInitialWorkers(min_worker_to_join_);

  // Adding main job to the job manager.
  AddMainJob();

  // Main Loop of the scheduler.
  while (true) {
    log_loop_.log_StartTimer();
    log_assign_.log_ResetTimer();
    log_server_.log_ResetTimer();
    log_job_manager_.log_ResetTimer();
    log_data_manager_.log_ResetTimer();
    log_version_manager_.log_ResetTimer();
    log_load_balancer_.log_ResetTimer();

    RegisterPendingWorkers();
    ProcessQueuedSchedulerCommands((size_t)MAX_BATCH_COMMAND_NUM);
    AssignReadyJobs();
    RemoveObsoleteJobEntries();
    TerminationProcedure();

    log_loop_.log_StopTimer();
    if (log_loop_.timer() >= .001) {
      char buff[LOG_MAX_BUFF_SIZE];
      snprintf(buff, sizeof(buff),
          "loop: %2.5lf  assign: %2.5lf server: %2.5lf job_manager: %2.5lf data_manager: %2.5lf version_manager: %2.5lf load_balancer: %2.5lf compute_count: %d  copy_count: %d ldo_count: %d time: %6.5lf", // NOLINT
          log_loop_.timer(),
          log_assign_.timer(),
          log_server_.timer(),
          log_job_manager_.timer(),
          log_data_manager_.timer(),
          log_version_manager_.timer(),
          log_load_balancer_.timer(),
          compute_count_,
          copy_count_,
          data_manager_->ldo_map_p()->size(),
          log_loop_.GetTime());

      log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);
    }
  }
}

void Scheduler::ProcessQueuedSchedulerCommands(size_t max_num) {
  SchedulerCommandList storage;
  if (server_->ReceiveCommands(&storage, max_num)) {
    SchedulerCommandList::iterator iter = storage.begin();
    for (; iter != storage.end(); iter++) {
      SchedulerCommand* comm = *iter;
      dbg(DBG_SCHED, "Processing command: %s.\n", comm->ToString().c_str());
      ProcessSchedulerCommand(comm);
      delete comm;
    }
  }
}

void Scheduler::ProcessSchedulerCommand(SchedulerCommand* cm) {
  switch (cm->type()) {
    case SchedulerCommand::SPAWN_COMPUTE:
      ProcessSpawnComputeJobCommand(reinterpret_cast<SpawnComputeJobCommand*>(cm));
      break;
    case SchedulerCommand::SPAWN_COPY:
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
    case SchedulerCommand::TERMINATE:
      ProcessTerminateCommand(reinterpret_cast<TerminateCommand*>(cm));
      break;
    case SchedulerCommand::PROFILE:
      break;
    default:
      dbg(DBG_ERROR, "ERROR: %s have not been implemented in ProcessSchedulerCommand yet.\n",
          cm->ToNetworkData().c_str());
  }
}

void Scheduler::ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm) {
  log_job_manager_.log_ResumeTimer();
  job_manager_->AddComputeJobEntry(cm->job_name(),
                                   cm->job_id().elem(),
                                   cm->read_set(),
                                   cm->write_set(),
                                   cm->before_set(),
                                   cm->after_set(),
                                   cm->parent_job_id().elem(),
                                   cm->future_job_id().elem(),
                                   cm->sterile(),
                                   cm->params());
  log_job_manager_.log_StopTimer();
}

void Scheduler::ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm) {
  std::string job_name = "copyjob";
  IDSet<logical_data_id_t> read_set;
  read_set.insert(cm->from_logical_id().elem());
  IDSet<logical_data_id_t> write_set;
  write_set.insert(cm->to_logical_id().elem());

  // TODO(omid): we need to add support for copy jobs.
  log_job_manager_.log_ResumeTimer();
  job_manager_->AddExplicitCopyJobEntry();
  log_job_manager_.log_StopTimer();
}

void Scheduler::ProcessDefineDataCommand(DefineDataCommand* cm) {
  log_data_manager_.log_ResumeTimer();
  bool success = data_manager_->AddLogicalObject(cm->logical_data_id().elem(),
                                 cm->data_name(),
                                 cm->partition_id().elem());
  log_data_manager_.log_StopTimer();

  if (success) {
    log_data_manager_.log_ResumeTimer();
    LdoAddCommand command(data_manager_->FindLogicalObject(cm->logical_data_id().elem()));
    log_data_manager_.log_StopTimer();
    server_->BroadcastCommand(&command);
  }

  log_job_manager_.log_ResumeTimer();
  job_manager_->DefineData(cm->parent_job_id().elem(),
                          cm->logical_data_id().elem());
  log_job_manager_.log_StopTimer();
}

void Scheduler::ProcessDefinePartitionCommand(DefinePartitionCommand* cm) {
  GeometricRegion r = *(cm->region());
  log_data_manager_.log_ResumeTimer();
  data_manager_->AddPartition(cm->partition_id().elem(), r);
  log_data_manager_.log_StopTimer();
  PartitionAddCommand command(cm->partition_id(), r);
  server_->BroadcastCommand(&command);
}

void Scheduler::ProcessHandshakeCommand(HandshakeCommand* cm) {
  SchedulerWorkerList::iterator iter = server_->workers()->begin();
  for (; iter != server_->workers()->end(); iter++) {
    if ((*iter)->worker_id() == cm->worker_id().elem()) {
      if ((*iter)->handshake_done()) {
        dbg(DBG_SCHED, "Worker already registered, id: %lu IP: %s port: %lu.\n",
            (*iter)->worker_id(), (*iter)->ip().c_str(), (*iter)->port());
      } else {
        std::string ip;
        if (cm->ip() == NIMBUS_RECEIVER_KNOWN_IP) {
          ip =
            (*iter)->connection()->socket()->remote_endpoint().address().to_string();
        } else {
          ip = cm->ip();
        }
        (*iter)->set_ip(ip);
        (*iter)->set_port(cm->port().elem());
        (*iter)->set_handshake_done(true);
        ++registered_worker_num_;
        dbg(DBG_SCHED, "Registered new worker, id: %lu IP: %s port: %lu.\n",
            (*iter)->worker_id(), (*iter)->ip().c_str(), (*iter)->port());
        log_load_balancer_.log_ResumeTimer();
        load_balancer_->NotifyRegisteredWorker(*iter);
        log_load_balancer_.log_StopTimer();
      }
      break;
    }
  }
}

void Scheduler::ProcessJobDoneCommand(JobDoneCommand* cm) {
  job_id_t job_id = cm->job_id().elem();

  JobEntry *job;
  log_job_manager_.log_ResumeTimer();
  if (job_manager_->GetJobEntry(job_id, job)) {
    job_manager_->NotifyJobDone(job);
    log_load_balancer_.log_ResumeTimer();
    load_balancer_->NotifyJobDone(job);
    log_load_balancer_.log_StopTimer();
  }
  log_job_manager_.log_StopTimer();

  std::string jname = job->job_name();
  if (jname == "loop_iteration") {
    log_.log_StartTimer();
    stamp_state_ = 0;
    compute_count_ = 0;
    copy_count_ = 0;
  }

  SchedulerWorkerList::iterator iter = server_->workers()->begin();
  for (; iter != server_->workers()->end(); iter++) {
    server_->SendCommand(*iter, cm);
  }
}

void Scheduler::ProcessTerminateCommand(TerminateCommand* cm) {
  terminate_application_flag_ = true;
  terminate_application_status_ = cm->exit_status().elem();
}

void Scheduler::TerminationProcedure() {
  if (terminate_application_flag_) {
    if (job_manager_->AllJobsAreDone()) {
      SchedulerCommand* command =
        new TerminateCommand(ID<exit_status_t>(terminate_application_status_));
      server_->BroadcastCommand(command);
      delete command;
      exit(NIMBUS_TERMINATE_SUCCESS);
    }
  }
}

void Scheduler::AddMainJob() {
  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 1);
  log_job_manager_.log_ResumeTimer();
  job_manager_->AddMainJobEntry(j[0]);
  log_job_manager_.log_StopTimer();
}

size_t Scheduler::RemoveObsoleteJobEntries() {
  log_job_manager_.log_ResumeTimer();
  size_t count = job_manager_->RemoveObsoleteJobEntries();
  log_job_manager_.log_StopTimer();

  return count;
}

size_t Scheduler::AssignReadyJobs() {
  return load_balancer_->AssignReadyJobs();
}

void Scheduler::RegisterInitialWorkers(size_t min_to_join) {
  while (registered_worker_num_ < min_to_join) {
    dbg(DBG_SCHED, "%d worker(s) are registered, waiting for %d more workers to join ...\n",
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
      dbg(DBG_SCHED, "Sending command: %s.\n", cm.ToString().c_str());
      server_->SendCommand(*iter, &cm);
    }
  }
  return registered_num;
}

void Scheduler::CreateIDMaker() {
  id_maker_ = new IDMaker();
}

void Scheduler::CreateSchedulerServer() {
  server_ = new SchedulerServer(listening_port_);
}

void Scheduler::CreateDataManager() {
  data_manager_ = new DataManager();
}

void Scheduler::CreateJobManager() {
  job_manager_ = new JobManager();
}

void Scheduler::CreateLoadBalancer() {
  load_balancer_ = new LoadBalancer();
}

void Scheduler::CreateJobAssigner() {
  job_assigner_ = new JobAssigner();
}

void Scheduler::SetupIDMaker() {
  id_maker_->Initialize(0);
}

void Scheduler::SetupSchedulerServer() {
  LoadWorkerCommands();
  server_->set_worker_command_table(&worker_command_table_);
  scheduler_server_thread_ = new boost::thread(boost::bind(&SchedulerServer::Run, server_));
}

void Scheduler::SetupDataManager() {
}

void Scheduler::SetupJobManager() {
  job_manager_->set_ldo_map_p(data_manager_->ldo_map_p());
}

void Scheduler::SetupLoadBalancer() {
  load_balancer_->set_job_manager(job_manager_);
  load_balancer_->set_data_manager(data_manager_);
  load_balancer_->set_job_assigner(job_assigner_);
  load_balancer_->set_max_job_to_assign(MAX_JOB_TO_ASSIGN);
  load_balancer_thread_ = new boost::thread(boost::bind(&LoadBalancer::Run, load_balancer_));
}

void Scheduler::SetupJobAssigner() {
  job_assigner_->set_id_maker(id_maker_);
  job_assigner_->set_server(server_);
  job_assigner_->set_job_manager(job_manager_);
  job_assigner_->set_data_manager(data_manager_);
  job_assigner_->set_load_balancer(load_balancer_);
  job_assigner_->set_thread_num(JOB_ASSIGNER_THREAD_NUM);
  job_assigner_->Run();
}

void Scheduler::SetupUserInterface() {
  LoadUserCommands();
  user_interface_thread_ = new boost::thread(
      boost::bind(&Scheduler::GetUserCommand, this));
}

void Scheduler::LoadWorkerCommands() {
  worker_command_table_[SchedulerCommand::SPAWN_COMPUTE]     = new SpawnComputeJobCommand();
  worker_command_table_[SchedulerCommand::SPAWN_COPY]        = new SpawnCopyJobCommand();
  worker_command_table_[SchedulerCommand::DEFINE_DATA]       = new DefineDataCommand();
  worker_command_table_[SchedulerCommand::HANDSHAKE]         = new HandshakeCommand();
  worker_command_table_[SchedulerCommand::JOB_DONE]          = new JobDoneCommand();
  worker_command_table_[SchedulerCommand::DEFINE_PARTITION]  = new DefinePartitionCommand();
  worker_command_table_[SchedulerCommand::TERMINATE]         = new TerminateCommand();
  worker_command_table_[SchedulerCommand::PROFILE]           = new ProfileCommand();
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
