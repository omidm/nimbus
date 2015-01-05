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

#define HANDSHAKE_COMMAND_NUM 1

namespace nimbus {

#define DEFAULT_MAX_JOB_TO_ASSIGN 1
#define DEFAULT_MAX_JOB_TO_REMOVE 1
#define DEFAULT_MIN_WORKER_TO_JOIN 2
#define DEFAULT_JOB_ASSIGNER_THREAD_NUM 0
#define DEFAULT_MAX_COMMAND_PROCESS_NUM 10000
#define DEFAULT_MAX_JOB_DONE_COMMAND_PROCESS_NUM 10000

Scheduler::Scheduler(port_t port) {
  server_ = NULL;
  id_maker_ = NULL;
  job_manager_ = NULL;
  template_manager_ = NULL;
  data_manager_ = NULL;
  job_assigner_ = NULL;
  load_balancer_ = NULL;
  listening_port_ = port;
  registered_worker_num_ = 0;
  terminate_application_flag_ = false;
  cleaner_thread_active_ = false;
  max_job_to_assign_ = DEFAULT_MAX_JOB_TO_ASSIGN;
  max_job_to_remove_ = DEFAULT_MAX_JOB_TO_REMOVE;
  min_worker_to_join_ = DEFAULT_MIN_WORKER_TO_JOIN;
  job_assigner_thread_num_ = DEFAULT_JOB_ASSIGNER_THREAD_NUM;
  max_command_process_num_ = DEFAULT_MAX_COMMAND_PROCESS_NUM;
  max_job_done_command_process_num_ = DEFAULT_MAX_JOB_DONE_COMMAND_PROCESS_NUM;
  log_.set_file_name("log_scheduler");
  log_receive_stamp_.set_file_name("log_receive_stamp");
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
  CreateAfterMap();
  CreateSchedulerServer();
  CreateDataManager();
  CreateJobManager();
  CreateTemplateManager();
  CreateLoadBalancer();
  CreateJobAssigner();

  SetupIDMaker();
  SetupSchedulerServer();
  SetupDataManager();
  SetupJobManager();
  SetupTemplateManager();
  SetupLoadBalancer();
  SetupJobAssigner();
  SetupJobDoneBouncer();
  SetupCleaner();
  // SetupUserInterface();

  SchedulerCoreProcessor();
}

void Scheduler::set_max_job_to_remove(size_t num) {
  max_job_to_remove_ = num;
}

void Scheduler::set_max_job_to_assign(size_t num) {
  max_job_to_assign_ = num;
}

void Scheduler::set_min_worker_to_join(size_t num) {
  min_worker_to_join_ = num;
}

void Scheduler::set_job_assigner_thread_num(size_t num) {
  job_assigner_thread_num_ = num;
}

void Scheduler::set_max_command_process_num(size_t num) {
  max_command_process_num_ = num;
}

void Scheduler::set_max_job_done_command_process_num(size_t num) {
  max_job_done_command_process_num_ = num;
}

void Scheduler::SchedulerCoreProcessor() {
  // Worker registration phase before starting the main job.
  RegisterInitialWorkers();

  // Adding main job to the job manager.
  AddMainJob();

  // Main Loop of the scheduler.
  while (true) {
    log_.log_StartTimer();

    RegisterPendingWorkers();

    log_process_.StartTimer();
    size_t command_num = ProcessQueuedSchedulerCommands();
    log_process_.StopTimer();

    log_assign_.StartTimer();
    AssignReadyJobs();
    log_assign_.StopTimer();

    RemoveObsoleteJobEntries();

    TerminationProcedure();

    log_.log_StopTimer();
    if ((log_.timer() >= .001) || (command_num > 0)) {
      char buff[LOG_MAX_BUFF_SIZE];
      snprintf(buff, sizeof(buff), "%10.9lf l: %2.5lf p: %2.5lf c: %5.0lu a: %2.5lf.",
          Log::GetRawTime(), log_.timer(), log_process_.timer(), command_num, log_assign_.timer());
      log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);
    }
  }
}

size_t Scheduler::ProcessQueuedSchedulerCommands() {
  size_t count = 0;
  SchedulerCommandList storage;
  if (server_->ReceiveCommands(&storage, max_command_process_num_)) {
    bool flush = false;
    SchedulerCommandList::iterator iter = storage.begin();
    for (; iter != storage.end(); iter++) {
      SchedulerCommand* comm = *iter;
      SchedulerCommand::Type type = comm->type();

      if (!flush || (type == SchedulerCommand::WORKER_DOWN)) {
        dbg(DBG_SCHED, "Processing command: %s.\n", comm->ToString().c_str());
        ProcessSchedulerCommand(comm);
      } else {
        dbg(DBG_SCHED, "Flushed command: %s.\n", comm->ToString().c_str());
      }

      if (type == SchedulerCommand::WORKER_DOWN) {
        flush = true;
      }

      delete comm;
      ++count;
    }
  }
  return count;
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
    case SchedulerCommand::SPAWN_TEMPLATE:
      ProcessSpawnTemplateCommand(reinterpret_cast<SpawnTemplateCommand*>(cm));
      break;
    case SchedulerCommand::START_TEMPLATE:
      ProcessStartTemplateCommand(reinterpret_cast<StartTemplateCommand*>(cm));
      break;
    case SchedulerCommand::END_TEMPLATE:
      ProcessEndTemplateCommand(reinterpret_cast<EndTemplateCommand*>(cm));
      break;
    case SchedulerCommand::SAVE_DATA_JOB_DONE:
      ProcessSaveDataJobDoneCommand(reinterpret_cast<SaveDataJobDoneCommand*>(cm));
      break;
    case SchedulerCommand::WORKER_DOWN:
      ProcessWorkerDownCommand(reinterpret_cast<WorkerDownCommand*>(cm));
      break;
    case SchedulerCommand::PREPARE_REWIND:
      ProcessPrepareRewindCommand(reinterpret_cast<PrepareRewindCommand*>(cm));
      break;
    default:
      dbg(DBG_ERROR, "ERROR: %s have not been implemented in ProcessSchedulerCommand yet.\n",
          cm->ToNetworkData().c_str());
      exit(-1);
  }
}

void Scheduler::ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm) {
  char buff[LOG_MAX_BUFF_SIZE];

  snprintf(buff, sizeof(buff), "%10.9f JR id: %lu n: %s.",
      Log::GetRawTime(), cm->job_id().elem(), cm->job_name().c_str());
  log_receive_stamp_.log_WriteToFile(std::string(buff));

  job_manager_->AddComputeJobEntry(cm->job_name(),
                                   cm->job_id().elem(),
                                   cm->read_set(),
                                   cm->write_set(),
                                   cm->before_set(),
                                   cm->after_set(),
                                   cm->parent_job_id().elem(),
                                   cm->future_job_id().elem(),
                                   cm->sterile(),
                                   cm->region(),
                                   cm->params());

  std::map<job_id_t, std::string>::iterator iter =
    template_spawner_map_.find(cm->parent_job_id().elem());
  if (iter != template_spawner_map_.end()) {
    template_manager_->AddComputeJobToTemplate(iter->second,
                                               cm->job_name(),
                                               cm->job_id().elem(),
                                               cm->read_set(),
                                               cm->write_set(),
                                               cm->before_set(),
                                               cm->after_set(),
                                               cm->parent_job_id().elem(),
                                               cm->future_job_id().elem(),
                                               cm->sterile(),
                                               cm->region());
  }

  snprintf(buff, sizeof(buff), "%10.9f JE id: %lu n: %s.",
      Log::GetRawTime(), cm->job_id().elem(), cm->job_name().c_str());
  log_receive_stamp_.log_WriteToFile(std::string(buff));
}

void Scheduler::ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm) {
  std::string job_name = "copyjob";
  IDSet<logical_data_id_t> read_set;
  read_set.insert(cm->from_logical_id().elem());
  IDSet<logical_data_id_t> write_set;
  write_set.insert(cm->to_logical_id().elem());

  // TODO(omid): we need to add support for copy jobs.
  job_manager_->AddExplicitCopyJobEntry();
}

void Scheduler::ProcessDefineDataCommand(DefineDataCommand* cm) {
  bool success = data_manager_->AddLogicalObject(cm->logical_data_id().elem(),
                                 cm->data_name(),
                                 cm->partition_id().elem());

  if (success) {
    LdoAddCommand command(data_manager_->FindLogicalObject(cm->logical_data_id().elem()));
    server_->BroadcastCommand(&command);
  }

  job_manager_->DefineData(cm->parent_job_id().elem(),
                          cm->logical_data_id().elem());
}

void Scheduler::ProcessDefinePartitionCommand(DefinePartitionCommand* cm) {
  GeometricRegion r = *(cm->region());
  data_manager_->AddPartition(cm->partition_id().elem(), r);
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

        char buff[LOG_MAX_BUFF_SIZE];
        snprintf(buff, sizeof(buff), "Registered new worker, id: %u IP: %s port: %u time: %2.2lf.", // NOLINT
            (*iter)->worker_id(), (*iter)->ip().c_str(), (*iter)->port(), log_.GetTime());
        log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);

        dbg(DBG_SCHED, "Registered new worker, id: %u IP: %s port: %u.\n",
            (*iter)->worker_id(), (*iter)->ip().c_str(), (*iter)->port());

        load_balancer_->NotifyRegisteredWorker(*iter);
      }
      break;
    }
  }
}

void Scheduler::ProcessSaveDataJobDoneCommand(SaveDataJobDoneCommand* cm) {
  job_manager_->NotifySaveDataJobDoneForCheckpoint(cm->checkpoint_id().elem(),
                                                   cm->job_id().elem(),
                                                   cm->handle());
}

void Scheduler::ProcessWorkerDownCommand(WorkerDownCommand* cm) {
  // Cleanup the state and rewind from a checkpoint.
  worker_id_t worker_id = cm->worker_id().elem();

  load_balancer_->NotifyDownWorker(worker_id);

  server_->RemoveWorker(worker_id);

  data_manager_->RemoveAllInstanceByWorker(worker_id);
  data_manager_->ResetAllInstances();

  checkpoint_id_t checkpoint_id = 0;
  job_manager_->RewindFromLastCheckpoint(&checkpoint_id);

  job_assigner_->set_checkpoint_id(checkpoint_id);

  PrepareRewindCommand command =
    PrepareRewindCommand(ID<worker_id_t>(1),
                         ID<checkpoint_id_t>(checkpoint_id));
  server_->BroadcastCommand(&command);
  WaitForAllPrepareRewindResponses();
}

void Scheduler::WaitForAllPrepareRewindResponses() {
  std::set<worker_id_t> pending_workers;
  SchedulerWorkerList::iterator it = server_->workers()->begin();
  for (; it != server_->workers()->end(); ++it) {
    pending_workers.insert((*it)->worker_id());
  }
  assert(pending_workers.size() > 0);

  while (pending_workers.size() > 0) {
    size_t count = 0;
    SchedulerCommandList storage;
    if (server_->ReceiveCommands(&storage, max_command_process_num_)) {
      SchedulerCommandList::iterator iter = storage.begin();
      for (; iter != storage.end(); iter++) {
        SchedulerCommand* comm = *iter;
        dbg(DBG_SCHED, "Wait For All Prepare Rewind Responses: considering command: %s.\n",
            comm->ToString().c_str());
        switch (comm->type()) {
          case SchedulerCommand::PREPARE_REWIND:
            pending_workers.erase(reinterpret_cast<PrepareRewindCommand*>(comm)->worker_id().elem()); // NOLINT
            break;
          default:
            dbg(DBG_WARN, "WARNING: Ignored command %s in rewind period.\n",
                comm->ToString().c_str());
        }
        delete comm;
        ++count;
      }
    }
    if (count == 0) {
      usleep(10);
    }
  }
}

void Scheduler::ProcessPrepareRewindCommand(PrepareRewindCommand* cm) {
  dbg(DBG_WARN, "WARNING: unexpected PrepareRewind command from worker %lu.\n",
      cm->worker_id().elem());
  exit(-1);
}

void Scheduler::ProcessJobDoneCommand(JobDoneCommand* cm) {
  job_id_t job_id = cm->job_id().elem();

  char buff[LOG_MAX_BUFF_SIZE];

  snprintf(buff, sizeof(buff), "%10.9f DR id: %lu.",
      Log::GetRawTime(), job_id);
  log_receive_stamp_.log_WriteToFile(std::string(buff));

  if (!id_maker_->SchedulerProducedJobID(job_id)) {
    std::list<SchedulerWorker*> waiting_list;
    after_map_->GetWorkersWaitingOnJob(job_id, &waiting_list);

    cm->set_final(true);
    std::list<SchedulerWorker*>::iterator iter = waiting_list.begin();
    for (; iter != waiting_list.end(); ++iter) {
      server_->SendCommand(*iter, cm);
    }
  }

  // AfterMap has internal locking.
  after_map_->RemoveJobRecords(job_id);

  JobEntry *job;
  if (job_manager_->GetJobEntry(job_id, job)) {
    job->set_done(true);
    load_balancer_->NotifyJobDone(job);
    job_manager_->NotifyJobDone(job);
  }

  snprintf(buff, sizeof(buff), "%10.9f DE id: %lu.",
      Log::GetRawTime(), job_id);
  log_receive_stamp_.log_WriteToFile(std::string(buff));
}

void Scheduler::ProcessTerminateCommand(TerminateCommand* cm) {
  terminate_application_flag_ = true;
  terminate_application_status_ = cm->exit_status().elem();
}

void Scheduler::ProcessSpawnTemplateCommand(SpawnTemplateCommand* cm) {
  // TODO(omidm): Implement!
  if (NIMBUS_TEMPLATES_ACTIVE) {
    Log log(Log::NO_FILE);
    log.StartTimer();
    template_manager_->InstantiateTemplate(cm->job_graph_name(),
        cm->inner_job_ids(),
        cm->outer_job_ids(),
        cm->parameters(),
        cm->parent_job_id().elem());
    log.StopTimer();
    std::cout << "OMID SPAWN TEMPLATE. " << log.timer() << std::endl;
  } else if (NIMBUS_NEW_TEMPLATES_ACTIVE) {
    ComplexJobEntry* complex_job;
    if (!template_manager_->GetComplexJobEntryForTemplate(complex_job,
                                                          cm->job_graph_name(),
                                                          cm->parent_job_id().elem(),
                                                          cm->inner_job_ids(),
                                                          cm->outer_job_ids(),
                                                          cm->parameters())) {
      assert(false);
    }
    if (!job_manager_->AddComplexJobEntry(complex_job)) {
      assert(false);
    }
  }
}

void Scheduler::ProcessStartTemplateCommand(StartTemplateCommand* cm) {
  // TODO(omidm): Implement!
  if (NIMBUS_TEMPLATES_ACTIVE || NIMBUS_NEW_TEMPLATES_ACTIVE) {
    template_manager_->DetectNewTemplate(cm->job_graph_name());
    template_spawner_map_[cm->parent_job_id().elem()] = cm->job_graph_name();
    std::cout << "OMID START TEMPLATE\n.";
  }
}

void Scheduler::ProcessEndTemplateCommand(EndTemplateCommand* cm) {
  // TODO(omidm): Implement!
  if (NIMBUS_TEMPLATES_ACTIVE || NIMBUS_NEW_TEMPLATES_ACTIVE) {
    template_manager_->FinalizeNewTemplate(cm->job_graph_name());
    DefinedTemplateCommand command(cm->job_graph_name());
    server_->BroadcastCommand(&command);
    std::cout << "OMID END TEMPLATE\n.";
  }
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
  job_manager_->AddMainJobEntry(j[0]);
}

size_t Scheduler::RemoveObsoleteJobEntries() {
  if (cleaner_thread_active_) {
    return 0;
  }

  return job_manager_->RemoveObsoleteJobEntries(max_job_to_remove_);
}

size_t Scheduler::AssignReadyJobs() {
  return load_balancer_->AssignReadyJobs();
}

void Scheduler::RegisterInitialWorkers() {
  while (registered_worker_num_ < min_worker_to_join_) {
    dbg(DBG_SCHED, "%d worker(s) are registered, waiting for %d more workers to join ...\n",
        registered_worker_num_, min_worker_to_join_ - registered_worker_num_);
    ProcessQueuedSchedulerCommands();
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
      // std::string ip("you-know");
      std::string ip =
        (*iter)->connection()->socket()->remote_endpoint().address().to_string();
      ID<port_t> port(0);
      {
        HandshakeCommand cm(worker_id, ip, port, Log::GetRawTime());
        dbg(DBG_SCHED, "Sending command: %s.\n", cm.ToString().c_str());
        server_->SendCommand(*iter, &cm);
      }
      for (int i = 0; i < (HANDSHAKE_COMMAND_NUM - 1); ++i) {
        sleep(1);
        HandshakeCommand cm(worker_id, ip, port, Log::GetRawTime());
        dbg(DBG_SCHED, "Sending command: %s.\n", cm.ToString().c_str());
        server_->SendCommand(*iter, &cm);
      }
    }
  }
  return registered_num;
}

void Scheduler::CreateIDMaker() {
  id_maker_ = new IDMaker();
}

void Scheduler::CreateAfterMap() {
  after_map_ = new AfterMap();
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

void Scheduler::CreateTemplateManager() {
  template_manager_ = new TemplateManager();
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
  job_manager_->set_after_map(after_map_);
  job_manager_->set_ldo_map_p(data_manager_->ldo_map_p());
}

void Scheduler::SetupTemplateManager() {
  template_manager_->set_job_manager(job_manager_);
  template_manager_->set_id_maker(id_maker_);
}

void Scheduler::SetupLoadBalancer() {
  load_balancer_->set_job_manager(job_manager_);
  load_balancer_->set_data_manager(data_manager_);
  load_balancer_->set_job_assigner(job_assigner_);
  load_balancer_->set_max_job_to_assign(max_job_to_assign_);
  load_balancer_thread_ = new boost::thread(boost::bind(&LoadBalancer::Run, load_balancer_));
}

void Scheduler::SetupJobAssigner() {
  job_assigner_->set_id_maker(id_maker_);
  job_assigner_->set_server(server_);
  job_assigner_->set_job_manager(job_manager_);
  job_assigner_->set_data_manager(data_manager_);
  job_assigner_->set_load_balancer(load_balancer_);
  job_assigner_->set_thread_num(job_assigner_thread_num_);
  job_assigner_->Run();
}

void Scheduler::SetupJobDoneBouncer() {
  job_done_bouncer_thread_ = new boost::thread(
      boost::bind(&Scheduler::JobDoneBouncerThread, this));
}

void Scheduler::SetupCleaner() {
  cleaner_thread_active_ = true;
  cleaner_thread_ = new boost::thread(
      boost::bind(&Scheduler::CleanerThread, this));
}

void Scheduler::SetupUserInterface() {
  LoadUserCommands();
  user_interface_thread_ = new boost::thread(
      boost::bind(&Scheduler::GetUserCommandThread, this));
}

void Scheduler::LoadWorkerCommands() {
  worker_command_table_[SchedulerCommand::SPAWN_COMPUTE]      = new SpawnComputeJobCommand();
  worker_command_table_[SchedulerCommand::SPAWN_COPY]         = new SpawnCopyJobCommand();
  worker_command_table_[SchedulerCommand::DEFINE_DATA]        = new DefineDataCommand();
  worker_command_table_[SchedulerCommand::HANDSHAKE]          = new HandshakeCommand();
  worker_command_table_[SchedulerCommand::JOB_DONE]           = new JobDoneCommand();
  worker_command_table_[SchedulerCommand::DEFINE_PARTITION]   = new DefinePartitionCommand();
  worker_command_table_[SchedulerCommand::TERMINATE]          = new TerminateCommand();
  worker_command_table_[SchedulerCommand::PROFILE]            = new ProfileCommand();
  worker_command_table_[SchedulerCommand::SPAWN_TEMPLATE]     = new SpawnTemplateCommand();
  worker_command_table_[SchedulerCommand::START_TEMPLATE]     = new StartTemplateCommand();
  worker_command_table_[SchedulerCommand::END_TEMPLATE]       = new EndTemplateCommand();
  worker_command_table_[SchedulerCommand::SAVE_DATA_JOB_DONE] = new SaveDataJobDoneCommand();
  worker_command_table_[SchedulerCommand::WORKER_DOWN]        = new WorkerDownCommand();
  worker_command_table_[SchedulerCommand::PREPARE_REWIND]     = new PrepareRewindCommand();
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

void Scheduler::CleanerThread() {
  while (true) {
    // TODO(omid): remove the busy loop!

    job_manager_->RemoveObsoleteJobEntries(max_job_to_remove_);
  }
}

void Scheduler::JobDoneBouncerThread() {
  while (true) {
    // TODO(omid): remove the busy loop!

    // Deal with new job dones.
    JobDoneCommandList storage;
    server_->ReceiveJobDoneCommands(&storage, max_job_done_command_process_num_);

    JobDoneCommandList::iterator iter = storage.begin();
    for (; iter != storage.end(); ++iter) {
      JobDoneCommand* comm = *iter;
      dbg(DBG_SCHED, "Bouncing job done command: %s.\n", comm->ToString().c_str());
      job_id_t job_id = comm->job_id().elem();

      after_map_->NotifyJobDone(job_id);

      if (!id_maker_->SchedulerProducedJobID(job_id)) {
        std::list<SchedulerWorker*> waiting_list;
        after_map_->GetWorkersWaitingOnJob(job_id, &waiting_list);

        comm->set_final(false);
        std::list<SchedulerWorker*>::iterator iter = waiting_list.begin();
        for (; iter != waiting_list.end(); ++iter) {
          server_->SendCommand(*iter, comm);
        }
      }
      delete comm;
    }

    // Deal with missed job dones.
    AfterMap::Map *late_map;
    if (after_map_->PullLateMap(late_map)) {
      AfterMap::Iter iter = late_map->begin();
      for (; iter != late_map->end(); ++iter) {
        if (!id_maker_->SchedulerProducedJobID(iter->first)) {
          ID<job_id_t> job_id(iter->first);
          JobDoneCommand comm(job_id);
          comm.set_final(false);

          AfterMap::Pool *pool = iter->second;
          AfterMap::Pool::iterator it = pool->begin();
          for (; it != pool->end(); ++it) {
            server_->SendCommand(*it, &comm);
          }
        }
      }

      AfterMap::DestroyMap(late_map);
    }
  }
}

void Scheduler::GetUserCommandThread() {
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
