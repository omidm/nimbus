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

#include "src/scheduler/scheduler.h"
#include "src/scheduler/dynamic_load_balancer.h"

#define HANDSHAKE_COMMAND_NUM 1

namespace nimbus {

#define DEFAULT_CLEANER_THREAD_ACTIVE true
#define DEFAULT_BOUNCER_THREAD_ACTIVE true

#define DEFAULT_CONTROLLER_TEMPLATE_ACTIVE true
#define DEFAULT_COMPLEX_MEMOIZATION_ACTIVE true
#define DEFAULT_BINDING_MEMOIZATION_ACTIVE true
#define DEFAULT_WORKER_TEMPLATE_ACTIVE true
#define DEFAULT_MEGA_RCR_JOB_ACTIVE true
#define DEFAULT_CASCADED_BINDING_ACTIVE true

#define DEFAULT_LOAD_BALANCING_ACTIVE false
#define DEFAULT_LOAD_BALANCING_PERIOD 30
#define DEFAULT_FAULT_TOLERANCE_ACTIVE false
#define DEFAULT_CHECKPOINT_CREATION_PERIOD 30


#define DEFAULT_ASSIGN_BATCH_SIZE 1
#define DEFAULT_REMOVE_BATCH_SIZE 1
#define DEFAULT_INIT_WORKER_NUM 2
#define DEFAULT_JOB_ASSIGNER_THREAD_NUM 0
#define DEFAULT_COMMAND_BATCH_SIZE 10000
#define DEFAULT_JOB_DONE_BATCH_SIZE 10000


Scheduler::Scheduler(port_t port) {
  server_ = NULL;
  id_maker_ = NULL;
  job_manager_ = NULL;
  template_manager_ = NULL;
  data_manager_ = NULL;
  job_assigner_ = NULL;
  load_balancer_ = NULL;

  listening_port_ = port;
  query_stat_id_ = 0;
  responded_worker_num_ = 0;
  registered_worker_num_ = 0;
  terminate_application_flag_ = false;
  last_query_stat_time_ = (int64_t)(Log::GetRawTime());

  split_ = std::vector<size_t>(3);
  sub_split_ = std::vector<size_t>(3);
  global_region_ = GeometricRegion();

  cleaner_thread_active_ = DEFAULT_CLEANER_THREAD_ACTIVE;
  bouncer_thread_active_ = DEFAULT_BOUNCER_THREAD_ACTIVE;

  controller_template_active_ = DEFAULT_CONTROLLER_TEMPLATE_ACTIVE;
  complex_memoization_active_ = DEFAULT_COMPLEX_MEMOIZATION_ACTIVE;
  binding_memoization_active_ = DEFAULT_BINDING_MEMOIZATION_ACTIVE;
  cascaded_binding_active_ = DEFAULT_CASCADED_BINDING_ACTIVE;
  worker_template_active_ = DEFAULT_WORKER_TEMPLATE_ACTIVE;
  mega_rcr_job_active_ = DEFAULT_MEGA_RCR_JOB_ACTIVE;

  load_balancing_active_ = DEFAULT_LOAD_BALANCING_ACTIVE;
  load_balancing_period_ = DEFAULT_LOAD_BALANCING_PERIOD;
  fault_tolerance_active_ = DEFAULT_FAULT_TOLERANCE_ACTIVE;
  checkpoint_creation_period_ = DEFAULT_CHECKPOINT_CREATION_PERIOD;

  init_worker_num_ = DEFAULT_INIT_WORKER_NUM;
  assign_batch_size_ = DEFAULT_ASSIGN_BATCH_SIZE;
  remove_batch_size_ = DEFAULT_REMOVE_BATCH_SIZE;
  job_assigner_thread_num_ = DEFAULT_JOB_ASSIGNER_THREAD_NUM;
  command_batch_size_ = DEFAULT_COMMAND_BATCH_SIZE;
  job_done_batch_size_ = DEFAULT_JOB_DONE_BATCH_SIZE;

  log_.set_file_name("log_scheduler");
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

void Scheduler::set_cleaner_thread_active(bool flag) {
  cleaner_thread_active_ = flag;
}

void Scheduler::set_bouncer_thread_active(bool flag) {
  bouncer_thread_active_ = flag;
}

void Scheduler::set_controller_template_active(bool flag) {
  controller_template_active_ = flag;
}

void Scheduler::set_complex_memoization_active(bool flag) {
  complex_memoization_active_ = flag;
}

void Scheduler::set_binding_memoization_active(bool flag) {
  binding_memoization_active_ = flag;
}

void Scheduler::set_cascaded_binding_active(bool flag) {
  cascaded_binding_active_ = flag;
}

void Scheduler::set_worker_template_active(bool flag) {
  worker_template_active_ = flag;
}

void Scheduler::set_mega_rcr_job_active(bool flag) {
  mega_rcr_job_active_ = flag;
}

void Scheduler::set_load_balancing_active(bool flag) {
  load_balancing_active_ = flag;
}

void Scheduler::set_fault_tolerance_active(bool flag) {
  fault_tolerance_active_ = flag;
}

void Scheduler::set_load_balancing_period(int64_t period) {
  load_balancing_period_ = period;
}

void Scheduler::set_checkpoint_creation_period(int64_t period) {
  checkpoint_creation_period_ = period;
}

void Scheduler::set_remove_batch_size(size_t num) {
  remove_batch_size_ = num;
}

void Scheduler::set_assign_batch_size(size_t num) {
  assign_batch_size_ = num;
}

void Scheduler::set_init_worker_num(size_t num) {
  init_worker_num_ = num;
}

void Scheduler::set_job_assigner_thread_num(size_t num) {
  job_assigner_thread_num_ = num;
}

void Scheduler::set_command_batch_size(size_t num) {
  command_batch_size_ = num;
}

void Scheduler::set_job_done_batch_size(size_t num) {
  job_done_batch_size_ = num;
}

void Scheduler::set_split_dimensions(const std::vector<size_t>& split) {
  split_ = split;
}

void Scheduler::set_sub_split_dimensions(const std::vector<size_t>& sub_split) {
  sub_split_ = sub_split;
}

void Scheduler::set_global_region(const GeometricRegion& region) {
  global_region_ = region;
}

void Scheduler::SchedulerCoreProcessor() {
  // Worker registration phase before starting the main job.
  RegisterInitialWorkers();

  // Adding main job to the job manager.
  AddMainJob();

  // reinit the last_query_stat_time_
  last_query_stat_time_ = (int64_t)(Log::GetRawTime());

  log_.log_StartTimer();
  processed_command_num_ = 0;
  // Main Loop of the scheduler.
  while (true) {
    RegisterPendingWorkers();

    QueryWorkerStats();

    log_process_.log_StartTimer();
    processed_command_num_ += ProcessQueuedSchedulerCommands();
    log_process_.log_StopTimer();

    bool overhead = (job_manager_->NumJobsReadyToAssign() > 0);
    if (overhead) {
      log_overhead_.log_AddToTimer(log_process_.timer());
    }

    if (overhead) {
      log_overhead_.ResumeTimer();
    }
    AssignReadyJobs();
    log_overhead_.StopTimer();

    RemoveObsoleteJobEntries();

    TerminationProcedure();

    if ((log_.timer() >= 5.0) || (processed_command_num_ > 10000)) {
      log_.log_StopTimer();
      char buff[LOG_MAX_BUFF_SIZE];
      snprintf(buff, sizeof(buff), "%10.9lf l: %2.2lf c: %5.0lu.",
          Log::GetRawTime(), log_.timer(), processed_command_num_);
      log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);
      log_.log_StartTimer();
      processed_command_num_ = 0;
    }
  }
}

size_t Scheduler::ProcessQueuedSchedulerCommands() {
  size_t count = 0;
  SchedulerCommandList storage;
  if (server_->ReceiveCommands(&storage, command_batch_size_)) {
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
    case SchedulerCommand::MEGA_JOB_DONE:
      ProcessMegaJobDoneCommand(reinterpret_cast<MegaJobDoneCommand*>(cm));
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
    case SchedulerCommand::RESPOND_STAT:
      ProcessRespondStatCommand(reinterpret_cast<RespondStatCommand*>(cm));
      break;
    default:
      dbg(DBG_ERROR, "ERROR: %s have not been implemented in ProcessSchedulerCommand yet.\n",
          cm->ToNetworkData().c_str());
      assert(false);
  }
}

void Scheduler::ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm) {
  JobEntry * job =
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
    TemplateJobEntry* template_job =
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
    job->set_memoize(true);
    job->set_template_job(template_job);
  }
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
        snprintf(buff, sizeof(buff), "%10.9lf Registered new worker, id: %u IP: %s port: %u time: %2.2lf.", // NOLINT
            Log::GetRawTime(), (*iter)->worker_id(),
            (*iter)->ip().c_str(), (*iter)->port(), log_.GetTime());
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
  --registered_worker_num_;

  data_manager_->RemoveAllInstanceByWorker(worker_id);
  data_manager_->ResetAllInstances();

  checkpoint_id_t checkpoint_id = 0;
  job_manager_->RewindFromLastCheckpoint(&checkpoint_id);

  job_assigner_->Reinitialize(checkpoint_id);

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
    if (server_->ReceiveCommands(&storage, command_batch_size_)) {
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
  assert(false);
}

void Scheduler::ProcessRespondStatCommand(RespondStatCommand* cm) {
  responded_worker_num_++;
  load_balancer_->AddWorkerStat(cm->query_id(),
                                cm->worker_id(),
                                cm->run_time(),
                                cm->block_time(),
                                cm->idle_time());

  std::cout << "\n** Stat Worker: " << cm->worker_id() << std::endl;
  std::cout << "Stat Run:   " << cm->run_time() / (double)1000000000 << std::endl; // NOLINT
  std::cout << "Stat Block: " << cm->block_time() / (double)1000000000 << std::endl; // NOLINT
  std::cout << "Stat Idle:  " << cm->idle_time() / (double)1000000000 << std::endl; // NOLINT
  std::cout << "*******************************" << cm->worker_id() << std::endl; // NOLINT
}

void Scheduler::ProcessMegaJobDoneCommand(MegaJobDoneCommand* cm) {
  std::vector<job_id_t>::const_iterator iter = cm->job_ids_p()->begin();
  for (; iter != cm->job_ids_p()->end(); ++iter) {
    job_manager_->NotifyJobDone(*iter);
  }
}

void Scheduler::ProcessJobDoneCommand(JobDoneCommand* cm) {
  job_id_t job_id = cm->job_id().elem();

  // Since we are flooding now, make sure that just compute jobs are dealt
  // here. workers should filter copy job dones. -omidm
  assert(!IDMaker::SchedulerProducedJobID(job_id) ||
         job_id == NIMBUS_KERNEL_JOB_ID + 1);

  if (!id_maker_->SchedulerProducedJobID(job_id) ||
      job_id == NIMBUS_KERNEL_JOB_ID + 1) {
    // No need for flooding with binding templates anymore. After map still
    // does not work with binding template, but ther is no explicit job done
    // needed with binding templates, anymore. -omidm
    // if (binding_memoization_active_) {
    //   cm->set_final(true);
    //   SchedulerWorkerList::iterator iter = server_->workers()->begin();
    //   for (; iter != server_->workers()->end(); ++iter) {
    //     server_->SendCommand(*iter, cm);
    //   }
    // Flood job main done, so that workers starts counting idle time. -omidm
    if (job_id == NIMBUS_KERNEL_JOB_ID + 1) {
      cm->set_final(true);
      SchedulerWorkerList::iterator iter = server_->workers()->begin();
      for (; iter != server_->workers()->end(); ++iter) {
        server_->SendCommand(*iter, cm);
      }
    } else {
      std::list<SchedulerWorker*> waiting_list;
      after_map_->GetWorkersWaitingOnJob(job_id, &waiting_list);

      cm->set_final(true);
      std::list<SchedulerWorker*>::iterator iter = waiting_list.begin();
      for (; iter != waiting_list.end(); ++iter) {
        server_->SendCommand(*iter, cm);
      }
    }
  }

  // AfterMap has internal locking.
  after_map_->RemoveJobRecords(job_id);

  // TODO(omidm): load balancer does not get notified for now!
  job_manager_->NotifyJobDone(job_id);
//  JobEntry *job;
//  if (job_manager_->GetJobEntry(job_id, job)) {
//    job->set_done(true);
//    load_balancer_->NotifyJobDone(job);
//    job_manager_->NotifyJobDone(job);
//  }
}

void Scheduler::ProcessTerminateCommand(TerminateCommand* cm) {
  terminate_application_flag_ = true;
  terminate_application_status_ = cm->exit_status().elem();
}

void Scheduler::ProcessSpawnTemplateCommand(SpawnTemplateCommand* cm) {
  if (complex_memoization_active_) {
    Log log(Log::NO_FILE);
    log.StartTimer();
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
    log.StopTimer();
    std::cout << "TEMPLATE: COMPLEX SPAWN: " << cm->job_graph_name() << " " << log.timer() << std::endl; // NOLINT
  } else {
    assert(controller_template_active_);
    Log log(Log::NO_FILE);
    log.StartTimer();
    template_manager_->InstantiateTemplate(cm->job_graph_name(),
                                           cm->inner_job_ids(),
                                           cm->outer_job_ids(),
                                           cm->parameters(),
                                           cm->parent_job_id().elem());
    log.StopTimer();
    std::cout << "TEMPLATE: SPAWN: " << cm->job_graph_name() << " " << log.timer() << std::endl;
  }
}

void Scheduler::ProcessStartTemplateCommand(StartTemplateCommand* cm) {
  if (controller_template_active_) {
    std::string template_name = cm->job_graph_name();
    job_id_t job_id = cm->parent_job_id().elem();

    if (!template_manager_->DetectNewTemplate(template_name)) {
      dbg(DBG_ERROR, "ERROR: could not detect new template %s.\n", template_name.c_str());
      assert(false);
      return;
    }

    template_spawner_map_[job_id] = template_name;

    // Memoize the parent version map as the base version map of the template.
    boost::shared_ptr<VersionMap> vmap_base;
    if (!job_manager_->GetBaseVersionMapFromJob(job_id, vmap_base)) {
      dbg(DBG_ERROR, "ERROR: could not get base version map from parent of template with job id %lu.\n", job_id); // NOLINT
      assert(false);
    }
    template_manager_->SetBaseVersionMapForTemplate(template_name,
                                                    vmap_base);
    std::cout << "TEMPLATE: START " << template_name << std::endl;
  }
}

void Scheduler::ProcessEndTemplateCommand(EndTemplateCommand* cm) {
  if (controller_template_active_) {
    template_manager_->FinalizeNewTemplate(cm->job_graph_name());
    DefinedTemplateCommand command(cm->job_graph_name());
    server_->BroadcastCommand(&command);
    std::cout << "TEMPLATE: END " << cm->job_graph_name() << std::endl;
  }
}

void Scheduler::TerminationProcedure() {
  if (terminate_application_flag_) {
    if (job_manager_->AllJobsAreDone()) {
      SchedulerCommand* command =
        new TerminateCommand(ID<exit_status_t>(terminate_application_status_));
      server_->BroadcastCommand(command);
      delete command;
      char buff[LOG_MAX_BUFF_SIZE];
      snprintf(buff, sizeof(buff), "%10.9lf Simulation Terminated.", Log::GetRawTime());
      log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);
      snprintf(buff, sizeof(buff), "Controller overhead: %2.5lf seconds.", log_overhead_.timer());
      log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);
      snprintf(buff, sizeof(buff), "Controller bytes sent: %4.2lf MB.",
          static_cast<double>(server_->total_bytes_sent()) / 1e6); // NOLINT
      log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);
      snprintf(buff, sizeof(buff), "Controller bytes received: %4.2lf MB.",
          static_cast<double>(server_->total_bytes_received()) / 1e6);
      log_.log_WriteToOutputStream(std::string(buff), LOG_INFO);
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

  return job_manager_->RemoveObsoleteJobEntries(remove_batch_size_);
}

size_t Scheduler::AssignReadyJobs() {
  return load_balancer_->AssignReadyJobs();
}

void Scheduler::QueryWorkerStats() {
  if (!load_balancing_active_) {
    return;
  }

  if (!load_balancer_->safe_to_load_balance()) {
    return;
  }

  int64_t time = (int64_t)(Log::GetRawTime());

  if ((time - last_query_stat_time_) >= load_balancing_period_) {
    std::cout << "\n\n***** Begin Stat Query *****\n";
    RequestStatCommand command(query_stat_id_);
    server_->BroadcastCommand(&command);

    while (responded_worker_num_ < registered_worker_num_) {
      ProcessQueuedSchedulerCommands();
    }

    load_balancer_->BalanceLoad(query_stat_id_);

    ++query_stat_id_;
    responded_worker_num_ = 0;
    last_query_stat_time_ = (int64_t)(Log::GetRawTime());
  }
}



void Scheduler::RegisterInitialWorkers() {
  while (registered_worker_num_ < init_worker_num_) {
    dbg(DBG_SCHED, "%d worker(s) are registered, waiting for %d more workers to join ...\n",
        registered_worker_num_, init_worker_num_ - registered_worker_num_);
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
  load_balancer_ = new DynamicLoadBalancer();
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
  server_->set_bouncer_thread_active(bouncer_thread_active_);
  scheduler_server_thread_ = new boost::thread(boost::bind(&SchedulerServer::Run, server_));
}

void Scheduler::SetupDataManager() {
}

void Scheduler::SetupJobManager() {
  job_manager_->set_after_map(after_map_);
  job_manager_->set_ldo_map_p(data_manager_->ldo_map_p());
  job_manager_->set_binding_memoization_active(binding_memoization_active_);
  job_manager_->set_cascaded_binding_active(cascaded_binding_active_);
  job_manager_->set_fault_tolerance_active(fault_tolerance_active_);
  job_manager_->set_checkpoint_creation_period(checkpoint_creation_period_);
}

void Scheduler::SetupTemplateManager() {
  template_manager_->set_job_manager(job_manager_);
  template_manager_->set_id_maker(id_maker_);
  template_manager_->set_worker_template_active(worker_template_active_);
  template_manager_->set_mega_rcr_job_active(mega_rcr_job_active_);
}

void Scheduler::SetupLoadBalancer() {
  load_balancer_->set_job_manager(job_manager_);
  load_balancer_->set_data_manager(data_manager_);
  load_balancer_->set_job_assigner(job_assigner_);
  load_balancer_->set_assign_batch_size(assign_batch_size_);
  load_balancer_->set_split_dimensions(split_);
  load_balancer_->set_sub_split_dimensions(sub_split_);
  load_balancer_->set_global_region(global_region_);
  load_balancer_thread_ = new boost::thread(boost::bind(&LoadBalancer::Run, load_balancer_));
}

void Scheduler::SetupJobAssigner() {
  job_assigner_->set_id_maker(id_maker_);
  job_assigner_->set_server(server_);
  job_assigner_->set_job_manager(job_manager_);
  job_assigner_->set_data_manager(data_manager_);
  job_assigner_->set_load_balancer(load_balancer_);
  job_assigner_->set_thread_num(job_assigner_thread_num_);
  job_assigner_->set_fault_tolerance_active(fault_tolerance_active_);
  job_assigner_->set_scheduler(this);
  job_assigner_->set_ldo_map_p(data_manager_->ldo_map_p());
  job_assigner_->set_log_overhead(&log_overhead_);
  job_assigner_->Run();
}

void Scheduler::SetupJobDoneBouncer() {
  if (bouncer_thread_active_) {
    job_done_bouncer_thread_ = new boost::thread(
        boost::bind(&Scheduler::JobDoneBouncerThread, this));
  }
}

void Scheduler::SetupCleaner() {
  if (cleaner_thread_active_) {
    cleaner_thread_ = new boost::thread(
        boost::bind(&Scheduler::CleanerThread, this));
  }
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
  worker_command_table_[SchedulerCommand::MEGA_JOB_DONE]      = new MegaJobDoneCommand();
  worker_command_table_[SchedulerCommand::DEFINE_PARTITION]   = new DefinePartitionCommand();
  worker_command_table_[SchedulerCommand::TERMINATE]          = new TerminateCommand();
  worker_command_table_[SchedulerCommand::PROFILE]            = new ProfileCommand();
  worker_command_table_[SchedulerCommand::SPAWN_TEMPLATE]     = new SpawnTemplateCommand();
  worker_command_table_[SchedulerCommand::START_TEMPLATE]     = new StartTemplateCommand();
  worker_command_table_[SchedulerCommand::END_TEMPLATE]       = new EndTemplateCommand();
  worker_command_table_[SchedulerCommand::SAVE_DATA_JOB_DONE] = new SaveDataJobDoneCommand();
  worker_command_table_[SchedulerCommand::WORKER_DOWN]        = new WorkerDownCommand();
  worker_command_table_[SchedulerCommand::PREPARE_REWIND]     = new PrepareRewindCommand();
  worker_command_table_[SchedulerCommand::RESPOND_STAT]       = new RespondStatCommand();
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

    job_manager_->RemoveObsoleteJobEntries(remove_batch_size_);
  }
}

void Scheduler::JobDoneBouncerThread() {
  while (true) {
    // TODO(omid): remove the busy loop!

    // Deal with new job dones.
    JobDoneCommandList storage;
    server_->ReceiveJobDoneCommands(&storage, job_done_batch_size_);

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

void Scheduler::PrintStats() {
  static FILE* file = fopen("controller_stats.txt", "w");
  static uint64_t l_sent = 0, l_received = 0;
  static double l_overhead = 0;
  static size_t counter = 1;

  uint64_t c_sent = server_->total_bytes_sent();
  uint64_t c_received = server_->total_bytes_received();
  double c_overhead = log_overhead_.timer();

  fprintf(file, "%10.9lf %3.1lu sent(MB): %.4f received(MB): %.4f overhead(s): %.4f\n",
      Log::GetRawTime(),
      counter,
      static_cast<double>(c_sent - l_sent) / 1e6,
      static_cast<double>(c_received - l_received) / 1e6,
      c_overhead - l_overhead);
  fflush(file);

  l_sent = c_sent;
  l_received = c_received;
  l_overhead = c_overhead;
  counter++;
}


}  // namespace nimbus
