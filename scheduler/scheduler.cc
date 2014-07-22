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
#define MAX_JOB_TO_ASSIGN 1

Scheduler::Scheduler(port_t p)
: listening_port_(p) {
  appId_ = 0;
  registered_worker_num_ = 0;
  data_manager_ = NULL;
  job_manager_ = NULL;
  min_worker_to_join_ = DEFAULT_MIN_WORKER_TO_JOIN;
  terminate_application_flag_ = false;
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

  SetupWorkerInterface();

  // SetupUserInterface();

  // First data manager should be instantiated then job manager, because of
  // ldo_map pointer that job manager needs to get from data manager.
  // Also load balancer is set up after both of them since it needs a pointer
  // to job manager and data manager.

  SetupDataManager();
  SetupJobManager();
  SetupLoadBalancer();

  id_maker_.Initialize(0);

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
    log_loop_.StartTimer();
    log_assign_.ResetTimer();
    log_job_manager_.ResetTimer();
    log_data_manager_.ResetTimer();
    log_version_manager_.ResetTimer();
    log_load_balancer_.ResetTimer();

    RegisterPendingWorkers();
    ProcessQueuedSchedulerCommands((size_t)MAX_BATCH_COMMAND_NUM);
    AssignReadyJobs();
    RemoveObsoleteJobEntries();
    TerminationProcedure();

    log_loop_.StopTimer();
    if (log_loop_.timer() >= .00001) {
      char buff[LOG_MAX_BUFF_SIZE];
      snprintf(buff, sizeof(buff),
          "loop: %2.5lf  assign: %2.5lf job_manager: %2.5lf data_manager: %2.5lf version_manager: %2.5lf load_balancer: %2.5lf time: %6.5lf", // NOLINT
          log_loop_.timer(),
          log_assign_.timer(),
          log_job_manager_.timer(),
          log_data_manager_.timer(),
          log_version_manager_.timer(),
          log_load_balancer_.timer(),
          log_loop_.GetTime());

      log_.WriteToOutputStream(std::string(buff), LOG_INFO);
    }
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
    case SchedulerCommand::TERMINATE:
      ProcessTerminateCommand(reinterpret_cast<TerminateCommand*>(cm));
      break;
    case SchedulerCommand::PROFILE:
      break;
    default:
      dbg(DBG_ERROR, "ERROR: %s have not been implemented in ProcessSchedulerCommand yet.\n",
          cm->toString().c_str());
  }
}

void Scheduler::ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm) {
  log_job_manager_.ResumeTimer();
  job_manager_->AddComputeJobEntry(
        cm->job_name(),
        cm->job_id().elem(),
        cm->read_set(),
        cm->write_set(),
        cm->before_set(),
        cm->after_set(),
        cm->parent_job_id().elem(),
        cm->params(),
        cm->sterile());
  log_job_manager_.StopTimer();
}

void Scheduler::ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm) {
  std::string job_name = "copyjob";
  IDSet<logical_data_id_t> read_set;
  read_set.insert(cm->from_logical_id().elem());
  IDSet<logical_data_id_t> write_set;
  write_set.insert(cm->to_logical_id().elem());

  // TODO(omid): we need to add support for copy jobs.
  log_job_manager_.ResumeTimer();
  job_manager_->AddExplicitCopyJobEntry();
  log_job_manager_.StopTimer();
}

void Scheduler::ProcessDefineDataCommand(DefineDataCommand* cm) {
  log_data_manager_.ResumeTimer();
  bool success = data_manager_->AddLogicalObject(cm->logical_data_id().elem(),
                                 cm->data_name(),
                                 cm->partition_id().elem());
  log_data_manager_.StopTimer();

  if (success) {
    log_data_manager_.ResumeTimer();
    LdoAddCommand command(data_manager_->FindLogicalObject(cm->logical_data_id().elem()));
    log_data_manager_.StopTimer();
    server_->BroadcastCommand(&command);
  }

  log_job_manager_.ResumeTimer();
  job_manager_->DefineData(cm->parent_job_id().elem(),
                          cm->logical_data_id().elem());
  log_job_manager_.StopTimer();
}

void Scheduler::ProcessDefinePartitionCommand(DefinePartitionCommand* cm) {
  GeometricRegion r = *(cm->region());
  log_data_manager_.ResumeTimer();
  data_manager_->AddPartition(cm->partition_id().elem(), r);
  log_data_manager_.StopTimer();
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
      }
      break;
    }
  }
}

void Scheduler::ProcessJobDoneCommand(JobDoneCommand* cm) {
  job_id_t job_id = cm->job_id().elem();

  log_job_manager_.ResumeTimer();
  job_manager_->NotifyJobDone(job_id);
  log_job_manager_.StopTimer();

  JobEntry *job;
  log_job_manager_.ResumeTimer();
  if (job_manager_->GetJobEntry(job_id, job)) {
    log_load_balancer_.ResumeTimer();
    load_balancer_->NotifyJobDone(job);
    log_load_balancer_.StopTimer();
  }
  log_job_manager_.StopTimer();

  std::string jname = job->job_name();
  if (jname == "loop_iteration") {
    log_.StartTimer();
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
  id_maker_.GetNewJobID(&j, 1);
  log_job_manager_.ResumeTimer();
  job_manager_->AddMainJobEntry(j[0]);
  log_job_manager_.StopTimer();
}

bool Scheduler::GetWorkerToAssignJob(JobEntry* job, SchedulerWorker*& worker) {
  // Simply just assign the job to first worker.
  dbg(DBG_WARN, "WARNING: Base scheduler only assigns jobs to first worker, override the GetWorkerToAssignJob to have more complicated logic.\n"); // NOLINT
  worker =  *(server_->workers()->begin());
  return true;
}

size_t Scheduler::RemoveObsoleteJobEntries() {
  log_job_manager_.ResumeTimer();
  size_t count = job_manager_->RemoveObsoleteJobEntries();
  log_job_manager_.StopTimer();

  return count;
}

bool Scheduler::AllocateLdoInstanceToJob(JobEntry* job,
    LogicalDataObject* ldo, PhysicalData pd) {
  assert(job->versioned());
  // IDSet<job_id_t> before_set = job->before_set();
  PhysicalData pd_new = pd;

  // data_version_t v_in, // v_out;
  // job->vtable_in()->que// ry_entry(ldo->id(), &v_in);
  // job->vtable_out()->qu// ery_entry(ldo->id(), &v_out);

  // Because of the clear_list_job_read the order of if blocks are important.
  if (job->write_set_p()->contains(ldo->id())) {
    // pd_new.set_version(job->version_table_out_query(ldo->id()));
    data_version_t v_out;
    // job->vmap_write_out()->query_entry(ldo->id(), &v_out);
    job->vmap_write()->query_entry(ldo->id(), &v_out);
    pd_new.set_version(v_out);
    pd_new.set_last_job_write(job->job_id());
    pd_new.clear_list_job_read();
    // before_set.insert(pd.list_job_read());
    job->before_set_p()->insert(pd.list_job_read());
    // become ancestor of last job write if not already ,so that read access
    // is safe for the jobs in list job read cleared by the previous write job - omidm
    // before_set.insert(pd.last_job_write());
    job->before_set_p()->insert(pd.last_job_write());
  }

  if (job->read_set_p()->contains(ldo->id())) {
    // assert(job->version_table_in_query(ldo->id()) == pd.version());
    data_version_t v_in;
    // job->vmap_read_in()->query_entry(ldo->id(), &v_in);
    job->vmap_read()->query_entry(ldo->id(), &v_in);
    assert(v_in == pd.version());
    pd_new.add_to_list_job_read(job->job_id());
    // before_set.insert(pd.last_job_write());
    job->before_set_p()->insert(pd.last_job_write());
  }

  job->set_physical_table_entry(ldo->id(), pd.id());
  // job->set_before_set(before_set);

  log_data_manager_.ResumeTimer();
  data_manager_->RemovePhysicalInstance(ldo, pd);
  data_manager_->AddPhysicalInstance(ldo, pd_new);
  log_data_manager_.StopTimer();

  return true;
}

size_t Scheduler::GetObsoleteLdoInstancesAtWorker(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalDataVector* dest) {
  size_t count = 0;
  dest->clear();
  PhysicalDataVector pv;
  log_data_manager_.ResumeTimer();
  data_manager_->InstancesByWorker(ldo, worker->worker_id(), &pv);
  log_data_manager_.StopTimer();
  PhysicalDataVector::iterator iter = pv.begin();
  for (; iter != pv.end(); ++iter) {
    JobEntryList list;
    VersionedLogicalData vld(ldo->id(), iter->version());
    log_job_manager_.ResumeTimer();
    log_version_manager_.ResumeTimer();
    if (job_manager_->GetJobsNeedDataVersion(&list, vld) == 0) {
      dest->push_back(*iter);
      ++count;
    }
    log_job_manager_.StopTimer();
    log_version_manager_.StopTimer();
  }
  return count;
}

bool Scheduler::CreateDataAtWorker(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalData* created_data) {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  std::vector<physical_data_id_t> d;
  id_maker_.GetNewPhysicalDataID(&d, 1);
  IDSet<job_id_t> before, after;

  // Update the job table.
  log_job_manager_.ResumeTimer();
  job_manager_->AddCreateDataJobEntry(j[0]);
  log_job_manager_.StopTimer();

  // Update data table.
  IDSet<job_id_t> list_job_read;
  list_job_read.insert(j[0]);  // if other job wants to write, waits for creation.
  PhysicalData p(d[0], worker->worker_id(), NIMBUS_INIT_DATA_VERSION, list_job_read, j[0]);
  log_data_manager_.ResumeTimer();
  data_manager_->AddPhysicalInstance(ldo, p);
  log_data_manager_.StopTimer();

  // send the create command to worker.
  log_job_manager_.ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.StopTimer();
  CreateDataCommand cm(ID<job_id_t>(j[0]), ldo->variable(),
      ID<logical_data_id_t>(ldo->id()),
      ID<physical_data_id_t>(d[0]), before, after);
  server_->SendCommand(worker, &cm);

  *created_data = p;

  return true;
}

bool Scheduler::RemoteCopyData(SchedulerWorker* from_worker,
    SchedulerWorker* to_worker, LogicalDataObject* ldo,
    PhysicalData* from_data, PhysicalData* to_data) {
  assert(from_worker->worker_id() == from_data->worker());
  assert(to_worker->worker_id() == to_data->worker());

  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 2);
  job_id_t receive_id = j[0];
  job_id_t send_id = j[1];
  IDSet<job_id_t> before, after;

  // Receive part

  // Update the job table.
  log_job_manager_.ResumeTimer();
  job_manager_->AddRemoteCopyReceiveJobEntry(receive_id);
  log_job_manager_.StopTimer();

  // Update data table.
  PhysicalData to_data_new = *to_data;
  to_data_new.set_version(from_data->version());
  to_data_new.set_last_job_write(receive_id);
  to_data_new.clear_list_job_read();
  log_data_manager_.ResumeTimer();
  data_manager_->RemovePhysicalInstance(ldo, *to_data);
  data_manager_->AddPhysicalInstance(ldo, to_data_new);
  log_data_manager_.StopTimer();

  // send remote copy receive job to worker.
  before.clear();
  before.insert(to_data->list_job_read());
  before.insert(to_data->last_job_write());
  log_job_manager_.ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.StopTimer();
  RemoteCopyReceiveCommand cm_r(ID<job_id_t>(receive_id),
      ID<physical_data_id_t>(to_data->id()), before, after);
  server_->SendCommand(to_worker, &cm_r);


  // Send Part.

  // Update the job table.
  log_job_manager_.ResumeTimer();
  job_manager_->AddRemoteCopySendJobEntry(send_id);
  log_job_manager_.StopTimer();

  // Update data table.
  PhysicalData from_data_new = *from_data;
  from_data_new.add_to_list_job_read(send_id);
  log_data_manager_.ResumeTimer();
  data_manager_->RemovePhysicalInstance(ldo, *from_data);
  data_manager_->AddPhysicalInstance(ldo, from_data_new);
  log_data_manager_.StopTimer();

  // send remote copy send command to worker.
  before.clear();
  before.insert(from_data->last_job_write());
  log_job_manager_.ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.StopTimer();
  RemoteCopySendCommand cm_s(ID<job_id_t>(send_id),
      ID<job_id_t>(receive_id), ID<physical_data_id_t>(from_data->id()),
      ID<worker_id_t>(to_worker->worker_id()),
      to_worker->ip(), ID<port_t>(to_worker->port()),
      before, after);
  server_->SendCommand(from_worker, &cm_s);


  *from_data = from_data_new;
  *to_data = to_data_new;

  return true;
}

bool Scheduler::LocalCopyData(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalData* from_data, PhysicalData* to_data) {
  assert(worker->worker_id() == from_data->worker());
  assert(worker->worker_id() == to_data->worker());

  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  IDSet<job_id_t> before, after;

  // Update the job table.
  log_job_manager_.ResumeTimer();
  job_manager_->AddLocalCopyJobEntry(j[0]);
  log_job_manager_.StopTimer();

  // Update data table.
  PhysicalData from_data_new = *from_data;
  from_data_new.add_to_list_job_read(j[0]);
  log_data_manager_.ResumeTimer();
  data_manager_->RemovePhysicalInstance(ldo, *from_data);
  data_manager_->AddPhysicalInstance(ldo, from_data_new);
  log_data_manager_.StopTimer();

  PhysicalData to_data_new = *to_data;
  to_data_new.set_version(from_data->version());
  to_data_new.set_last_job_write(j[0]);
  to_data_new.clear_list_job_read();
  log_data_manager_.ResumeTimer();
  data_manager_->RemovePhysicalInstance(ldo, *to_data);
  data_manager_->AddPhysicalInstance(ldo, to_data_new);
  log_data_manager_.StopTimer();

  // send local copy command to worker.
  before.insert(to_data->list_job_read());
  before.insert(to_data->last_job_write());
  before.insert(from_data->last_job_write());
  log_job_manager_.ResumeTimer();
  job_manager_->UpdateBeforeSet(&before);
  log_job_manager_.StopTimer();
  LocalCopyCommand cm_c(ID<job_id_t>(j[0]),
      ID<physical_data_id_t>(from_data->id()),
      ID<physical_data_id_t>(to_data->id()), before, after);
  server_->SendCommand(worker, &cm_c);

  *from_data = from_data_new;
  *to_data = to_data_new;

  return true;
}

bool Scheduler::GetFreeDataAtWorker(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalData* free_data) {
  PhysicalDataVector obsolete_instances;
  if (GetObsoleteLdoInstancesAtWorker(worker, ldo, &obsolete_instances) > 0) {
    *free_data = obsolete_instances[0];
    return true;
  }

  // Since there are no obsoletes, go ahead and create a new one.
  return CreateDataAtWorker(worker, ldo, free_data);
}


bool Scheduler::PrepareDataForJobAtWorker(JobEntry* job,
    SchedulerWorker* worker, logical_data_id_t l_id) {
  bool reading = job->read_set_p()->contains(l_id);
  bool writing = job->write_set_p()->contains(l_id);
  assert(reading || writing);

  log_data_manager_.ResumeTimer();
  LogicalDataObject* ldo =
    const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(l_id));
  log_data_manager_.StopTimer();


  // data_version_t version = job->version_table_in_query(l_id);

  // data_version_t version;
  // if (!job->vtable_in()->query_entry(l_id, &version)) {
  //   dbg(DBG_ERROR, "ERROR: logical id %lu is not versioned in the context of %s.\n",
  //       l_id, job->job_name().c_str());
  // }

//  data_version_t version;
//  if (reading) {
//    if (!job->vmap_read_in()->query_entry(l_id, &version)) {
//      dbg(DBG_ERROR, "ERROR: logical id %lu is not versioned in the read context of %s.\n",
//          l_id, job->job_name().c_str());
//      exit(-1);
//    }
//  }
//
//  // Just for checking
//  data_version_t unused_version;
//  if (writing) {
//    if (!job->vmap_write_out()->query_entry(l_id, &unused_version)) {
//      dbg(DBG_ERROR, "ERROR: logical id %lu is not versioned in the write context of %s.\n",
//          l_id, job->job_name().c_str());
//      exit(-1);
//    }
//  }



  data_version_t version;
  if (reading) {
    if (!job->vmap_read()->query_entry(l_id, &version)) {
      dbg(DBG_ERROR, "ERROR: logical id %lu is not versioned in the read context of %s.\n",
          l_id, job->job_name().c_str());
      exit(-1);
    }
  }

  // Just for checking
  data_version_t unused_version;
  if (writing) {
    if (!job->vmap_write()->query_entry(l_id, &unused_version)) {
      dbg(DBG_ERROR, "ERROR: logical id %lu is not versioned in the write context of %s.\n",
          l_id, job->job_name().c_str());
      exit(-1);
    }
  }



  // Checking correctness of the new versioning system
  // data_version_t version_c;
  // if (reading) {
  //   assert(job->vmap_read_in()->query_entry(l_id, &version_c));
  //   assert(version_c == version);
  // }
  // if (writing) {
  //   assert(job->vmap_write_out()->query_entry(l_id, &version_c));
  //   assert(version_c == (version + 1));
  // }


  if (!reading) {
    PhysicalData target_instance;
    GetFreeDataAtWorker(worker, ldo, &target_instance);

    log_job_manager_.ResumeTimer();
    if (job_manager_->CausingUnwantedSerialization(job, l_id, target_instance)) {
      dbg(DBG_SCHED, "Causing unwanted serialization for data %lu.\n", l_id);
    }
    log_job_manager_.StopTimer();

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }

  PhysicalDataVector instances_at_worker;
  PhysicalDataVector instances_in_system;
  log_data_manager_.ResumeTimer();
  data_manager_->InstancesByWorkerAndVersion(
      ldo, worker->worker_id(), version, &instances_at_worker);
  data_manager_->InstancesByVersion(ldo, version, &instances_in_system);
  log_data_manager_.StopTimer();

  JobEntryList list;
  VersionedLogicalData vld(l_id, version);
  log_job_manager_.ResumeTimer();
  log_version_manager_.ResumeTimer();
  job_manager_->GetJobsNeedDataVersion(&list, vld);
  log_job_manager_.StopTimer();
  log_version_manager_.StopTimer();
  assert(list.size() >= 1);
  bool writing_needed_version = (list.size() > 1) && writing;


  if (instances_at_worker.size() > 1) {
    PhysicalData target_instance;

    bool found = false;
    PhysicalDataVector::iterator iter;
    for (iter = instances_at_worker.begin(); iter != instances_at_worker.end(); iter++) {
      log_job_manager_.ResumeTimer();
      if (!job_manager_->CausingUnwantedSerialization(job, l_id, *iter)) {
        target_instance = *iter;
        found = true;
        break;
      }
      log_job_manager_.StopTimer();
    }

    if (!found) {
      dbg(DBG_SCHED, "Avoiding unwanted serialization for data %lu (1).\n", l_id);
      GetFreeDataAtWorker(worker, ldo, &target_instance);
      LocalCopyData(worker, ldo, &instances_at_worker[0], &target_instance);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }


  if ((instances_at_worker.size() == 1) && !writing_needed_version) {
    PhysicalData target_instance;

    log_job_manager_.ResumeTimer();
    if (!job_manager_->CausingUnwantedSerialization(job, l_id, instances_at_worker[0])) {
      target_instance = instances_at_worker[0];
      log_job_manager_.StopTimer();
    } else {
      log_job_manager_.StopTimer();
      dbg(DBG_SCHED, "Avoiding unwanted serialization for data %lu (2).\n", l_id);
      GetFreeDataAtWorker(worker, ldo, &target_instance);
      LocalCopyData(worker, ldo, &instances_at_worker[0], &target_instance);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }


  if ((instances_at_worker.size() == 1) && writing_needed_version) {
    PhysicalData target_instance;

    log_job_manager_.ResumeTimer();
    if (!job_manager_->CausingUnwantedSerialization(job, l_id, instances_at_worker[0])) {
      log_job_manager_.StopTimer();
      target_instance = instances_at_worker[0];
      PhysicalData copy_data;
      GetFreeDataAtWorker(worker, ldo, &copy_data);
      LocalCopyData(worker, ldo, &target_instance, &copy_data);
    } else {
      log_job_manager_.StopTimer();
      dbg(DBG_SCHED, "Avoiding unwanted serialization for data %lu (3).\n", l_id);
      GetFreeDataAtWorker(worker, ldo, &target_instance);
      LocalCopyData(worker, ldo, &instances_at_worker[0], &target_instance);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }


  if ((instances_at_worker.size() == 0) && (version == NIMBUS_INIT_DATA_VERSION)) {
    PhysicalData created_data;
    CreateDataAtWorker(worker, ldo, &created_data);

    AllocateLdoInstanceToJob(job, ldo, created_data);
    return true;
  }


  if ((instances_at_worker.size() == 0) && (instances_in_system.size() >= 1)) {
    PhysicalData from_instance = instances_in_system[0];
    worker_id_t sender_id = from_instance.worker();
    SchedulerWorker* worker_sender;
    if (!server_->GetSchedulerWorkerById(worker_sender, sender_id)) {
      dbg(DBG_ERROR, "ERROR: could not find worker with id %lu.\n", sender_id);
      exit(-1);
    }

    PhysicalData target_instance;
    GetFreeDataAtWorker(worker, ldo, &target_instance);
    RemoteCopyData(worker_sender, worker, ldo, &from_instance, &target_instance);

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }

  dbg(DBG_ERROR, "ERROR: the version (%lu) of logical data %s (%lu) needed for job %s (%lu) does not exist.\n", // NOLINT
      version, ldo->variable().c_str(), l_id, job->job_name().c_str(), job->job_id());
  assert(instances_in_system.size() >= 1);

  return false;
}

bool Scheduler::SendComputeJobToWorker(SchedulerWorker* worker, JobEntry* job) {
  if (job->job_type() == JOB_COMP) {
    ID<job_id_t> id(job->job_id());
    IDSet<physical_data_id_t> read_set, write_set;
    // TODO(omidm): check the return value of the following methods.
    job->GetPhysicalReadSet(&read_set);
    job->GetPhysicalWriteSet(&write_set);
    ComputeJobCommand cm(job->job_name(), id,
        read_set, write_set, job->before_set(), job->after_set(), job->params(), job->sterile());
    dbg(DBG_SCHED, "Sending compute job %lu to worker %lu.\n", job->job_id(), worker->worker_id());
    server_->SendCommand(worker, &cm);
    return true;
  } else {
    dbg(DBG_ERROR, "Job with id %lu is not a compute job.\n", job->job_id());
    return false;
  }
}

bool Scheduler::SendCreateJobToWorker(SchedulerWorker* worker,
    const std::string& data_name, const logical_data_id_t& logical_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    job_id_t* job_id, physical_data_id_t* physical_data_id) {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  *job_id = j[0];
  std::vector<physical_data_id_t> d;
  id_maker_.GetNewPhysicalDataID(&d, 1);
  *physical_data_id = d[0];
  CreateDataCommand cm(ID<job_id_t>(j[0]), data_name,
      ID<logical_data_id_t>(logical_data_id), ID<physical_data_id_t>(d[0]), before, after);
  dbg(DBG_SCHED, "Sending create job %lu to worker %lu.\n", j[0], worker->worker_id());
  server_->SendCommand(worker, &cm);
  log_job_manager_.ResumeTimer();
  job_manager_->AddCreateDataJobEntry(j[0]);
  log_job_manager_.StopTimer();
  return true;
}

bool Scheduler::SendLocalCopyJobToWorker(SchedulerWorker* worker,
    const ID<physical_data_id_t>& from_physical_data_id,
    const ID<physical_data_id_t>& to_physical_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    job_id_t* job_id) {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  *job_id = j[0];
  LocalCopyCommand cm_c(ID<job_id_t>(j[0]),
      ID<physical_data_id_t>(from_physical_data_id),
      ID<physical_data_id_t>(to_physical_data_id), before, after);
  dbg(DBG_SCHED, "Sending local copy job %lu to worker %lu.\n", j[0], worker->worker_id());
  server_->SendCommand(worker, &cm_c);
  return true;
}

bool Scheduler::SendCopyReceiveJobToWorker(SchedulerWorker* worker,
    const physical_data_id_t& physical_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    job_id_t* job_id) {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  *job_id = j[0];
  RemoteCopyReceiveCommand cm_r(ID<job_id_t>(j[0]),
      ID<physical_data_id_t>(physical_data_id), before, after);
  dbg(DBG_SCHED, "Sending remote copy receive job %lu to worker %lu.\n", j[0], worker->worker_id());
  server_->SendCommand(worker, &cm_r);
  return true;
}


bool Scheduler::SendCopySendJobToWorker(SchedulerWorker* worker,
    const job_id_t& receive_job_id, const physical_data_id_t& physical_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    job_id_t* job_id) {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  *job_id = j[0];
  RemoteCopySendCommand cm_s(ID<job_id_t>(j[0]),
      ID<job_id_t>(receive_job_id), ID<physical_data_id_t>(physical_data_id),
      ID<worker_id_t>(worker->worker_id()),
      worker->ip(), ID<port_t>(worker->port()),
      before, after);
  server_->SendCommand(worker, &cm_s);
  return true;
}

bool Scheduler::AssignJob(JobEntry* job) {
  log_assign_.ResumeTimer();

  SchedulerWorker* worker;
  GetWorkerToAssignJob(job, worker);

//   {
//     Log log;
//     log.StartTimer();
//     IDSet<logical_data_id_t>::ConstIter it = job->union_set_p()->begin();
//     for (; it != job->union_set_p()->end(); ++it) {
//       JobEntryList list;
//       data_version_t version;
//       job->vmap_read()->query_entry(*it, &version);
//       VersionedLogicalData vld(*it, version);
//       job_manager_->GetJobsNeedDataVersion(&list, vld);
//     }
//     log.StopTimer();
//     std::cout << "OMID: " << job->job_name() << " " << log.timer() << std::endl;
//   }

  bool prepared_data = true;
  IDSet<logical_data_id_t>::ConstIter it;
  for (it = job->union_set_p()->begin(); it != job->union_set_p()->end(); ++it) {
    if (!PrepareDataForJobAtWorker(job, worker, *it)) {
      prepared_data = false;
      break;
    }
  }
  if (prepared_data) {
    log_job_manager_.ResumeTimer();
    job_manager_->UpdateJobBeforeSet(job);
    log_job_manager_.StopTimer();
    SendComputeJobToWorker(worker, job);

    job_manager_->NotifyJobAssignment(job);

    static bool got_stamp = false;
    std::string jname = job->job_name();
    if (jname == "update_ghost_velocities" && !got_stamp) {
      log_.StopTimer();
      std::cout << "STAMP FIRST ASSIGNMENT LATENCY: " << log_.timer() << std::endl;
      got_stamp = true;
    }

    load_balancer_->NotifyJobAssignment(job, worker);

    log_assign_.StopTimer();
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: could not assign job %s (id: %lu).\n", job->job_name().c_str(), job->job_id()); // NOLINT
    log_assign_.StopTimer();
    exit(-1);
    return false;
  }
}

size_t Scheduler::AssignReadyJobs() {
  size_t count = 0;
  JobEntryList list;
  log_job_manager_.ResumeTimer();
  job_manager_->GetJobsReadyToAssign(&list, (size_t)(MAX_JOB_TO_ASSIGN));
  log_job_manager_.StopTimer();
  JobEntryList::iterator iter;
  for (iter = list.begin(); iter != list.end(); ++iter) {
    JobEntry* job = *iter;
    if (AssignJob(job)) {
      ++count;
    }
  }
  return count;
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
  data_manager_ = new DataManager();
}

void Scheduler::SetupJobManager() {
  job_manager_ = new JobManager();
  job_manager_->set_ldo_map_p(data_manager_->ldo_map_p());
}

void Scheduler::SetupLoadBalancer() {
  load_balancer_ = new LoadBalancer();
  load_balancer_->set_job_manager(job_manager_);
  load_balancer_->set_data_manager(data_manager_);
}




void Scheduler::LoadWorkerCommands() {
  // std::stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");   // NOLINT
  worker_command_table_.push_back(new SpawnComputeJobCommand());
  worker_command_table_.push_back(new SpawnCopyJobCommand());
  worker_command_table_.push_back(new DefineDataCommand());
  worker_command_table_.push_back(new HandshakeCommand());
  worker_command_table_.push_back(new JobDoneCommand());
  worker_command_table_.push_back(new DefinePartitionCommand());
  worker_command_table_.push_back(new TerminateCommand());
  worker_command_table_.push_back(new ProfileCommand());
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
