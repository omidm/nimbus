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
  * This is the base class that serves the scheduler regarding job assignment
  * queries in the cluster. It tries to minimize the completion of simulation
  * by reducing the cost of communication by locality aware data placement and
  * mitigate the effect of stragglers in the system by adapting the job
  * assignment strategies to the dynamic changes of the cloud.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "scheduler/load_balancer.h"

#define LB_UPDATE_RATE 100
#define JOB_ASSIGNER_THREAD_NUM 5

namespace nimbus {

LoadBalancer::LoadBalancer() {
  Initialize();
}

LoadBalancer::LoadBalancer(ClusterMap* cluster_map) {
  Initialize();
  cluster_map_ = cluster_map;
}

void LoadBalancer::Initialize() {
  worker_num_ = 0;
  global_region_ = GeometricRegion(0, 0, 0, 0, 0, 0);
  update_ = false;
  init_phase_ = true;
  blame_counter_ = 0;
  server_ = NULL;
  id_maker_ = NULL;
  cluster_map_ = NULL;
  job_manager_ = NULL;
  data_manager_ = NULL;
  log_.set_file_name("load_balancer_log");
}

LoadBalancer::~LoadBalancer() {
}

ClusterMap* LoadBalancer::cluster_map() {
  return cluster_map_;
}

void LoadBalancer::set_cluster_map(ClusterMap* cluster_map) {
  cluster_map_ = cluster_map;
}

void LoadBalancer::set_server(SchedulerServer *server) {
  server_ = server;
}

void LoadBalancer::set_id_maker(IDMaker *id_maker) {
  id_maker_ = id_maker;
}

void LoadBalancer::set_job_manager(JobManager *job_manager) {
  job_manager_ = job_manager;
}

void LoadBalancer::set_data_manager(DataManager *data_manager) {
  data_manager_ = data_manager;
}

void LoadBalancer::Run() {
  for (size_t i = 0; i < JOB_ASSIGNER_THREAD_NUM; ++i) {
    job_assigner_threads_.push_back(
        new boost::thread(boost::bind(&LoadBalancer::JobAssignerThread, this)));
  }

  while (true) {
    boost::unique_lock<boost::recursive_mutex> update_lock(update_mutex_);
    while (!update_) {
      update_cond_.wait(update_lock);
    }
    update_ = false;
    update_cond_.notify_all();

    boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
    boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

    if (worker_num_ != region_map_.table_size() ||
        global_region_ != data_manager_->global_bounding_region()) {
      if (!data_manager_->initialized_global_bounding_region()) {
        continue;
      } else {
        InitializeRegionMap();
      }
    } else {
      UpdateRegionMap();
    }
  }
}


void LoadBalancer::JobAssignerThread() {
  while (true) {
    JobEntry *job;
    {
      boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);

      while (job_queue_.size() == 0) {
        job_queue_cond_.wait(job_queue_lock);
      }

      JobEntryList::iterator iter = job_queue_.begin();
      job = *iter;
      job_queue_.erase(iter);
    }

    if (!AssignJob(job)) {
      dbg(DBG_ERROR, "ERROR: LoadBalancer: could not assign job %lu.\n", job->job_id());
      exit(-1);
    }

    job_queue_cond_.notify_all();
  }
}

void LoadBalancer::AssignJobs(const JobEntryList& list) {
  boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);
  assert(job_queue_.size() == 0);
  job_queue_ = list;
  job_queue_cond_.notify_all();

  while (job_queue_.size() > 0) {
    job_queue_cond_.wait(job_queue_lock);
  }
}

bool LoadBalancer::AssignJob(JobEntry *job) {
  SchedulerWorker* worker;
  GetWorkerToAssignJob(job, worker);

  job_manager_->ResolveJobDataVersions(job);

  bool prepared_data = true;
  IDSet<logical_data_id_t>::ConstIter it;
  for (it = job->union_set_p()->begin(); it != job->union_set_p()->end(); ++it) {
    if (!PrepareDataForJobAtWorker(job, worker, *it)) {
      prepared_data = false;
      break;
    }
  }

  if (prepared_data) {
    job_manager_->UpdateJobBeforeSet(job);
    SendComputeJobToWorker(worker, job);

    job_manager_->NotifyJobAssignment(job, worker);

    static bool loop_begins = true;
    std::string jname = job->job_name();
    if (jname == "update_ghost_velocities" && loop_begins) {
      std::cout << "STAMP: FIRST ASSIGNMENT LATENCY: " << log_.timer() << std::endl;
      loop_begins = false;
    }
    if (jname == "projection_main" && !loop_begins) {
      std::cout << "STAMP: ALL ASSIGNMENT LATENCY: " << log_.timer() << std::endl;
      loop_begins = true;
    }

    NotifyJobAssignment(job, worker);

    return true;
  }

  return false;
}


bool LoadBalancer::PrepareDataForJobAtWorker(JobEntry* job,
    SchedulerWorker* worker, logical_data_id_t l_id) {
  bool reading = job->read_set_p()->contains(l_id);
  bool writing = job->write_set_p()->contains(l_id);
  assert(reading || writing);

  LogicalDataObject* ldo =
    const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(l_id));

  boost::unique_lock<boost::mutex> lock(ldo->mutex());

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

  if (!reading) {
    PhysicalData target_instance;
    GetFreeDataAtWorker(worker, ldo, &target_instance);

    if (job_manager_->CausingUnwantedSerialization(job, l_id, target_instance)) {
      dbg(DBG_SCHED, "Causing unwanted serialization for data %lu.\n", l_id);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }

  PhysicalDataVector instances_at_worker;
  PhysicalDataVector instances_in_system;
  data_manager_->InstancesByWorkerAndVersion(
      ldo, worker->worker_id(), version, &instances_at_worker);
  data_manager_->InstancesByVersion(ldo, version, &instances_in_system);

  JobEntryList list;
  VersionedLogicalData vld(l_id, version);
  job_manager_->GetJobsNeedDataVersion(&list, vld);
  assert(list.size() >= 1);
  bool writing_needed_version = (list.size() > 1) && writing;


  if (instances_at_worker.size() > 1) {
    PhysicalData target_instance;

    bool found = false;
    PhysicalDataVector::iterator iter;
    for (iter = instances_at_worker.begin(); iter != instances_at_worker.end(); iter++) {
      if (!job_manager_->CausingUnwantedSerialization(job, l_id, *iter)) {
        target_instance = *iter;
        found = true;
        break;
      }
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

    if (!job_manager_->CausingUnwantedSerialization(job, l_id, instances_at_worker[0])) {
      target_instance = instances_at_worker[0];
    } else {
      dbg(DBG_SCHED, "Avoiding unwanted serialization for data %lu (2).\n", l_id);
      GetFreeDataAtWorker(worker, ldo, &target_instance);
      LocalCopyData(worker, ldo, &instances_at_worker[0], &target_instance);
    }

    AllocateLdoInstanceToJob(job, ldo, target_instance);
    return true;
  }


  if ((instances_at_worker.size() == 1) && writing_needed_version) {
    PhysicalData target_instance;

    if (!job_manager_->CausingUnwantedSerialization(job, l_id, instances_at_worker[0])) {
      target_instance = instances_at_worker[0];
      PhysicalData copy_data;
      GetFreeDataAtWorker(worker, ldo, &copy_data);
      LocalCopyData(worker, ldo, &target_instance, &copy_data);
    } else {
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


bool LoadBalancer::AllocateLdoInstanceToJob(JobEntry* job,
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

  data_manager_->RemovePhysicalInstance(ldo, pd);
  data_manager_->AddPhysicalInstance(ldo, pd_new);

  return true;
}

size_t LoadBalancer::GetObsoleteLdoInstancesAtWorker(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalDataVector* dest) {
  size_t count = 0;
  dest->clear();
  PhysicalDataVector pv;
  data_manager_->InstancesByWorker(ldo, worker->worker_id(), &pv);
  PhysicalDataVector::iterator iter = pv.begin();
  for (; iter != pv.end(); ++iter) {
    JobEntryList list;
    VersionedLogicalData vld(ldo->id(), iter->version());
    if (job_manager_->GetJobsNeedDataVersion(&list, vld) == 0) {
      dest->push_back(*iter);
      ++count;
    }
  }
  return count;
}

bool LoadBalancer::CreateDataAtWorker(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalData* created_data) {
  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 1);
  std::vector<physical_data_id_t> d;
  id_maker_->GetNewPhysicalDataID(&d, 1);
  IDSet<job_id_t> before;

  // Update the job table.
  job_manager_->AddCreateDataJobEntry(j[0]);

  // Update data table.
  IDSet<job_id_t> list_job_read;
  list_job_read.insert(j[0]);  // if other job wants to write, waits for creation.
  PhysicalData p(d[0], worker->worker_id(), NIMBUS_INIT_DATA_VERSION, list_job_read, j[0]);
  data_manager_->AddPhysicalInstance(ldo, p);

  // send the create command to worker.
  job_manager_->UpdateBeforeSet(&before);
  CreateDataCommand cm(ID<job_id_t>(j[0]),
                       ldo->variable(),
                       ID<logical_data_id_t>(ldo->id()),
                       ID<physical_data_id_t>(d[0]),
                       before);
  server_->SendCommand(worker, &cm);

  *created_data = p;

  return true;
}

bool LoadBalancer::RemoteCopyData(SchedulerWorker* from_worker,
    SchedulerWorker* to_worker, LogicalDataObject* ldo,
    PhysicalData* from_data, PhysicalData* to_data) {
  assert(from_worker->worker_id() == from_data->worker());
  assert(to_worker->worker_id() == to_data->worker());

  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 2);
  job_id_t receive_id = j[0];
  job_id_t send_id = j[1];
  IDSet<job_id_t> before;

  // Receive part

  // Update the job table.
  job_manager_->AddRemoteCopyReceiveJobEntry(receive_id);

  // Update data table.
  PhysicalData to_data_new = *to_data;
  to_data_new.set_version(from_data->version());
  to_data_new.set_last_job_write(receive_id);
  to_data_new.clear_list_job_read();
  data_manager_->RemovePhysicalInstance(ldo, *to_data);
  data_manager_->AddPhysicalInstance(ldo, to_data_new);

  // send remote copy receive job to worker.
  before.clear();
  before.insert(to_data->list_job_read());
  before.insert(to_data->last_job_write());
  job_manager_->UpdateBeforeSet(&before);
  RemoteCopyReceiveCommand cm_r(ID<job_id_t>(receive_id),
                                ID<physical_data_id_t>(to_data->id()),
                                before);
  server_->SendCommand(to_worker, &cm_r);


  // Send Part.

  // Update the job table.
  job_manager_->AddRemoteCopySendJobEntry(send_id);

  // Update data table.
  PhysicalData from_data_new = *from_data;
  from_data_new.add_to_list_job_read(send_id);
  data_manager_->RemovePhysicalInstance(ldo, *from_data);
  data_manager_->AddPhysicalInstance(ldo, from_data_new);

  // send remote copy send command to worker.
  before.clear();
  before.insert(from_data->last_job_write());
  job_manager_->UpdateBeforeSet(&before);
  RemoteCopySendCommand cm_s(ID<job_id_t>(send_id),
                             ID<job_id_t>(receive_id),
                             ID<physical_data_id_t>(from_data->id()),
                             ID<worker_id_t>(to_worker->worker_id()),
                             to_worker->ip(),
                             ID<port_t>(to_worker->port()),
                             before);
  server_->SendCommand(from_worker, &cm_s);


  *from_data = from_data_new;
  *to_data = to_data_new;

  return true;
}

bool LoadBalancer::LocalCopyData(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalData* from_data, PhysicalData* to_data) {
  assert(worker->worker_id() == from_data->worker());
  assert(worker->worker_id() == to_data->worker());

  std::vector<job_id_t> j;
  id_maker_->GetNewJobID(&j, 1);
  IDSet<job_id_t> before;

  // Update the job table.
  job_manager_->AddLocalCopyJobEntry(j[0]);

  // Update data table.
  PhysicalData from_data_new = *from_data;
  from_data_new.add_to_list_job_read(j[0]);
  data_manager_->RemovePhysicalInstance(ldo, *from_data);
  data_manager_->AddPhysicalInstance(ldo, from_data_new);

  PhysicalData to_data_new = *to_data;
  to_data_new.set_version(from_data->version());
  to_data_new.set_last_job_write(j[0]);
  to_data_new.clear_list_job_read();
  data_manager_->RemovePhysicalInstance(ldo, *to_data);
  data_manager_->AddPhysicalInstance(ldo, to_data_new);

  // send local copy command to worker.
  before.insert(to_data->list_job_read());
  before.insert(to_data->last_job_write());
  before.insert(from_data->last_job_write());
  job_manager_->UpdateBeforeSet(&before);
  LocalCopyCommand cm_c(ID<job_id_t>(j[0]),
                        ID<physical_data_id_t>(from_data->id()),
                        ID<physical_data_id_t>(to_data->id()),
                        before);
  server_->SendCommand(worker, &cm_c);

  *from_data = from_data_new;
  *to_data = to_data_new;

  return true;
}

bool LoadBalancer::GetFreeDataAtWorker(SchedulerWorker* worker,
    LogicalDataObject* ldo, PhysicalData* free_data) {
  PhysicalDataVector obsolete_instances;
  if (GetObsoleteLdoInstancesAtWorker(worker, ldo, &obsolete_instances) > 0) {
    *free_data = obsolete_instances[0];
    return true;
  }

  // Since there are no obsoletes, go ahead and create a new one.
  return CreateDataAtWorker(worker, ldo, free_data);
}







bool LoadBalancer::SendComputeJobToWorker(SchedulerWorker* worker, JobEntry* job) {
  if (job->job_type() == JOB_COMP) {
    ID<job_id_t> job_id(job->job_id());
    ID<job_id_t> future_job_id(job->future_job_id());
    IDSet<physical_data_id_t> read_set, write_set;
    // TODO(omidm): check the return value of the following methods.
    job->GetPhysicalReadSet(&read_set);
    job->GetPhysicalWriteSet(&write_set);
    ComputeJobCommand cm(job->job_name(),
                         job_id,
                         read_set,
                         write_set,
                         job->before_set(),
                         job->after_set(),
                         future_job_id,
                         job->sterile(),
                         job->params());
    dbg(DBG_SCHED, "Sending compute job %lu to worker %lu.\n", job->job_id(), worker->worker_id());
    server_->SendCommand(worker, &cm);
    return true;
  } else {
    dbg(DBG_ERROR, "Job with id %lu is not a compute job.\n", job->job_id());
    return false;
  }
}






bool LoadBalancer::GetWorkerToAssignJob(
    JobEntry *job, SchedulerWorker*& worker) {
  Log log;
  log.StartTimer();

  boost::unique_lock<boost::recursive_mutex> update_lock(update_mutex_);
  while (update_) {
    update_cond_.wait(update_lock);
  }

  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

  if ((worker_num_ != region_map_.table_size()) ||
      (global_region_ != data_manager_->global_bounding_region())) {
    if (data_manager_->initialized_global_bounding_region()) {
      InitializeRegionMap();
    }
  }

  assert(worker_map_.size() > 0);
  assert(worker_num_ > 0);

  GeometricRegion region;
  bool got_region = job->GetWriteSetRegion(data_manager_, &region);

  if (init_phase_ || !got_region) {
    worker = worker_map_.begin()->second;
  } else {
    assert(worker_num_ == region_map_.table_size());
    worker_id_t w_id;

    if (!region_map_.QueryWorkerWithMostOverlap(&region, &w_id)) {
      dbg(DBG_ERROR, "ERROR: LoadBalancer: could not find worker for assigning job %lu.\n", job->job_id()); // NOLINT
      return false;
    }

    worker = worker_map_[w_id];
  }

  log.StopTimer();
  std::cout
    << "Picked worker: " << worker->worker_id()
    << " for job: " << job->job_name()
    << " took: " << log.timer()
    << " for union set size of: " << job->union_set_p()->size()
    << std::endl;


  return true;
}


void LoadBalancer::NotifyJobAssignment(
    const JobEntry *job, const SchedulerWorker* worker) {
  double time = log_.GetTime();

  if (job->job_type() != JOB_COMP) {
    return;
  }

  JobProfile *job_profile =
    new JobProfile(
        job->job_type(),
        job->job_name(),
        job->job_id(),
        job->parent_job_id(),
        worker->worker_id(),
        job->sterile());

  job_profile->set_assign_time(time);
  job_profile->set_assigned(true);


  Vertex<JobEntry, job_id_t>* vertex;
  job_manager_->job_graph_p()->GetVertex(job->job_id(), &vertex);

  typename Edge<JobEntry, job_id_t>::Iter iter;
  for (iter = vertex->incoming_edges()->begin(); iter != vertex->incoming_edges()->end(); ++iter) {
    JobEntry *j = iter->second->start_vertex()->entry();
    if (!j->done()) {
      job_profile->waiting_set_p()->insert(j->job_id());
    } else {
      /*
      JobHistory::iterator it = job_history_.find(j->job_id());
      if (it != job_history_.end()) {
        JobProfile *jp = it->second;
        job_profile->add_log_entry(
            jp->worker_id(), jp->job_id(), jp->done_time());
      } else {
        dbg(DBG_WARN, "WARNING: Load balancer, could not find done job in job history.");
        exit(-1);
      }
      */
    }
  }

  if (job_profile->waiting_set_p()->size() == 0) {
    job_profile->set_ready_time(time);
    job_profile->set_ready(true);
  }

  boost::unique_lock<boost::recursive_mutex> lock(job_history_mutex_);
  job_history_[job->job_id()] = job_profile;
}

void LoadBalancer::NotifyJobDone(const JobEntry *job) {
  double time = log_.GetTime();

  if (job->job_type() != JOB_COMP) {
    return;
  }

  assert(job->done());
  done_jobs_.push_back(job->job_id());

  JobHistory::iterator it = job_history_.find(job->job_id());
  assert(it != job_history_.end());
  JobProfile *job_profile = it->second;

  job_profile->set_done_time(time);
  job_profile->set_done(true);
  job_profile->set_execute_duration(time - job_profile->ready_time());

  Vertex<JobEntry, job_id_t>* vertex;
  job_manager_->job_graph_p()->GetVertex(job->job_id(), &vertex);

  typename Edge<JobEntry, job_id_t>::Iter iter;
  for (iter = vertex->outgoing_edges()->begin(); iter != vertex->outgoing_edges()->end(); ++iter) {
    JobEntry *j = iter->second->end_vertex()->entry();
    it = job_history_.find(j->job_id());
    if (it != job_history_.end()) {
      assert(j->assigned());
      JobProfile *jp = it->second;
      assert(jp->assigned());
      jp->add_log_entry(
          job_profile->worker_id(), job->job_id(), job->job_name(), time);
      jp->waiting_set_p()->remove(job->job_id());
      if (jp->waiting_set_p()->size() == 0) {
        jp->set_ready_time(time);
        jp->set_ready(true);
      }
    }
  }

  // log_.WriteToFile(job_profile->Print());

  worker_id_t blamed_worker_id;
  if (job_profile->FindBlamedWorker(&blamed_worker_id)) {
    boost::unique_lock<boost::recursive_mutex> straggler_map_lock(straggler_map_mutex_);
    straggler_map_.AddRecord(job->assigned_worker(), blamed_worker_id);
    std::cout << "STRAGGLER ADD RECORD: job name: " << job->job_name()
              << " worker: " << job->assigned_worker()
              << " blamed: " << blamed_worker_id << std::endl;

    ++blame_map_[blamed_worker_id];

    blame_counter_++;
    if (blame_counter_ > LB_UPDATE_RATE) {
      blame_counter_ = 0;
      update_ = true;
      update_cond_.notify_all();
    }
  }
}


void LoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  boost::unique_lock<boost::recursive_mutex> update_lock(update_mutex_);
  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);

  worker_id_t worker_id = worker->worker_id();
  WorkerMapIter iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
    worker_num_ = worker_map_.size();
    update_ = true;
    update_cond_.notify_all();
  } else {
    dbg(DBG_ERROR, "ERROR: LoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
  }
}

void LoadBalancer::InitializeRegionMap() {
  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

  std::vector<worker_id_t> worker_ids;
  WorkerMapIter iter = worker_map_.begin();
  for (; iter != worker_map_.end(); ++iter) {
    worker_ids.push_back(iter->first);
  }
  assert(worker_num_ > 0);
  global_region_ = data_manager_->global_bounding_region();

  region_map_.Initialize(worker_ids, global_region_);
  log_.WriteToFile(region_map_.Print());

  init_phase_ = false;
}

void LoadBalancer::UpdateRegionMap() {
  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> straggler_map_lock(straggler_map_mutex_);

  worker_id_t fast, slow;
  if (straggler_map_.GetMostImbalanceWorkers(&fast, &slow)) {
    std::cout << "LOAD BALANCER: fast worker: " << fast
              << ", slow worker: " << slow << std::endl;
    if (region_map_.BalanceRegions(fast, slow)) {
      straggler_map_.ClearRecords();
      log_.WriteToFile(region_map_.Print());
    }
  }

  worker_id_t worst_worker = 0;
  size_t count = 0;
  std::map<worker_id_t, size_t>::iterator iter = blame_map_.begin();
  for (; iter != blame_map_.end(); ++iter) {
    if (iter->second > count) {
      count = iter->second;
      worst_worker = iter->first;
    }
  }
  std::cout << "LOAD BALANCER: worst worker: " << worst_worker << std::endl;
  blame_map_.clear();
}

}  // namespace nimbus
