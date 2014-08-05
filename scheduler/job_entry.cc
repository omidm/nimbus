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
  * Job entry in the job table of the job manager. Each entry holds the
  * meta data of the compute and copy jobs received by the scheduler.   
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/job_entry.h"

using namespace nimbus; // NOLINT

JobEntry::JobEntry() {
  Initialize();
  job_depth_ = NIMBUS_INIT_JOB_DEPTH;
  sterile_ = false;
  partial_versioned_ = false;
  versioned_ = false;
  assigned_ = false;
  done_ = false;
  future_ = false;
  region_valid_ = false;
}

void JobEntry::Initialize() {
    vmap_read_ = boost::shared_ptr<VersionMap>(new VersionMap());
    vmap_write_ = boost::shared_ptr<VersionMap>(new VersionMap());
    vmap_partial_ = boost::shared_ptr<VersionMap>(new VersionMap());
    meta_before_set_ =  boost::shared_ptr<MetaBeforeSet>(new MetaBeforeSet());
/*
  static boost::shared_ptr<VersionMap> empty_vmap =
    boost::shared_ptr<VersionMap>(new VersionMap());
  vmap_read_ = empty_vmap;
  vmap_write_ = empty_vmap;


  static boost::shared_ptr<MetaBeforeSet> empty_meta_before_set =
    boost::shared_ptr<MetaBeforeSet> (new MetaBeforeSet());
  meta_before_set_ = empty_meta_before_set;
*/
}

JobEntry::~JobEntry() {
}

JobType JobEntry::job_type() const {
  return job_type_;
}

std::string JobEntry::job_name() const {
  return job_name_;
}

job_id_t JobEntry::job_id() const {
  return job_id_;
}

IDSet<logical_data_id_t> JobEntry::read_set() const {
  return read_set_;
}

IDSet<logical_data_id_t> JobEntry::write_set() const {
  return write_set_;
}
IDSet<logical_data_id_t> JobEntry::union_set() const {
  return union_set_;
}

IDSet<job_id_t> JobEntry::before_set() const {
  return before_set_;
}

IDSet<job_id_t> JobEntry::after_set() const {
  return after_set_;
}

job_id_t JobEntry::parent_job_id() const {
  return parent_job_id_;
}

Parameter JobEntry::params() const {
  return params_;
}

boost::shared_ptr<VersionMap> JobEntry::vmap_read() const {
  return vmap_read_;
}

boost::shared_ptr<VersionMap> JobEntry::vmap_write() const {
  return vmap_write_;
}

boost::shared_ptr<VersionMap> JobEntry::vmap_partial() const {
  return vmap_partial_;
}

boost::shared_ptr<MetaBeforeSet> JobEntry::meta_before_set() const {
  return meta_before_set_;
}

job_depth_t JobEntry::job_depth() const {
  return job_depth_;
}

JobEntry::PhysicalTable JobEntry::physical_table() const {
  return physical_table_;
}

IDSet<job_id_t> JobEntry::jobs_passed_versions() const {
  return jobs_passed_versions_;
}

IDSet<job_id_t> JobEntry::need_set() const {
  IDSet<job_id_t> need = before_set_;
  need.insert(parent_job_id_);
  need.remove(jobs_passed_versions_);
  return need;
}

bool JobEntry::sterile() const {
  return sterile_;
}

bool JobEntry::partial_versioned() const {
  return partial_versioned_;
}

bool JobEntry::versioned() const {
  return versioned_;
}

bool JobEntry::assigned() const {
  return assigned_;
}

bool JobEntry::done() const {
  return done_;
}

bool JobEntry::future() const {
  return future_;
}

const IDSet<logical_data_id_t>* JobEntry::read_set_p() const {
  return &read_set_;
}

const IDSet<logical_data_id_t>* JobEntry::write_set_p() const {
  return &write_set_;
}

const IDSet<logical_data_id_t>* JobEntry::union_set_p() const {
  return &union_set_;
}

const IDSet<job_id_t>* JobEntry::before_set_p() const {
  return &before_set_;
}

IDSet<job_id_t>* JobEntry::before_set_p() {
  return &before_set_;
}

void JobEntry::set_job_type(JobType job_type) {
  job_type_ = job_type;
}

void JobEntry::set_job_name(std::string job_name) {
  job_name_ = job_name;
}

void JobEntry::set_job_id(job_id_t job_id) {
  job_id_ = job_id;
}

void JobEntry::set_read_set(IDSet<logical_data_id_t> read_set) {
  union_set_.remove(read_set_);
  union_set_.insert(read_set);
  region_valid_ = false;
  read_set_ = read_set;
}

void JobEntry::set_write_set(IDSet<logical_data_id_t> write_set) {
  union_set_.remove(write_set_);
  union_set_.insert(write_set);
  region_valid_ = false;
  write_set_ = write_set;
}

void JobEntry::set_before_set(IDSet<job_id_t> before_set, bool update_dependencies) {
  if (update_dependencies) {
    assignment_dependencies_.remove(before_set_);
    assignment_dependencies_.insert(before_set);
    versioning_dependencies_ = assignment_dependencies_;
  }
  before_set_ = before_set;
}

void JobEntry::set_after_set(IDSet<job_id_t> after_set) {
  after_set_ = after_set;
}

void JobEntry::set_parent_job_id(job_id_t parent_job_id, bool update_dependencies) {
  if (update_dependencies) {
    assignment_dependencies_.remove(parent_job_id_);
    assignment_dependencies_.insert(parent_job_id);
    versioning_dependencies_ = assignment_dependencies_;
  }
  parent_job_id_ = parent_job_id;
}

void JobEntry::set_params(Parameter params) {
  params_ = params;
}

void JobEntry::set_vmap_read(boost::shared_ptr<VersionMap> vmap_read) {
  vmap_read_ = vmap_read;
}

void JobEntry::set_vmap_write(boost::shared_ptr<VersionMap> vmap_write) {
  vmap_write_ = vmap_write;
}

void JobEntry::set_meta_before_set(boost::shared_ptr<MetaBeforeSet> meta_before_set) {
  meta_before_set_ = meta_before_set;
}

void JobEntry::set_job_depth(job_depth_t job_depth) {
  job_depth_ = job_depth;
  meta_before_set_->set_job_depth(job_depth);
}

void JobEntry::set_physical_table(PhysicalTable physical_table) {
  physical_table_ = physical_table;
}

void JobEntry::set_physical_table_entry(logical_data_id_t l_id, physical_data_id_t p_id) {
  physical_table_[l_id] = p_id;
}

void JobEntry::set_jobs_passed_versions(IDSet<job_id_t> jobs) {
  jobs_passed_versions_ = jobs;
}

void JobEntry::add_job_passed_versions(job_id_t job_id) {
  jobs_passed_versions_.insert(job_id);
}


void JobEntry::set_sterile(bool flag) {
  sterile_ = flag;
}

void JobEntry::set_partial_versioned(bool flag) {
  partial_versioned_ = flag;
}

void JobEntry::set_versioned(bool flag) {
  versioned_ = flag;
}

void JobEntry::set_assigned(bool flag) {
  assigned_ = flag;
}

void JobEntry::set_done(bool flag) {
  done_ = flag;
}

void JobEntry::set_future(bool flag) {
  future_ = flag;
}

bool JobEntry::GetPhysicalReadSet(IDSet<physical_data_id_t>* set) {
  set->clear();
  IDSet<logical_data_id_t>::IDSetIter it;
  for (it = read_set_.begin(); it != read_set_.end(); ++it) {
    if (physical_table_.count(*it) == 0)
      return false;
    set->insert(physical_table_[*it]);
  }
  return true;
}

bool JobEntry::GetPhysicalWriteSet(IDSet<physical_data_id_t>* set) {
  set->clear();
  IDSet<logical_data_id_t>::IDSetIter it;
  for (it = write_set_.begin(); it != write_set_.end(); ++it) {
    if (physical_table_.count(*it) == 0)
      return false;
    set->insert(physical_table_[*it]);
  }
  return true;
}

bool JobEntry::IsReadyToAssign() {
  return (assignment_dependencies_.size() == 0);
}

void JobEntry::remove_assignment_dependency(job_id_t job_id) {
  assignment_dependencies_.remove(job_id);
}

bool JobEntry::IsReadyForCompleteVersioning() {
  return (versioning_dependencies_.size() == 0 && !future_);
}

void JobEntry::remove_versioning_dependency(job_id_t job_id) {
  versioning_dependencies_.remove(job_id);
  meta_before_set_->InvalidateNegativeQueryCache();
}

bool JobEntry::GetRegion(DataManager *data_manager, GeometricRegion *region) {
  if (region_valid_) {
    *region = region_;
    return true;
  } else {
    if (union_set_.size() == 0) {
      return false;
    } else {
      const LogicalDataObject* ldo;
      IDSet<logical_data_id_t>::IDSetIter iter = union_set_.begin();
      ldo = data_manager->FindLogicalObject(*iter);
      region_ = GeometricRegion::GetBoundingBox(region_, *ldo->region());
      ++iter;
      for (; iter != union_set_.end(); ++iter) {
        ldo = data_manager->FindLogicalObject(*iter);
        region_ = GeometricRegion::GetBoundingBox(region_, *ldo->region());
      }
      *region = region_;
      region_valid_ = true;
      return true;
    }
  }
}


ComputeJobEntry::ComputeJobEntry(
    const std::string& job_name,
    const job_id_t& job_id,
    const IDSet<logical_data_id_t>& read_set,
    const IDSet<logical_data_id_t>& write_set,
    const IDSet<job_id_t>& before_set,
    const IDSet<job_id_t>& after_set,
    const job_id_t& parent_job_id,
    const Parameter& params,
    const bool& sterile) {
    job_type_ = JOB_COMP;
    job_name_ = job_name;
    job_id_ = job_id;
    read_set_ = read_set;
    write_set_ = write_set;
    before_set_ = before_set;
    after_set_ = after_set;
    parent_job_id_ = parent_job_id;
    params_ = params;
    sterile_ = sterile;

    union_set_.insert(read_set_);
    union_set_.insert(write_set_);

    assignment_dependencies_ = before_set;
    assignment_dependencies_.insert(parent_job_id);
    versioning_dependencies_ = assignment_dependencies_;
}

ComputeJobEntry::~ComputeJobEntry() {
}


KernelJobEntry::KernelJobEntry() {
  job_type_ = JOB_SCHED;
  job_name_ = NIMBUS_KERNEL_JOB_NAME;
  job_id_ = NIMBUS_KERNEL_JOB_ID;
  parent_job_id_ = NIMBUS_KERNEL_JOB_ID;
  sterile_ = true;
  versioned_ = true;
  assigned_ = true;

  meta_before_set_ = boost::shared_ptr<MetaBeforeSet> (new MetaBeforeSet());
  job_depth_ = NIMBUS_INIT_JOB_DEPTH;
}

KernelJobEntry::~KernelJobEntry() {
}


MainJobEntry::MainJobEntry(const job_id_t& job_id) {
  job_type_ = JOB_COMP;
  job_name_ = NIMBUS_MAIN_JOB_NAME;
  job_id_ = job_id;
  parent_job_id_ = NIMBUS_KERNEL_JOB_ID;
}

MainJobEntry::~MainJobEntry() {
}


FutureJobEntry::FutureJobEntry(const job_id_t& job_id) {
  job_type_ = JOB_FUTURE;
  job_id_ = job_id;
  future_ = true;
}

FutureJobEntry::~FutureJobEntry() {
}


LocalCopyJobEntry::LocalCopyJobEntry(const job_id_t& job_id) {
  job_type_ = JOB_COPY;
  job_name_ = NIMBUS_LOCAL_COPY_JOB_NAME;
  job_id_ = job_id;
  parent_job_id_ = NIMBUS_KERNEL_JOB_ID;
  sterile_ = true;
  versioned_ = true;
  assigned_ = true;
}

LocalCopyJobEntry::~LocalCopyJobEntry() {
}


CreateDataJobEntry::CreateDataJobEntry(const job_id_t& job_id) {
  job_type_ = JOB_CREATE;
  job_name_ = NIMBUS_CREATE_DATA_JOB_NAME;
  job_id_ = job_id;
  parent_job_id_ = NIMBUS_KERNEL_JOB_ID;
  sterile_ = true;
  versioned_ = true;
  assigned_ = true;
}

CreateDataJobEntry::~CreateDataJobEntry() {
}


RemoteCopyReceiveJobEntry::RemoteCopyReceiveJobEntry(const job_id_t& job_id) {
  job_type_ = JOB_COPY;
  job_name_ = NIMBUS_REMOTE_COPY_RECEIVE_JOB_NAME;
  job_id_ = job_id;
  parent_job_id_ = NIMBUS_KERNEL_JOB_ID;
  sterile_ = true;
  versioned_ = true;
  assigned_ = true;
}

RemoteCopyReceiveJobEntry::~RemoteCopyReceiveJobEntry() {
}


RemoteCopySendJobEntry::RemoteCopySendJobEntry(const job_id_t& job_id) {
  job_type_ = JOB_COPY;
  job_name_ = NIMBUS_REMOTE_COPY_SEND_JOB_NAME;
  job_id_ = job_id;
  parent_job_id_ = NIMBUS_KERNEL_JOB_ID;
  sterile_ = true;
  versioned_ = true;
  assigned_ = true;
}

RemoteCopySendJobEntry::~RemoteCopySendJobEntry() {
}

ID<job_id_t> RemoteCopySendJobEntry::receive_job_id() {
  return receive_job_id_;
}

ID<worker_id_t> RemoteCopySendJobEntry::to_worker_id() {
  return to_worker_id_;
}

std::string RemoteCopySendJobEntry::to_ip() {
  return to_ip_;
}

ID<port_t> RemoteCopySendJobEntry::to_port() {
  return to_port_;
}

void RemoteCopySendJobEntry::set_receive_job_id(ID<job_id_t> receive_job_id) {
  receive_job_id_ = receive_job_id;
}

void RemoteCopySendJobEntry::set_to_worker_id(ID<worker_id_t> worker_id) {
  to_worker_id_ = worker_id;
}

void RemoteCopySendJobEntry::set_to_ip(std::string ip) {
  to_ip_ = ip;
}

void RemoteCopySendJobEntry::set_to_port(ID<port_t> port) {
  to_port_ = port;
}


