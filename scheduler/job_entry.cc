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
}

JobEntry::JobEntry(const JobType& job_type,
    const std::string& job_name,
    const job_id_t& job_id,
    const IDSet<logical_data_id_t>& read_set,
    const IDSet<logical_data_id_t>& write_set,
    const IDSet<job_id_t>& before_set,
    const IDSet<job_id_t>& after_set,
    const job_id_t& parent_job_id,
    const Parameter& params,
    const bool& sterile)
  : job_type_(job_type),
  job_name_(job_name), job_id_(job_id),
  read_set_(read_set), write_set_(write_set),
  before_set_(before_set), after_set_(after_set),
  parent_job_id_(parent_job_id), params_(params),
  sterile_(sterile) {
    Initialize();
    union_set_.insert(read_set_);
    union_set_.insert(write_set_);
    partial_versioned_ = false;
    versioned_ = false;
    assigned_ = false;
    done_ = false;
    future_ = false;
}

JobEntry::JobEntry(const JobType& job_type,
    const std::string& job_name,
    const job_id_t& job_id,
    const job_id_t& parent_job_id,
    const bool& sterile,
    const bool& versioned,
    const bool& assigned)
  : job_type_(job_type),
  job_name_(job_name), job_id_(job_id),
  parent_job_id_(parent_job_id),
  sterile_(sterile), versioned_(versioned),
  assigned_(assigned) {
    Initialize();
    partial_versioned_ = versioned;
    done_ = false;
    future_ = false;
}

JobEntry::JobEntry(const job_id_t& job_id)
  : job_id_(job_id) {
    Initialize();
    sterile_ = false;
    partial_versioned_ = false;
    versioned_ = false;
    assigned_ = false;
    done_ = false;
    future_ = true;
}

void JobEntry::Initialize() {
//   static boost::shared_ptr<nimbus::VersionTable> empty_vtable =
//     boost::shared_ptr<nimbus::VersionTable>(new
//         nimbus::VersionTable(NIMBUS_EMPTY_VERSION_TABLE_ID));
//   vtable_in_ = empty_vtable;
//   vtable_out_ = empty_vtable;

  static boost::shared_ptr<AncestorChain> empty_chain =
    boost::shared_ptr<AncestorChain>(new AncestorChain());

  ancestor_chain_ = empty_chain;
  ancestor_chain_to_pass_ = empty_chain;
}

JobEntry::~JobEntry() {
}

JobType JobEntry::job_type() {
  return job_type_;
}

std::string JobEntry::job_name() {
  return job_name_;
}

job_id_t JobEntry::job_id() {
  return job_id_;
}

IDSet<logical_data_id_t> JobEntry::read_set() {
  return read_set_;
}

IDSet<logical_data_id_t> JobEntry::write_set() {
  return write_set_;
}

const IDSet<logical_data_id_t>* JobEntry::read_set_p() {
  return &read_set_;
}

const IDSet<logical_data_id_t>* JobEntry::write_set_p() {
  return &write_set_;
}

IDSet<logical_data_id_t> JobEntry::union_set() {
  return union_set_;
}

IDSet<job_id_t> JobEntry::before_set() {
  return before_set_;
}

IDSet<job_id_t> JobEntry::after_set() {
  return after_set_;
}

job_id_t JobEntry::parent_job_id() {
  return parent_job_id_;
}

Parameter JobEntry::params() {
  return params_;
}

JobEntry::VersionTable JobEntry::version_table_in() {
  return version_table_in_;
}

JobEntry::VersionTable JobEntry::version_table_out() {
  return version_table_out_;
}

// boost::shared_ptr<nimbus::VersionTable> JobEntry::vtable_in() {
//   return vtable_in_;
// }

// boost::shared_ptr<nimbus::VersionTable> JobEntry::vtable_out() {
//   return vtable_out_;
// }

boost::shared_ptr<VersionMap> JobEntry::vmap_read_in() {
  return vmap_read_in_;
}

boost::shared_ptr<VersionMap> JobEntry::vmap_write_out() {
  return vmap_write_out_;
}

boost::shared_ptr<AncestorChain> JobEntry::ancestor_chain() {
  return ancestor_chain_;
}

boost::shared_ptr<AncestorChain> JobEntry::ancestor_chain_to_pass() {
  return ancestor_chain_to_pass_;
}

const JobEntry::VersionTable* JobEntry::version_table_in_p() {
  return &version_table_in_;
}

const JobEntry::VersionTable* JobEntry::version_table_out_p() {
  return &version_table_out_;
}

data_version_t JobEntry::version_table_in_query(logical_data_id_t l_id) {
  return version_table_in_[l_id];
}

data_version_t JobEntry::version_table_out_query(logical_data_id_t l_id) {
  return version_table_out_[l_id];
}

JobEntry::PhysicalTable JobEntry::physical_table() {
  return physical_table_;
}

IDSet<job_id_t> JobEntry::jobs_passed_versions() {
  return jobs_passed_versions_;
}

IDSet<job_id_t> JobEntry::need_set() {
  IDSet<job_id_t> need = before_set_;
  need.insert(parent_job_id_);
  need.remove(jobs_passed_versions_);
  return need;
}

bool JobEntry::sterile() {
  return sterile_;
}

bool JobEntry::partial_versioned() {
  return partial_versioned_;
}

bool JobEntry::versioned() {
  return versioned_;
}

bool JobEntry::assigned() {
  return assigned_;
}

bool JobEntry::done() {
  return done_;
}

bool JobEntry::future() {
  return future_;
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
  read_set_ = read_set;
}

void JobEntry::set_write_set(IDSet<logical_data_id_t> write_set) {
  union_set_.remove(write_set_);
  union_set_.insert(write_set);
  write_set_ = write_set;
}

void JobEntry::set_before_set(IDSet<job_id_t> before_set) {
  before_set_ = before_set;
}

void JobEntry::set_after_set(IDSet<job_id_t> after_set) {
  after_set_ = after_set;
}

void JobEntry::set_parent_job_id(job_id_t parent_job_id) {
  parent_job_id_ = parent_job_id;
}

void JobEntry::set_params(Parameter params) {
  params_ = params;
}

void JobEntry::set_version_table_in(VersionTable version_table) {
  version_table_in_ = version_table;
}

void JobEntry::set_version_table_in_entry(logical_data_id_t l_id, data_version_t version) {
  version_table_in_[l_id] = version;
}

void JobEntry::set_version_table_out(VersionTable version_table) {
  version_table_out_ = version_table;
}

void JobEntry::set_version_table_out_entry(logical_data_id_t l_id, data_version_t version) {
  version_table_out_[l_id] = version;
}

// void JobEntry::set_vtable_in(boost::shared_ptr<nimbus::VersionTable> vtable_in) {
//   vtable_in_ = vtable_in;
// }

// void JobEntry::set_vtable_out(boost::shared_ptr<nimbus::VersionTable> vtable_out) {
//   vtable_out_ = vtable_out;
// }


void JobEntry::set_vmap_read_in(boost::shared_ptr<VersionMap> vmap_read_in) {
  vmap_read_in_ = vmap_read_in;
}

void JobEntry::set_vmap_write_out(boost::shared_ptr<VersionMap> vmap_write_out) {
  vmap_write_out_ = vmap_write_out;
}

void JobEntry::set_ancestor_chain(boost::shared_ptr<AncestorChain> ancestor_chain) {
  ancestor_chain_ = ancestor_chain;
}

void JobEntry::set_ancestor_chain_to_pass(boost::shared_ptr<AncestorChain> ancestor_chain_to_pass) {
  ancestor_chain_to_pass_ = ancestor_chain_to_pass;
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

bool JobEntry::build_version_table_out_and_check() {
  version_table_out_ = version_table_in_;
  IDSet<logical_data_id_t>::ConstIter iter_data;

  for (iter_data = read_set_.begin(); iter_data != read_set_.end(); ++iter_data) {
    if (version_table_in_.count(*iter_data) == 0) {
      dbg(DBG_ERROR, "ERROR: parent and before set could not resolve read id %lu.\n", *iter_data);
      return false;
    }
  }
  for (iter_data = write_set_.begin(); iter_data != write_set_.end(); ++iter_data) {
    if (version_table_in_.count(*iter_data) == 0) {
      dbg(DBG_ERROR, "ERROR: parent and before set could not resolve write id %lu.\n", *iter_data);
      return false;
    } else {
      ++version_table_out_[*iter_data];
    }
  }

  return true;
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




