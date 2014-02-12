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
}

JobEntry::JobEntry(const JobType& job_type,
    const std::string& job_name,
    const job_id_t& job_id,
    const IDSet<logical_data_id_t>& read_set,
    const IDSet<logical_data_id_t>& write_set,
    const IDSet<job_id_t>& before_set,
    const IDSet<job_id_t>& after_set,
    const job_id_t& parent_job_id,
    const Parameter& params)
  : job_type_(job_type),
  job_name_(job_name), job_id_(job_id),
  read_set_(read_set), write_set_(write_set),
  before_set_(before_set), after_set_(after_set),
  parent_job_id_(parent_job_id), params_(params) {
    IDSet<logical_data_id_t>::IDSetIter it;
    for (it = read_set_.begin(); it != read_set_.end(); ++it) {
      union_set_.insert(*it);
    }
    for (it = write_set_.begin(); it != write_set_.end(); ++it) {
      union_set_.insert(*it);
    }
    versioned_ = false;
    assigned_ = false;
    done_ = false;
}

JobEntry::JobEntry(const JobType& job_type,
    const std::string& job_name,
    const job_id_t& job_id,
    const job_id_t& parent_job_id)
  : job_type_(job_type),
  job_name_(job_name), job_id_(job_id),
  parent_job_id_(parent_job_id) {
    versioned_ = false;
    assigned_ = false;
    done_ = false;
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

JobEntry::PhysicalTable JobEntry::physical_table() {
  return physical_table_;
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

void JobEntry::set_before_set(IDSet<job_id_t> before_set) {
  before_set_ = before_set;
}

void JobEntry::set_after_set(IDSet<job_id_t> after_set) {
  after_set_ = after_set;
}

void JobEntry::set_version_table_in(VersionTable version_table) {
  version_table_in_ = version_table;
}

void JobEntry::set_version_table_out(VersionTable version_table) {
  version_table_out_ = version_table;
}

void JobEntry::set_physical_table(PhysicalTable physical_table) {
  physical_table_ = physical_table;
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




