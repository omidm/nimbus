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
  * Job profile in the job graph of load balancer.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/job_profile.h"

using namespace nimbus; // NOLINT

JobProfile::JobProfile() {
  Initialize();
}

JobProfile::JobProfile(
    const JobType& job_type,
    const std::string& job_name,
    const job_id_t& job_id,
    const job_id_t& parent_job_id,
    const worker_id_t& worker_id,
    const bool& sterile)
  : job_type_(job_type),
  job_name_(job_name),
  job_id_(job_id),
  parent_job_id_(parent_job_id),
  worker_id_(worker_id),
  sterile_(sterile) {
    Initialize();
}

void JobProfile::Initialize() {
  assigned_ = false;
  ready_ = false;
  done_ = false;
}

JobProfile::~JobProfile() {
}

JobType JobProfile::job_type() {
  return job_type_;
}

std::string JobProfile::job_name() {
  return job_name_;
}

job_id_t JobProfile::job_id() {
  return job_id_;
}

job_id_t JobProfile::parent_job_id() {
  return parent_job_id_;
}

worker_id_t JobProfile::worker_id() {
  return worker_id_;
}

bool JobProfile::sterile() {
  return sterile_;
}

bool JobProfile::assigned() {
  return assigned_;
}

bool JobProfile::ready() {
  return ready_;
}

bool JobProfile::done() {
  return done_;
}

double JobProfile::assign_time() {
  return assign_time_;
}

double JobProfile::ready_time() {
  return ready_time_;
}

double JobProfile::done_time() {
  return done_time_;
}

double JobProfile::execute_duration() {
  return execute_duration_;
}

void JobProfile::set_job_type(JobType job_type) {
  job_type_ = job_type;
}

void JobProfile::set_job_name(std::string job_name) {
  job_name_ = job_name;
}

void JobProfile::set_job_id(job_id_t job_id) {
  job_id_ = job_id;
}

void JobProfile::set_parent_job_id(job_id_t parent_job_id) {
  parent_job_id_ = parent_job_id;
}

void JobProfile::set_worker_id(worker_id_t worker_id) {
  worker_id_ = worker_id;
}

void JobProfile::set_sterile(bool flag) {
  sterile_ = flag;
}

void JobProfile::set_assigned(bool flag) {
  assigned_ = flag;
}

void JobProfile::set_ready(bool flag) {
  ready_ = flag;
}

void JobProfile::set_done(bool flag) {
  done_ = flag;
}

void JobProfile::set_assign_time(double assign_time) {
  assign_time_ = assign_time;
}

void JobProfile::set_ready_time(double ready_time) {
  ready_time_ = ready_time;
}

void JobProfile::set_done_time(double done_time) {
  done_time_ = done_time;
}

void JobProfile::set_execute_duration(double execute_duration) {
  execute_duration_ = execute_duration;
}

IDSet<job_id_t>* JobProfile::waiting_set_p() {
  return &waiting_set_;
}

void JobProfile::add_log_entry(
    worker_id_t worker_id,
    job_id_t job_id,
    std::string job_name,
    double done_time) {
  dependency_log_.push_back(LogEntry(worker_id, job_id, job_name, done_time));
}


std::string JobProfile::PrintDependencyLog() {
  std::string rval;
  rval += "\n+++++++ Dpendency Log Begin +++++++\n";

  DependencyLog::iterator iter = dependency_log_.begin();
  for (; iter != dependency_log_.end(); ++iter) {
    rval += iter->ToString();
  }

  rval += "\n++++++++ Dpendency Log End ++++++++\n";
  return rval;
}

std::string JobProfile::Print() {
  std::string rval;
  std::ostringstream ss;

  rval += "\n*****************************************\n";
  rval += "*********** Job Profile Begin ***********\n";
  rval += "*****************************************\n";

  rval += "job_id: ";
  ss << job_id_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "job_name: ";
  rval += job_name_;
  rval += "\n";

  rval += "worker_id: ";
  ss << worker_id_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "parent_job_id: ";
  ss << parent_job_id_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "assign_time: ";
  ss << assign_time_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "ready_time: ";
  ss << ready_time_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "done_time: ";
  ss << done_time_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "execute_duration: ";
  ss << execute_duration_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += PrintDependencyLog();


  rval += "\n*****************************************\n";
  rval += "************ Job Profile End ************\n";
  rval += "*****************************************\n";
  return rval;
}






JobProfile::LogEntry::LogEntry(
    worker_id_t worker_id,
    job_id_t job_id,
    std::string job_name,
    double done_time) :
  worker_id_(worker_id),
  job_id_(job_id),
  job_name_(job_name),
  done_time_(done_time) {
}

JobProfile::LogEntry::~LogEntry() {
}

worker_id_t JobProfile::LogEntry::worker_id() {
  return worker_id_;
}

job_id_t JobProfile::LogEntry::job_id() {
  return job_id_;
}

std::string JobProfile::LogEntry::job_name() {
  return job_name_;
}

double JobProfile::LogEntry::done_time() {
  return done_time_;
}

std::string JobProfile::LogEntry::ToString() {
  std::string rval;
  std::ostringstream ss;

  rval += "\n=== Dpendency Entry Begin ===\n";

  rval += "worker_id: ";
  ss << worker_id_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "job_id: ";
  ss << job_id_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "job_name: ";
  rval += job_name_;
  rval += "\n";

  rval += "done_time: ";
  ss << done_time_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "\n==== Dpendency Entry End ====\n";

  return rval;
}



