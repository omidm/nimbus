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
    const IDSet<job_id_t>& before_set,
    const job_id_t& parent_job_id,
    const bool& sterile)
  : job_type_(job_type),
  job_name_(job_name),
  job_id_(job_id),
  before_set_(before_set),
  parent_job_id_(parent_job_id),
  sterile_(sterile) {
    Initialize();
}

void JobProfile::Initialize() {
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

IDSet<job_id_t> JobProfile::before_set() {
  return before_set_;
}

IDSet<job_id_t>* JobProfile::before_set_p() {
  return &before_set_;
}

const IDSet<job_id_t>* JobProfile::before_set_p() const {
  return &before_set_;
}

job_id_t JobProfile::parent_job_id() {
  return parent_job_id_;
}

bool JobProfile::sterile() {
  return sterile_;
}

bool JobProfile::done() {
  return done_;
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

void JobProfile::set_before_set(IDSet<job_id_t> before_set) {
  before_set_ = before_set;
}

void JobProfile::set_parent_job_id(job_id_t parent_job_id) {
  parent_job_id_ = parent_job_id;
}

void JobProfile::set_sterile(bool flag) {
  sterile_ = flag;
}

void JobProfile::set_done(bool flag) {
  done_ = flag;
}



