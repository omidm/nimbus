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
  * Spawn job job command used to spawn jobs from worker to scheduler.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/spawn_job_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SpawnJobCommand::SpawnJobCommand() {
  name_ = SPAWN_JOB_NAME;
  type_ = SPAWN_JOB;
}

SpawnJobCommand::SpawnJobCommand(const std::string& job_name,
    const IDSet<job_id_t>& job_id,
    const IDSet<logical_data_id_t>& read, const IDSet<logical_data_id_t>& write,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    const JobType& job_type, const Parameter& params)
: job_name_(job_name), job_id_(job_id),
  read_set_(read), write_set_(write),
  before_set_(before), after_set_(after),
  job_type_(job_type), params_(params) {
  name_ = SPAWN_JOB_NAME;
  type_ = SPAWN_JOB;
}

SpawnJobCommand::~SpawnJobCommand() {
}

SchedulerCommand* SpawnJobCommand::Clone() {
  return new SpawnJobCommand();
}

bool SpawnJobCommand::Parse(const std::string& params) {
  int num = 8;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: SpawnJobCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: SpawnJobCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  job_name_ = *iter;

  iter++;
  if (job_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid job id set." << std::endl;
    return false;
  }

  iter++;
  if (!read_set_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid read set." << std::endl;
    return false;
  }

  iter++;
  if (!write_set_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid write set." << std::endl;
    return false;
  }

  iter++;
  if (!before_set_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid before set." << std::endl;
    return false;
  }

  iter++;
  if (!after_set_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid after set." << std::endl;
    return false;
  }

  iter++;
  if (*iter == "COMP") {
    job_type_ = JOB_COMP;
  } else if (*iter == "COPY") {
    job_type_ = JOB_COPY;
  } else {
    std::cout << "ERROR: Unknown job type." << std::endl;
    return false;
  }

  iter++;
  if (params_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid parameter." << std::endl;
    return false;
  }

  return true;
}

std::string SpawnJobCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_name_ + " ");
  str += (job_id_.toString() + " ");
  str += (read_set_.toString() + " ");
  str += (write_set_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += (after_set_.toString() + " ");
  if (job_type_ == JOB_COMP)
    str += "COMP ";
  else
    str += "COPY ";
  str += params_.toString();

  return str;
}

std::string SpawnJobCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("read:" + read_set_.toString() + " ");
  str += ("write:" + write_set_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  if (job_type_ == JOB_COMP)
    str += "type:COMP ";
  else
    str += "type:COPY ";
  str += ("params:" + params_.toString());
  return str;
}

std::string SpawnJobCommand::job_name() {
  return job_name_;
}

JobType SpawnJobCommand::job_type() {
  return job_type_;
}

IDSet<job_id_t> SpawnJobCommand::job_id() {
  return job_id_;
}

IDSet<logical_data_id_t> SpawnJobCommand::read_set() {
  return read_set_;
}

IDSet<logical_data_id_t> SpawnJobCommand::write_set() {
  return write_set_;
}

IDSet<job_id_t> SpawnJobCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> SpawnJobCommand::before_set() {
  return before_set_;
}

Parameter SpawnJobCommand::params() {
  return params_;
}


