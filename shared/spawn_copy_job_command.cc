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
  * Spawn copy job command used to send copy jobs to the scheduler from worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/spawn_copy_job_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SpawnCopyJobCommand::SpawnCopyJobCommand() {
  name_ = SPAWN_COPY_JOB_NAME;
  type_ = SPAWN_COPY_JOB;
}

SpawnCopyJobCommand::SpawnCopyJobCommand(const ID<job_id_t>& job_id,
    const ID<logical_data_id_t>& from_logical_id,
    const ID<logical_data_id_t>& to_logical_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    const ID<job_id_t>& parent_job_id,
    const Parameter& params)
: job_id_(job_id),
  from_logical_id_(from_logical_id),
  to_logical_id_(to_logical_id),
  before_set_(before), after_set_(after),
  parent_job_id_(parent_job_id),
  params_(params) {
  name_ = SPAWN_COPY_JOB_NAME;
  type_ = SPAWN_COPY_JOB;
}

SpawnCopyJobCommand::~SpawnCopyJobCommand() {
}

SchedulerCommand* SpawnCopyJobCommand::Clone() {
  return new SpawnCopyJobCommand();
}

bool SpawnCopyJobCommand::Parse(const std::string& params) {
  int num = 7;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: SpawnCopyJobCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: SpawnCopyJobCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  if (!job_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid job id." << std::endl;
    return false;
  }

  iter++;
  if (!from_logical_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid from id." << std::endl;
    return false;
  }

  iter++;
  if (!to_logical_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid to id." << std::endl;
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
  if (!parent_job_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid parent job id." << std::endl;
    return false;
  }

  iter++;
  if (!params_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid parameter." << std::endl;
    return false;
  }

  return true;
}


std::string SpawnCopyJobCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (from_logical_id_.toString() + " ");
  str += (to_logical_id_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += (after_set_.toString() + " ");
  str += (parent_job_id_.toString() + " ");
  str += params_.toString();

  return str;
}

std::string SpawnCopyJobCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("from-logical:" + from_logical_id_.toString() + " ");
  str += ("to-logical:" + to_logical_id_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  str += ("parent-id:" + parent_job_id_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

ID<job_id_t> SpawnCopyJobCommand::job_id() {
  return job_id_;
}

ID<job_id_t> SpawnCopyJobCommand::parent_job_id() {
  return parent_job_id_;
}

ID<logical_data_id_t> SpawnCopyJobCommand::from_logical_id() {
  return from_logical_id_;
}

ID<logical_data_id_t> SpawnCopyJobCommand::to_logical_id() {
  return to_logical_id_;
}

IDSet<job_id_t> SpawnCopyJobCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> SpawnCopyJobCommand::before_set() {
  return before_set_;
}

Parameter SpawnCopyJobCommand::params() {
  return params_;
}



