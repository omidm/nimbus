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
  * Job done command to signal completion of a job.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/job_done_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

JobDoneCommand::JobDoneCommand() {
  name_ = JOB_DONE_NAME;
  type_ = JOB_DONE;
}

JobDoneCommand::JobDoneCommand(const ID<job_id_t>& job_id,
    const IDSet<job_id_t>& after_set,
    const Parameter& params)
: job_id_(job_id), after_set_(after_set), params_(params) {
  name_ = JOB_DONE_NAME;
  type_ = JOB_DONE;
  run_time_ = 0;
  wait_time_ = 0;
}

JobDoneCommand::JobDoneCommand(const ID<job_id_t>& job_id,
    const IDSet<job_id_t>& after_set,
    const Parameter& params,
    const double run_time,
    const double wait_time)
: job_id_(job_id), after_set_(after_set), params_(params),
  run_time_(run_time), wait_time_(wait_time) {
  name_ = JOB_DONE_NAME;
  type_ = JOB_DONE;
}

JobDoneCommand::~JobDoneCommand() {
}

SchedulerCommand* JobDoneCommand::Clone() {
  return new JobDoneCommand();
}

bool JobDoneCommand::Parse(const std::string& params) {
  int num = 5;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: JobDoneCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: JobDone has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  if (!job_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid job id." << std::endl;
    return false;
  }

  iter++;
  if (!after_set_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid after set." << std::endl;
    return false;
  }

  iter++;
  if (!params_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid parameter." << std::endl;
    return false;
  }

  /* Error handling when reading values of timing? */
  iter++;
  run_time_ = boost::lexical_cast<double>(*iter);

  iter++;
  wait_time_ = boost::lexical_cast<double>(*iter);

  return true;
}

std::string JobDoneCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (after_set_.toString() + " ");
  str += (params_.toString() + " ");
  str += (boost::lexical_cast<std::string>(run_time_) + " ");
  str += boost::lexical_cast<std::string>(wait_time_);
  return str;
}

std::string JobDoneCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  str += ("params:" + params_.toString() + " ");
  str += ("run_time: " + boost::lexical_cast<std::string>(run_time_) + " ");
  str += ("wait_time: " + boost::lexical_cast<std::string>(wait_time_));
  return str;
}

ID<job_id_t> JobDoneCommand::job_id() {
  return job_id_;
}

IDSet<job_id_t> JobDoneCommand::after_set() {
  return after_set_;
}

Parameter JobDoneCommand::params() {
  return params_;
}

double JobDoneCommand::run_time() {
  return run_time_;
}

double JobDoneCommand::wait_time() {
  return wait_time_;
}


