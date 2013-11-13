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
  * Create data command to create physical instances of logical data on the
  * worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/create_data_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

CreateDataCommand::CreateDataCommand() {
  name_ = CREATE_DATA_NAME;
  type_ = CREATE_DATA;
}

CreateDataCommand::CreateDataCommand(const ID<job_id_t>& job_id,
    const std::string& data_name,
    const ID<logical_data_id_t>& logical_data_id,
    const ID<physical_data_id_t>& physical_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after)
: job_id_(job_id),
  data_name_(data_name),
  logical_data_id_(logical_data_id),
  physical_data_id_(physical_data_id),
  before_set_(before), after_set_(after) {
  name_ = CREATE_DATA_NAME;
  type_ = CREATE_DATA;
}

CreateDataCommand::~CreateDataCommand() {
}

SchedulerCommand* CreateDataCommand::Clone() {
  return new CreateDataCommand();
}

bool CreateDataCommand::Parse(const std::string& params) {
  int num = 6;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: CreateDataCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: CreateDataCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  if (!job_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid job id." << std::endl;
    return false;
  }

  iter++;
  data_name_ = *iter;

  iter++;
  if (!logical_data_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid logical data id." << std::endl;
    return false;
  }

  iter++;
  if (!physical_data_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid physical data id." << std::endl;
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

  return true;
}

std::string CreateDataCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (data_name_ + " ");
  str += (logical_data_id_.toString() + " ");
  str += (physical_data_id_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += after_set_.toString();
  return str;
}

std::string CreateDataCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("job-id:" + job_id_.toString() + " ");
  str += ("name:" + data_name_ + " ");
  str += ("logical-data-id:" + logical_data_id_.toString() + " ");
  str += ("physical-data-id:" + physical_data_id_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString());
  return str;
}

ID<job_id_t> CreateDataCommand::job_id() {
  return job_id_;
}

std::string CreateDataCommand::data_name() {
  return data_name_;
}

ID<logical_data_id_t> CreateDataCommand::logical_data_id() {
  return logical_data_id_;
}

ID<physical_data_id_t> CreateDataCommand::physical_data_id() {
  return physical_data_id_;
}

IDSet<job_id_t> CreateDataCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> CreateDataCommand::before_set() {
  return before_set_;
}


