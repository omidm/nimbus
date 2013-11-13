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
  * Remote copy receive command to issue receiver side of the copy job to the
  * worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/remote_copy_receive_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;
RemoteCopyReceiveCommand::RemoteCopyReceiveCommand() {
  name_ = REMOTE_COPY_RECEIVE_NAME;
  type_ = REMOTE_COPY_RECEIVE;
}

RemoteCopyReceiveCommand::RemoteCopyReceiveCommand(const ID<job_id_t>& job_id,
    const ID<physical_data_id_t>& to_physical_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after)
: job_id_(job_id),
  to_physical_data_id_(to_physical_data_id),
  before_set_(before), after_set_(after) {
  name_ = REMOTE_COPY_RECEIVE_NAME;
  type_ = REMOTE_COPY_RECEIVE;
}

RemoteCopyReceiveCommand::~RemoteCopyReceiveCommand() {
}

SchedulerCommand* RemoteCopyReceiveCommand::Clone() {
  return new RemoteCopyReceiveCommand();
}

bool RemoteCopyReceiveCommand::Parse(const std::string& params) {
  int num = 4;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: RemoteCopyReceiveCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: RemoteCopyReceiveCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  if (!job_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid job id." << std::endl;
    return false;
  }

  iter++;
  if (!to_physical_data_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid to data id." << std::endl;
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

std::string RemoteCopyReceiveCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (to_physical_data_id_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += after_set_.toString();
  return str;
}

std::string RemoteCopyReceiveCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("job-id:" + job_id_.toString() + " ");
  str += ("to-physical-data-id:" + to_physical_data_id_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString());
  return str;
}

ID<job_id_t> RemoteCopyReceiveCommand::job_id() {
  return job_id_;
}

ID<physical_data_id_t> RemoteCopyReceiveCommand::to_physical_data_id() {
  return to_physical_data_id_;
}

IDSet<job_id_t> RemoteCopyReceiveCommand::after_set() {
  return after_set_;
}

IDSet<job_id_t> RemoteCopyReceiveCommand::before_set() {
  return before_set_;
}


