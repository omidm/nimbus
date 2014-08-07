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
  name_ = SPAWN_COPY_NAME;
  type_ = SPAWN_COPY;
}

SpawnCopyJobCommand::SpawnCopyJobCommand(const ID<job_id_t>& job_id,
                                         const ID<logical_data_id_t>& from_logical_id,
                                         const ID<logical_data_id_t>& to_logical_id,
                                         const IDSet<job_id_t>& before,
                                         const IDSet<job_id_t>& after,
                                         const ID<job_id_t>& parent_job_id)
  : job_id_(job_id),
    from_logical_id_(from_logical_id),
    to_logical_id_(to_logical_id),
    before_set_(before),
    after_set_(after),
    parent_job_id_(parent_job_id) {
  name_ = SPAWN_COPY_NAME;
  type_ = SPAWN_COPY;
}

SpawnCopyJobCommand::~SpawnCopyJobCommand() {
}

SchedulerCommand* SpawnCopyJobCommand::Clone() {
  return new SpawnCopyJobCommand();
}

bool SpawnCopyJobCommand::Parse(const std::string& data) {
  SubmitCopyJobPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    // Throw an error message
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool SpawnCopyJobCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_submit_copy()) {
    return false;
  } else {
    return ReadFromProtobuf(buf.submit_copy());
  }
}


std::string SpawnCopyJobCommand::toString() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_SPAWN_COPY);
  SubmitCopyJobPBuf* cmd = buf.mutable_submit_copy();
  WriteToProtobuf(cmd);
  buf.SerializeToString(&result);

  return result;
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

bool SpawnCopyJobCommand::ReadFromProtobuf(const SubmitCopyJobPBuf& cmd) {
  job_id_.set_elem(cmd.job_id());
  from_logical_id_.set_elem(cmd.from_id());
  to_logical_id_.set_elem(cmd.to_id());
  before_set_.ConvertFromRepeatedField(cmd.before_set().ids());
  after_set_.ConvertFromRepeatedField(cmd.after_set().ids());
  parent_job_id_.set_elem(cmd.parent_id());
  return true;
}

bool SpawnCopyJobCommand::WriteToProtobuf(SubmitCopyJobPBuf* cmd) {
  cmd->set_job_id(job_id().elem());
  cmd->set_from_id(from_logical_id().elem());
  cmd->set_to_id(to_logical_id().elem());
  before_set().ConvertToRepeatedField(cmd->mutable_before_set()->mutable_ids());
  after_set().ConvertToRepeatedField(cmd->mutable_after_set()->mutable_ids());
  cmd->set_parent_id(parent_job_id().elem());
  return true;
}


