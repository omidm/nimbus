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
#include "shared/protobuf_compiled/commands.pb.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SpawnJobCommand::SpawnJobCommand() {
  name_ = SPAWN_JOB_NAME;
  type_ = SPAWN_JOB;
}

SpawnJobCommand::SpawnJobCommand(const std::string& job_name,
                                 const ID<job_id_t>& job_id,
                                 const IDSet<logical_data_id_t>& read,
                                 const IDSet<logical_data_id_t>& write,
                                 const IDSet<job_id_t>& before,
                                 const IDSet<job_id_t>& after,
                                 const ID<job_id_t>& future_id,
                                 const JobType& job_type,
                                 const Parameter& params)
  : job_name_(job_name),
    job_id_(job_id),
    read_set_(read),
    write_set_(write),
    before_set_(before),
    after_set_(after),
    future_id_(future_id),
    job_type_(job_type),
    params_(params) {
  name_ = SPAWN_JOB_NAME;
  type_ = SPAWN_JOB;
}

SpawnJobCommand::~SpawnJobCommand() {
}

SchedulerCommand* SpawnJobCommand::Clone() {
  return new SpawnJobCommand();
}

bool SpawnJobCommand::Parse(const std::string& data) {
  SubmitJobCommand cmd;
  bool result = cmd.ParseFromString(data);

  if (!result) {
    std::cout << "ERROR: Failed to parse SubmitJobCommand." << std::endl;
    return false;
  }

  ReadFromMsg(&cmd);
  return true;
}

std::string SpawnJobCommand::toString() {
  std::string str;
  SubmitJobCommand cmd;
  WriteToMsg(&cmd);
  cmd.SerializeToString(&str);

  // Note asymmetry of parsing forces this copy.
  // That is, toString inserts the command name,
  // while Parse assumes the command name has been consumed.
  // Correspondingly, the protocol buffer doesn't have
  // the command name because its typing needs to follow
  // the higher-level parsers expectations of ASCII, but
  // we need to insert it in toString. Asymmetry is
  // a bad idea!
  std::string result = name_ + " ";
  result += str;
  return result;
}

void SpawnJobCommand::WriteToMsg(SubmitJobCommand* cmd) {
  cmd->set_name(job_name());
  cmd->set_job_id(job_id().elem());
  read_set().ConvertToRepeatedField(cmd->mutable_read_set()->mutable_ids());
  write_set().ConvertToRepeatedField(cmd->mutable_write_set()->mutable_ids());
  before_set().ConvertToRepeatedField(cmd->mutable_before_set()->mutable_ids());
  after_set().ConvertToRepeatedField(cmd->mutable_after_set()->mutable_ids());
  cmd->set_type((uint32_t)job_type());
  cmd->set_future_id(future_id().elem());
  cmd->set_params(params().ser_data().toString());
}

void SpawnJobCommand::ReadFromMsg(SubmitJobCommand* cmd) {
  job_name_ = cmd->name();
  job_id_.set_elem(cmd->job_id());
  read_set().ConvertFromRepeatedField(cmd->read_set().ids());
  write_set().ConvertFromRepeatedField(cmd->write_set().ids());
  before_set().ConvertFromRepeatedField(cmd->before_set().ids());
  after_set().ConvertFromRepeatedField(cmd->after_set().ids());
  job_type_ = (JobType)cmd->type();
  future_id_.set_elem(cmd->future_id());
  // Is this safe?
  SerializedData d(cmd->params());
  params_.set_ser_data(d);
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

ID<job_id_t> SpawnJobCommand::job_id() {
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


