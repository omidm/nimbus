/*
 * Copyright 2014 Stanford University.
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
  * A AddCopyJobCommand is a message sent from a worker to the controller to
  * add a copy job to a job graph template and gradually build it up by
  * individual vertices of the graph.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/add_copy_job_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

AddCopyJobCommand::AddCopyJobCommand() {
  name_ = ADD_COPY_NAME;
  type_ = ADD_COPY;
}

AddCopyJobCommand::AddCopyJobCommand(const ID<job_id_t>& job_id,
                                     const ID<logical_data_id_t>& from_logical_id,
                                     const ID<logical_data_id_t>& to_logical_id,
                                     const IDSet<job_id_t>& before,
                                     const IDSet<job_id_t>& after,
                                     const std::string& job_graph_name)
  : job_id_(job_id),
    from_logical_id_(from_logical_id),
    to_logical_id_(to_logical_id),
    before_set_(before),
    after_set_(after),
    job_graph_name_(job_graph_name) {
  name_ = ADD_COPY_NAME;
  type_ = ADD_COPY;
}

AddCopyJobCommand::~AddCopyJobCommand() {
}

SchedulerCommand* AddCopyJobCommand::Clone() {
  return new AddCopyJobCommand();
}


bool AddCopyJobCommand::Parse(const std::string& data) {
  AddCopyJobPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse AddCopyJobCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool AddCopyJobCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_add_copy()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse AddCopyJobCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.add_copy());
  }
}

std::string AddCopyJobCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the add copy field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_ADD_COPY);
  AddCopyJobPBuf* cmd = buf.mutable_add_copy();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string AddCopyJobCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + job_id_.ToNetworkData() + ",");
  str += ("from-logical:" + from_logical_id_.ToNetworkData() + ",");
  str += ("to-logical:" + to_logical_id_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData() + ",");
  str += ("after:" + after_set_.ToNetworkData() + ",");
  str += ("job-graph-name:" + job_graph_name_ + ",");
  return str;
}

ID<job_id_t> AddCopyJobCommand::job_id() {
  return job_id_;
}

ID<logical_data_id_t> AddCopyJobCommand::from_logical_id() {
  return from_logical_id_;
}

ID<logical_data_id_t> AddCopyJobCommand::to_logical_id() {
  return to_logical_id_;
}

IDSet<job_id_t> AddCopyJobCommand::after_set() {
  return after_set_;
}

IDSet<job_id_t> AddCopyJobCommand::before_set() {
  return before_set_;
}

std::string AddCopyJobCommand::job_graph_name() {
  return job_graph_name_;
}

bool AddCopyJobCommand::ReadFromProtobuf(const AddCopyJobPBuf& buf) {
  job_id_.set_elem(buf.job_id());
  from_logical_id_.set_elem(buf.from_logical_id());
  to_logical_id_.set_elem(buf.to_logical_id());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  after_set_.ConvertFromRepeatedField(buf.after_set().ids());
  job_graph_name_ = buf.job_graph_name();

  return true;
}

bool AddCopyJobCommand::WriteToProtobuf(AddCopyJobPBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_from_logical_id(from_logical_id().elem());
  buf->set_to_logical_id(to_logical_id().elem());
  before_set().ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  after_set().ConvertToRepeatedField(buf->mutable_after_set()->mutable_ids());
  buf->set_job_graph_name(job_graph_name());
  return true;
}
