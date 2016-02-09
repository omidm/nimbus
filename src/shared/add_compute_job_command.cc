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
  * A AddComputeJobCommand is a message sent from a worker to the controller to
  * add a compute job to a job graph template and gradually build it up by
  * individual vertices of the graph.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/add_compute_job_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

AddComputeJobCommand::AddComputeJobCommand() {
  name_ = ADD_COMPUTE_NAME;
  type_ = ADD_COMPUTE;
}

AddComputeJobCommand::AddComputeJobCommand(const std::string& job_name,
                                           const ID<job_id_t>& job_id,
                                           const IDSet<logical_data_id_t>& read,
                                           const IDSet<logical_data_id_t>& write,
                                           const IDSet<job_id_t>& before,
                                           const IDSet<job_id_t>& after,
                                           const bool& sterile,
                                           const GeometricRegion& region,
                                           const ID<job_id_t>& future_job_id,
                                           const std::string& job_graph_name)
  : job_name_(job_name),
    job_id_(job_id),
    read_set_(read),
    write_set_(write),
    before_set_(before),
    after_set_(after),
    sterile_(sterile),
    region_(region),
    future_job_id_(future_job_id),
    job_graph_name_(job_graph_name) {
  name_ = ADD_COMPUTE_NAME;
  type_ = ADD_COMPUTE;
}

AddComputeJobCommand::~AddComputeJobCommand() {
}

SchedulerCommand* AddComputeJobCommand::Clone() {
  return new AddComputeJobCommand();
}


bool AddComputeJobCommand::Parse(const std::string& data) {
  AddComputeJobPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse AddComputeJobCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool AddComputeJobCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_add_compute()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse AddComputeJobCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.add_compute());
  }
}

std::string AddComputeJobCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the add compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_ADD_COMPUTE);
  AddComputeJobPBuf* cmd = buf.mutable_add_compute();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string AddComputeJobCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_name_ + ",");
  str += ("id:" + job_id_.ToNetworkData() + ",");
  str += ("read:" + read_set_.ToNetworkData() + ",");
  str += ("write:" + write_set_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData() + ",");
  str += ("after:" + after_set_.ToNetworkData() + ",");
  str += ("future-id:" + future_job_id_.ToNetworkData() + ",");
  str += ("region:" + region_.ToNetworkData() + ",");
  str += ("job-graph-name:" + job_graph_name_ + ",");
  if (sterile_) {
    str += "sterile";
  } else {
    str += "not_sterile";
  }
  return str;
}

std::string AddComputeJobCommand::job_name() {
  return job_name_;
}

ID<job_id_t> AddComputeJobCommand::job_id() {
  return job_id_;
}

ID<job_id_t> AddComputeJobCommand::future_job_id() {
  return future_job_id_;
}

IDSet<logical_data_id_t> AddComputeJobCommand::read_set() {
  return read_set_;
}

IDSet<logical_data_id_t> AddComputeJobCommand::write_set() {
  return write_set_;
}

IDSet<job_id_t> AddComputeJobCommand::after_set() {
  return after_set_;
}

IDSet<job_id_t> AddComputeJobCommand::before_set() {
  return before_set_;
}

std::string AddComputeJobCommand::job_graph_name() {
  return job_graph_name_;
}

bool AddComputeJobCommand::sterile() {
  return sterile_;
}

GeometricRegion AddComputeJobCommand::region() {
  return region_;
}


bool AddComputeJobCommand::ReadFromProtobuf(const AddComputeJobPBuf& buf) {
  job_name_ = buf.job_name();
  job_id_.set_elem(buf.job_id());
  read_set_.ConvertFromRepeatedField(buf.read_set().ids());
  write_set_.ConvertFromRepeatedField(buf.write_set().ids());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  after_set_.ConvertFromRepeatedField(buf.after_set().ids());
  sterile_ = buf.sterile();
  region_.FillInValues(&buf.region());
  future_job_id_.set_elem(buf.future_id());
  job_graph_name_ = buf.job_graph_name();

  return true;
}

bool AddComputeJobCommand::WriteToProtobuf(AddComputeJobPBuf* buf) {
  buf->set_job_name(job_name());
  buf->set_job_id(job_id().elem());
  read_set().ConvertToRepeatedField(buf->mutable_read_set()->mutable_ids());
  write_set().ConvertToRepeatedField(buf->mutable_write_set()->mutable_ids());
  before_set().ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  after_set().ConvertToRepeatedField(buf->mutable_after_set()->mutable_ids());
  buf->set_sterile(sterile_);
  region_.FillInMessage(buf->mutable_region());
  buf->set_future_id(future_job_id().elem());
  buf->set_job_graph_name(job_graph_name());
  return true;
}
