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
  * Compute job command used to send compute jobs from scheduler to workers.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/compute_job_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

ComputeJobCommand::ComputeJobCommand() {
  name_ = EXECUTE_COMPUTE_NAME;
  type_ = EXECUTE_COMPUTE;
}

ComputeJobCommand::ComputeJobCommand(const std::string& job_name,
                                     const ID<job_id_t>& job_id,
                                     const IDSet<physical_data_id_t>& read,
                                     const IDSet<physical_data_id_t>& write, // NOLINT
                                     const IDSet<job_id_t>& before,
                                     const IDSet<job_id_t>& after,
                                     const ID<job_id_t>& future_job_id,
                                     const bool& sterile,
                                     const GeometricRegion& region,
                                     const Parameter& params)
  : job_name_(job_name),
    job_id_(job_id),
    read_set_(read),
    write_set_(write),
    before_set_(before),
    after_set_(after),
    future_job_id_(future_job_id),
    sterile_(sterile),
    region_(region),
    params_(params) {
  name_ = EXECUTE_COMPUTE_NAME;
  type_ = EXECUTE_COMPUTE;
}

ComputeJobCommand::~ComputeJobCommand() {
}

SchedulerCommand* ComputeJobCommand::Clone() {
  return new ComputeJobCommand();
}

bool ComputeJobCommand::Parse(const std::string& data) {
  ExecuteComputeJobPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse ComputeJobCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}


bool ComputeJobCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_execute_compute()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse ComputeJobCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.execute_compute());
  }
}


std::string ComputeJobCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_EXECUTE_COMPUTE);
  ExecuteComputeJobPBuf* cmd = buf.mutable_execute_compute();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string ComputeJobCommand::ToString() {
  std::string str;
  str += (name_ + ",");
  str += ("name:" + job_name_ + ",");
  str += ("id:" + job_id_.ToNetworkData() + ",");
  str += ("read:" + read_set_.ToNetworkData() + ",");
  str += ("write:" + write_set_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData() + ",");
  str += ("after:" + after_set_.ToNetworkData() + ",");
  str += ("future-job-id:" + future_job_id_.ToNetworkData() + ",");
  str += ("params:" + params_.ToNetworkData() + ",");
  str += ("region:" + region_.ToNetworkData() + ",");
  if (sterile_) {
    str += "sterile";
  } else {
    str += "not_sterile";
  }
  return str;
}

std::string ComputeJobCommand::job_name() {
  return job_name_;
}

ID<job_id_t> ComputeJobCommand::job_id() {
  return job_id_;
}

IDSet<physical_data_id_t> ComputeJobCommand::read_set() {
  return read_set_;
}

IDSet<physical_data_id_t>* ComputeJobCommand::read_set_p() {
  return &read_set_;
}

IDSet<physical_data_id_t> ComputeJobCommand::write_set() {
  return write_set_;
}

IDSet<physical_data_id_t>* ComputeJobCommand::write_set_p() {
  return &write_set_;
}

IDSet<job_id_t> ComputeJobCommand::after_set() {
  return after_set_;
}

IDSet<job_id_t>* ComputeJobCommand::after_set_p() {
  return &after_set_;
}

IDSet<job_id_t> ComputeJobCommand::before_set() {
  return before_set_;
}

IDSet<job_id_t>* ComputeJobCommand::before_set_p() {
  return &before_set_;
}

Parameter ComputeJobCommand::params() {
  return params_;
}

bool ComputeJobCommand::sterile() {
  return sterile_;
}

GeometricRegion ComputeJobCommand::region() {
  return region_;
}

ID<job_id_t> ComputeJobCommand::future_job_id() {
  return future_job_id_;
}

bool ComputeJobCommand::ReadFromProtobuf(const ExecuteComputeJobPBuf& buf) {
  job_name_ = buf.name();
  job_id_.set_elem(buf.job_id());
  read_set_.ConvertFromRepeatedField(buf.read_set().ids());
  write_set_.ConvertFromRepeatedField(buf.write_set().ids());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  after_set_.ConvertFromRepeatedField(buf.after_set().ids());
  future_job_id_.set_elem(buf.future_job_id());
  sterile_ = buf.sterile();
  region_.FillInValues(&buf.region());
  SerializedData d(buf.params());
  params_.set_ser_data(d);
  return true;
}

bool ComputeJobCommand::WriteToProtobuf(ExecuteComputeJobPBuf* buf) {
  buf->set_name(job_name());
  buf->set_job_id(job_id().elem());
  read_set().ConvertToRepeatedField(buf->mutable_read_set()->mutable_ids());
  write_set().ConvertToRepeatedField(buf->mutable_write_set()->mutable_ids());
  before_set().ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  after_set().ConvertToRepeatedField(buf->mutable_after_set()->mutable_ids());
  buf->set_future_job_id(future_job_id().elem());
  buf->set_sterile(sterile());
  region_.FillInMessage(buf->mutable_region());
  buf->set_params(params().ser_data().ToNetworkData());
  return true;
}
