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
  * A local copy is between two physical objects on a worker: a common
  * use is when data for timestep t is computed with data from
  * timestep t, or when one job needs to read data from timestep t
  * while another writes it for timestep t+1. Copies between workers
  * are different and called remote copies.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "src/shared/local_copy_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

LocalCopyCommand::LocalCopyCommand() {
  name_ = LOCAL_COPY_NAME;
  type_ = LOCAL_COPY;
}

LocalCopyCommand::LocalCopyCommand(const ID<job_id_t>& job_id,
                                   const ID<physical_data_id_t>& from_physical_data_id,
                                   const ID<physical_data_id_t>& to_physical_data_id,
                                   const IDSet<job_id_t>& before,
                                   const IDSet<job_id_t>& extra_dependency)
  : job_id_(job_id),
    from_physical_data_id_(from_physical_data_id),
    to_physical_data_id_(to_physical_data_id),
    before_set_(before),
    extra_dependency_(extra_dependency) {
  name_ = LOCAL_COPY_NAME;
  type_ = LOCAL_COPY;
}

LocalCopyCommand::~LocalCopyCommand() {
}

SchedulerCommand* LocalCopyCommand::Clone() {
  return new LocalCopyCommand();
}

bool LocalCopyCommand::Parse(const std::string& data) {
  LocalCopyPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse LocalCopyCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool LocalCopyCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_local_copy()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse LocalCopyCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.local_copy());
  }
}

std::string LocalCopyCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_LOCAL_COPY);
  LocalCopyPBuf* cmd = buf.mutable_local_copy();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string LocalCopyCommand::ToString() {
  std::string str;
  str += (name_ + ",");
  str += ("job-id:" + job_id_.ToNetworkData() + ",");
  str += ("from-physical-data-id:" + from_physical_data_id_.ToNetworkData() + ",");
  str += ("to-physical-data-id:" + to_physical_data_id_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData());
  str += ("extra_dependency:" + extra_dependency_.ToNetworkData());
  return str;
}

ID<job_id_t> LocalCopyCommand::job_id() {
  return job_id_;
}

ID<physical_data_id_t> LocalCopyCommand::from_physical_data_id() {
  return from_physical_data_id_;
}

ID<physical_data_id_t> LocalCopyCommand::to_physical_data_id() {
  return to_physical_data_id_;
}

IDSet<job_id_t> LocalCopyCommand::before_set() {
  return before_set_;
}

IDSet<job_id_t>* LocalCopyCommand::before_set_p() {
  return &before_set_;
}

IDSet<job_id_t> LocalCopyCommand::extra_dependency() {
  return extra_dependency_;
}

IDSet<job_id_t>* LocalCopyCommand::extra_dependency_p() {
  return &extra_dependency_;
}

bool LocalCopyCommand::ReadFromProtobuf(const LocalCopyPBuf& buf) {
  job_id_.set_elem(buf.job_id());
  from_physical_data_id_.set_elem(buf.from_physical_id());
  to_physical_data_id_.set_elem(buf.to_physical_id());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  extra_dependency_.ConvertFromRepeatedField(buf.extra_dependency().ids());
  return true;
}

bool LocalCopyCommand::WriteToProtobuf(LocalCopyPBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_from_physical_id(from_physical_data_id().elem());
  buf->set_to_physical_id(to_physical_data_id().elem());
  before_set_.ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  extra_dependency_.ConvertToRepeatedField(buf->mutable_extra_dependency()->mutable_ids());
  return true;
}
