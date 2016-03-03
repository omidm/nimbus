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
  * The controller sends combine jobs to workers to combine the scratch
  * instances of a logical id for reduction.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/combine_job_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

CombineJobCommand::CombineJobCommand() {
  name_ = EXECUTE_COMBINE_NAME;
  type_ = EXECUTE_COMBINE;
}

CombineJobCommand::CombineJobCommand(const std::string& job_name,
                                     const ID<job_id_t>& job_id,
                                     const IDSet<physical_data_id_t>& scratch,
                                     const IDSet<physical_data_id_t>& reduce,
                                     const IDSet<job_id_t>& before,
                                     const IDSet<job_id_t>& extra_dependency,
                                     const GeometricRegion& region)
  : job_name_(job_name),
    job_id_(job_id),
    scratch_set_(scratch),
    reduce_set_(reduce),
    before_set_(before),
    extra_dependency_(extra_dependency),
    region_(region) {
  name_ = EXECUTE_COMBINE_NAME;
  type_ = EXECUTE_COMBINE;
}

CombineJobCommand::~CombineJobCommand() {
}

SchedulerCommand* CombineJobCommand::Clone() {
  return new CombineJobCommand();
}

bool CombineJobCommand::Parse(const std::string& data) {
  ExecuteCombineJobPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse CombineJobCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}


bool CombineJobCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_execute_combine()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse CombineJobCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.execute_combine());
  }
}


std::string CombineJobCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn combine field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_EXECUTE_COMBINE);
  ExecuteCombineJobPBuf* cmd = buf.mutable_execute_combine();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string CombineJobCommand::ToString() {
  std::string str;
  str += (name_ + ",");
  str += ("name:" + job_name_ + ",");
  str += ("id:" + job_id_.ToNetworkData() + ",");
  str += ("scratch:" + scratch_set_.ToNetworkData() + ",");
  str += ("reduce:" + reduce_set_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData() + ",");
  str += ("extra-dependency:" + extra_dependency_.ToNetworkData() + ",");
  str += ("region:" + region_.ToNetworkData() + ",");
  return str;
}

std::string CombineJobCommand::job_name() {
  return job_name_;
}

ID<job_id_t> CombineJobCommand::job_id() {
  return job_id_;
}

IDSet<physical_data_id_t> CombineJobCommand::scratch_set() {
  return scratch_set_;
}

IDSet<physical_data_id_t>* CombineJobCommand::scratch_set_p() {
  return &scratch_set_;
}

IDSet<physical_data_id_t> CombineJobCommand::reduce_set() {
  return reduce_set_;
}

IDSet<physical_data_id_t>* CombineJobCommand::reduce_set_p() {
  return &reduce_set_;
}

IDSet<job_id_t> CombineJobCommand::before_set() {
  return before_set_;
}

IDSet<job_id_t> CombineJobCommand::extra_dependency() {
  return extra_dependency_;
}

IDSet<job_id_t>* CombineJobCommand::before_set_p() {
  return &before_set_;
}

IDSet<job_id_t>* CombineJobCommand::extra_dependency_p() {
  return &extra_dependency_;
}

GeometricRegion CombineJobCommand::region() {
  return region_;
}

bool CombineJobCommand::ReadFromProtobuf(const ExecuteCombineJobPBuf& buf) {
  job_name_ = buf.name();
  job_id_.set_elem(buf.job_id());
  scratch_set_.ConvertFromRepeatedField(buf.scratch_set().ids());
  reduce_set_.ConvertFromRepeatedField(buf.reduce_set().ids());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  extra_dependency_.ConvertFromRepeatedField(buf.extra_dependency().ids());
  region_.FillInValues(&buf.region());
  return true;
}

bool CombineJobCommand::WriteToProtobuf(ExecuteCombineJobPBuf* buf) {
  buf->set_name(job_name());
  buf->set_job_id(job_id().elem());
  scratch_set().ConvertToRepeatedField(buf->mutable_scratch_set()->mutable_ids());
  reduce_set().ConvertToRepeatedField(buf->mutable_reduce_set()->mutable_ids());
  before_set().ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  extra_dependency().ConvertToRepeatedField(buf->mutable_extra_dependency()->mutable_ids());
  region_.FillInMessage(buf->mutable_region());
  return true;
}
