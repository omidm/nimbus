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

#include "src/shared/create_data_command.h"

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
                                     const IDSet<job_id_t>& before)
  : job_id_(job_id),
    data_name_(data_name),
    logical_data_id_(logical_data_id),
    physical_data_id_(physical_data_id),
    before_set_(before) {
  name_ = CREATE_DATA_NAME;
  type_ = CREATE_DATA;
}

CreateDataCommand::~CreateDataCommand() {
}

SchedulerCommand* CreateDataCommand::Clone() {
  return new CreateDataCommand();
}

bool CreateDataCommand::Parse(const std::string& data) {
  CreateDataPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse CreateDataCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool CreateDataCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_create_data()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse CreateDataCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.create_data());
  }
}

std::string CreateDataCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_CREATE_DATA);
  CreateDataPBuf* cmd = buf.mutable_create_data();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string CreateDataCommand::ToString() {
  std::string str;
  str += (name_ + ",");
  str += ("job-id:" + job_id_.ToNetworkData() + ",");
  str += ("name:" + data_name_ + ",");
  str += ("logical-data-id:" + logical_data_id_.ToNetworkData() + ",");
  str += ("physical-data-id:" + physical_data_id_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData());
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

IDSet<job_id_t> CreateDataCommand::before_set() {
  return before_set_;
}

bool CreateDataCommand::ReadFromProtobuf(const CreateDataPBuf& buf) {
  job_id_.set_elem(buf.job_id());
  data_name_ = buf.name();
  logical_data_id_.set_elem(buf.logical_id());
  physical_data_id_.set_elem(buf.physical_id());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  return true;
}

bool CreateDataCommand::WriteToProtobuf(CreateDataPBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_name(data_name());
  buf->set_logical_id(logical_data_id().elem());
  buf->set_physical_id(physical_data_id().elem());
  before_set().ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  return true;
}

