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
  * A command sent from the controller to a worker to load a physical data from
  * non-volatile memory. This command is used to rewind from a  distributed
  * checkpoint in the system in case of failure.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "shared/load_data_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

LoadDataCommand::LoadDataCommand() {
  name_ = LOAD_DATA_NAME;
  type_ = LOAD_DATA;
}

LoadDataCommand::LoadDataCommand(const ID<job_id_t>& job_id,
                                 const ID<physical_data_id_t>& to_physical_data_id,
                                 const std::string& handle,
                                 const IDSet<job_id_t>& before)
  : job_id_(job_id),
    to_physical_data_id_(to_physical_data_id),
    handle_(handle),
    before_set_(before) {
  name_ = LOAD_DATA_NAME;
  type_ = LOAD_DATA;
}

LoadDataCommand::~LoadDataCommand() {
}

SchedulerCommand* LoadDataCommand::Clone() {
  return new LoadDataCommand();
}

bool LoadDataCommand::Parse(const std::string& data) {
  LoadDataPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse LoadDataCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool LoadDataCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_load_data()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse LoadDataCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.load_data());
  }
}

std::string LoadDataCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_LOAD_DATA);
  LoadDataPBuf* cmd = buf.mutable_load_data();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string LoadDataCommand::ToString() {
  std::string str;
  str += (name_ + ",");
  str += ("job-id:" + job_id_.ToNetworkData() + ",");
  str += ("to-physical-data-id:" + to_physical_data_id_.ToNetworkData() + ",");
  str += ("handle:" + handle_ + ",");
  str += ("before:" + before_set_.ToNetworkData());
  return str;
}

ID<job_id_t> LoadDataCommand::job_id() {
  return job_id_;
}

ID<physical_data_id_t> LoadDataCommand::to_physical_data_id() {
  return to_physical_data_id_;
}

std::string LoadDataCommand::handle() {
  return handle_;
}

IDSet<job_id_t> LoadDataCommand::before_set() {
  return before_set_;
}

bool LoadDataCommand::ReadFromProtobuf(const LoadDataPBuf& buf) {
  job_id_.set_elem(buf.job_id());
  to_physical_data_id_.set_elem(buf.to_physical_id());
  handle_ = buf.handle();
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  return true;
}

bool LoadDataCommand::WriteToProtobuf(LoadDataPBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_to_physical_id(to_physical_data_id().elem());
  buf->set_handle(handle());
  before_set_.ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  return true;
}

