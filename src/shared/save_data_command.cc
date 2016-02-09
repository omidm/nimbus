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
  * A command sent from the controller to a worker to save a physical data to
  * non-volatile memory. This command is used to create distributed checkpoint
  * in the system to rewind back to in case of failure.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "src/shared/save_data_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SaveDataCommand::SaveDataCommand() {
  name_ = SAVE_DATA_NAME;
  type_ = SAVE_DATA;
}

SaveDataCommand::SaveDataCommand(const ID<job_id_t>& job_id,
                                   const ID<physical_data_id_t>& from_physical_data_id,
                                   const ID<checkpoint_id_t>& checkpoint_id,
                                   const IDSet<job_id_t>& before)
  : job_id_(job_id),
    from_physical_data_id_(from_physical_data_id),
    checkpoint_id_(checkpoint_id),
    before_set_(before) {
  name_ = SAVE_DATA_NAME;
  type_ = SAVE_DATA;
}

SaveDataCommand::~SaveDataCommand() {
}

SchedulerCommand* SaveDataCommand::Clone() {
  return new SaveDataCommand();
}

bool SaveDataCommand::Parse(const std::string& data) {
  SaveDataPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SaveDataCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool SaveDataCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_save_data()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SaveDataCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.save_data());
  }
}

std::string SaveDataCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_SAVE_DATA);
  SaveDataPBuf* cmd = buf.mutable_save_data();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string SaveDataCommand::ToString() {
  std::string str;
  str += (name_ + ",");
  str += ("job-id:" + job_id_.ToNetworkData() + ",");
  str += ("from-physical-data-id:" + from_physical_data_id_.ToNetworkData() + ",");
  str += ("checkpoint-id:" + checkpoint_id_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData());
  return str;
}

ID<job_id_t> SaveDataCommand::job_id() {
  return job_id_;
}

ID<physical_data_id_t> SaveDataCommand::from_physical_data_id() {
  return from_physical_data_id_;
}

ID<checkpoint_id_t> SaveDataCommand::checkpoint_id() {
  return checkpoint_id_;
}

IDSet<job_id_t> SaveDataCommand::before_set() {
  return before_set_;
}

bool SaveDataCommand::ReadFromProtobuf(const SaveDataPBuf& buf) {
  job_id_.set_elem(buf.job_id());
  from_physical_data_id_.set_elem(buf.from_physical_id());
  checkpoint_id_.set_elem(buf.checkpoint_id());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  return true;
}

bool SaveDataCommand::WriteToProtobuf(SaveDataPBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_from_physical_id(from_physical_data_id().elem());
  buf->set_checkpoint_id(checkpoint_id().elem());
  before_set_.ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  return true;
}
