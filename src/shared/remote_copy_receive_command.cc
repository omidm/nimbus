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
  * A remote copy operation between two workers has two jobs: the
  * send and receive. This is the receive half.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "src/shared/remote_copy_receive_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

RemoteCopyReceiveCommand::RemoteCopyReceiveCommand() {
  name_ = REMOTE_RECEIVE_NAME;
  type_ = REMOTE_RECEIVE;
}

RemoteCopyReceiveCommand::RemoteCopyReceiveCommand(const ID<job_id_t>& job_id,
                                                   const ID<physical_data_id_t>& to_physical_data_id, // NOLINT
                                                   const IDSet<job_id_t>& before)
: job_id_(job_id),
  to_physical_data_id_(to_physical_data_id),
  before_set_(before) {
  name_ = REMOTE_RECEIVE_NAME;
  type_ = REMOTE_RECEIVE;
}

RemoteCopyReceiveCommand::~RemoteCopyReceiveCommand() {
}

SchedulerCommand* RemoteCopyReceiveCommand::Clone() {
  return new RemoteCopyReceiveCommand();
}

bool RemoteCopyReceiveCommand::Parse(const std::string& data) {
  RemoteCopyReceivePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RemoteCopyReceiveCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool RemoteCopyReceiveCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_remote_receive()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RemoteCopyReceiveCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.remote_receive());
  }
}

std::string RemoteCopyReceiveCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_REMOTE_RECEIVE);
  RemoteCopyReceivePBuf* cmd = buf.mutable_remote_receive();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string RemoteCopyReceiveCommand::ToString() {
  std::string str;
  str += (name_ + ",");
  str += ("job-id:" + job_id_.ToNetworkData() + ",");
  str += ("to-physical-data-id:" + to_physical_data_id_.ToNetworkData() + ",");
  str += ("before:" + before_set_.ToNetworkData());
  return str;
}

ID<job_id_t> RemoteCopyReceiveCommand::job_id() {
  return job_id_;
}

ID<physical_data_id_t> RemoteCopyReceiveCommand::to_physical_data_id() {
  return to_physical_data_id_;
}

IDSet<job_id_t> RemoteCopyReceiveCommand::before_set() {
  return before_set_;
}

IDSet<job_id_t>* RemoteCopyReceiveCommand::before_set_p() {
  return &before_set_;
}

bool RemoteCopyReceiveCommand::ReadFromProtobuf(const RemoteCopyReceivePBuf& buf) {
  job_id_.set_elem(buf.job_id());
  to_physical_data_id_.set_elem(buf.physical_id());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  return true;
}

bool RemoteCopyReceiveCommand::WriteToProtobuf(RemoteCopyReceivePBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_physical_id(to_physical_data_id().elem());
  before_set_.ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  return true;
}

