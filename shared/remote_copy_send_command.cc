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
  * send and receive. This is the send half.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "shared/remote_copy_send_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

RemoteCopySendCommand::RemoteCopySendCommand() {
  name_ = REMOTE_COPY_SEND_NAME;
  type_ = REMOTE_COPY_SEND;
}

RemoteCopySendCommand::RemoteCopySendCommand(const ID<job_id_t>& job_id,
                                             const ID<job_id_t>& receive_job_id,
                                             const ID<physical_data_id_t>& from_physical_data_id,
                                             const ID<worker_id_t>& to_worker_id,
                                             const std::string to_ip, const ID<port_t>& to_port,
                                             const IDSet<job_id_t>& before)
  : job_id_(job_id),
    receive_job_id_(receive_job_id),
    from_physical_data_id_(from_physical_data_id),
    to_worker_id_(to_worker_id),
    to_ip_(to_ip), to_port_(to_port),
    before_set_(before) {
  name_ = REMOTE_COPY_SEND_NAME;
  type_ = REMOTE_COPY_SEND;
}

RemoteCopySendCommand::~RemoteCopySendCommand() {
}

SchedulerCommand* RemoteCopySendCommand::Clone() {
  return new RemoteCopySendCommand();
}

bool RemoteCopySendCommand::Parse(const std::string& data) {
  RemoteCopySendPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RemoteCopySendCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool RemoteCopySendCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_remote_send()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RemoteCopySendCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.remote_send());
  }
}

std::string RemoteCopySendCommand::toString() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the remote copy send field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_REMOTE_SEND);
  RemoteCopySendPBuf* cmd = buf.mutable_remote_send();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string RemoteCopySendCommand::toStringWTags() {
  std::string str;
  str += (name_ + ",");
  str += ("job-id:" + job_id_.toString() + ",");
  str += ("receive-job-id:" + receive_job_id_.toString() + ",");
  str += ("from-physical-data-id:" + from_physical_data_id_.toString() + ",");
  str += ("to-worker-id:" + to_worker_id_.toString() + ",");
  str += ("to-ip:" + to_ip_ + ",");
  str += ("to-port:" + to_port_.toString() + ",");
  str += ("before:" + before_set_.toString());
  return str;
}

ID<job_id_t> RemoteCopySendCommand::job_id() {
  return job_id_;
}

ID<job_id_t> RemoteCopySendCommand::receive_job_id() {
  return receive_job_id_;
}

ID<physical_data_id_t> RemoteCopySendCommand::from_physical_data_id() {
  return from_physical_data_id_;
}

ID<worker_id_t> RemoteCopySendCommand::to_worker_id() {
  return to_worker_id_;
}

std::string RemoteCopySendCommand::to_ip() {
  return to_ip_;
}

ID<port_t> RemoteCopySendCommand::to_port() {
  return to_port_;
}

IDSet<job_id_t> RemoteCopySendCommand::before_set() {
  return before_set_;
}

bool RemoteCopySendCommand::ReadFromProtobuf(const RemoteCopySendPBuf& buf) {
  job_id_.set_elem(buf.job_id());
  receive_job_id_.set_elem(buf.receive_job_id());
  from_physical_data_id_.set_elem(buf.from_physical_id());
  to_worker_id_.set_elem(buf.to_worker_id());
  to_ip_ = buf.to_ip();
  to_port_.set_elem(buf.to_port());
  before_set_.ConvertFromRepeatedField(buf.before_set().ids());
  return true;
}

bool RemoteCopySendCommand::WriteToProtobuf(RemoteCopySendPBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_receive_job_id(receive_job_id().elem());
  buf->set_from_physical_id(from_physical_data_id().elem());
  buf->set_to_worker_id(to_worker_id().elem());
  buf->set_to_ip(to_ip());
  buf->set_to_port(to_port().elem());
  before_set_.ConvertToRepeatedField(buf->mutable_before_set()->mutable_ids());
  return true;
}
