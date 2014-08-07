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
  * Handshake command used to stablish connection between scheduler and worker
  * when worker initially joins the network. Scheduler registers the worker and
  * assigns unique worker id to the worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/handshake_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

HandshakeCommand::HandshakeCommand() {
  name_ = HANDSHAKE_NAME;
  type_ = HANDSHAKE;
}

HandshakeCommand::HandshakeCommand(const ID<worker_id_t>& worker_id,
                                   const std::string& ip,
                                   const ID<port_t>& port)
: worker_id_(worker_id),
  ip_(ip),
  port_(port) {
  name_ = HANDSHAKE_NAME;
  type_ = HANDSHAKE;
}

HandshakeCommand::~HandshakeCommand() {
}

SchedulerCommand* HandshakeCommand::Clone() {
  return new HandshakeCommand();
}

bool HandshakeCommand::Parse(const std::string& data) {
  HandshakePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SpawnComputeJobCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool HandshakeCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_handshake()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse HandshakeCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.handshake());
  }
}

std::string HandshakeCommand::toString() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the handshake field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_HANDSHAKE);
  HandshakePBuf* cmd = buf.mutable_handshake();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string HandshakeCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("worker-id:" + worker_id_.toString() + " ");
  str += ("ip:" + ip_ + " ");
  str += ("port:" + port_.toString());
  return str;
}


ID<worker_id_t> HandshakeCommand::worker_id() {
  return worker_id_;
}

std::string HandshakeCommand::ip() {
  return ip_;
}

ID<port_t> HandshakeCommand::port() {
  return port_;
}

bool HandshakeCommand::ReadFromProtobuf(const HandshakePBuf& buf) {
  worker_id_.set_elem(buf.worker_id());
  ip_ = buf.ip();
  port_.set_elem(buf.port());
  return true;
}

bool HandshakeCommand::WriteToProtobuf(HandshakePBuf* buf) {
  buf->set_worker_id(worker_id().elem());
  buf->set_ip(ip());
  buf->set_port(port().elem());
  return true;
}
