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
  * Worker down command to signal the controller from server that the worker is down.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/worker_down_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

WorkerDownCommand::WorkerDownCommand() {
  name_ = WORKER_DOWN_NAME;
  type_ = WORKER_DOWN;
  worker_id_ = ID<worker_id_t>(NIMBUS_SCHEDULER_ID);
}

WorkerDownCommand::WorkerDownCommand(const ID<worker_id_t>& worker_id)
  : worker_id_(worker_id) {
  name_ = WORKER_DOWN_NAME;
  type_ = WORKER_DOWN;
}

WorkerDownCommand::~WorkerDownCommand() {
}

SchedulerCommand* WorkerDownCommand::Clone() {
  return new WorkerDownCommand();
}

bool WorkerDownCommand::Parse(const std::string& data) {
  WorkerDownPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse WorkerDownCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool WorkerDownCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_worker_down()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse WorkerDownCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.worker_down());
  }
}

std::string WorkerDownCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the job done field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_WORKER_DOWN);
  WorkerDownPBuf* jdbuf = buf.mutable_worker_down();
  WriteToProtobuf(jdbuf);

  buf.SerializeToString(&result);

  return result;
}

std::string WorkerDownCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("worker_id:" + worker_id_.ToNetworkData());
  return str;
}

ID<worker_id_t> WorkerDownCommand::worker_id() {
  return worker_id_;
}

bool WorkerDownCommand::ReadFromProtobuf(const WorkerDownPBuf& buf) {
  worker_id_.set_elem(buf.worker_id());
  return true;
}

bool WorkerDownCommand::WriteToProtobuf(WorkerDownPBuf* buf) {
  buf->set_worker_id(worker_id().elem());
  return true;
}

