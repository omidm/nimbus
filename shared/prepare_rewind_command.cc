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
  * Job done command to signal completion of a job.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/prepare_rewind_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

PrepareRewindCommand::PrepareRewindCommand() {
  name_ = PREPARE_REWIND_NAME;
  type_ = PREPARE_REWIND;
  worker_id_ = ID<worker_id_t>(NIMBUS_SCHEDULER_ID);
  checkpoint_id_ = ID<checkpoint_id_t>(NIMBUS_INIT_CHECKPOINT_ID);
}

PrepareRewindCommand::PrepareRewindCommand(const ID<worker_id_t>& worker_id,
                                           const ID<checkpoint_id_t>& checkpoint_id)
  : worker_id_(worker_id),
    checkpoint_id_(checkpoint_id) {
  name_ = PREPARE_REWIND_NAME;
  type_ = PREPARE_REWIND;
}

PrepareRewindCommand::~PrepareRewindCommand() {
}

SchedulerCommand* PrepareRewindCommand::Clone() {
  return new PrepareRewindCommand();
}

bool PrepareRewindCommand::Parse(const std::string& data) {
  PrepareRewindPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse PrepareRewindCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool PrepareRewindCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_prepare_rewind()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse PrepareRewindCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.prepare_rewind());
  }
}

std::string PrepareRewindCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the job done field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_PREPARE_REWIND);
  PrepareRewindPBuf* jdbuf = buf.mutable_prepare_rewind();
  WriteToProtobuf(jdbuf);

  buf.SerializeToString(&result);

  return result;
}

std::string PrepareRewindCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("worker_id:" + worker_id_.ToNetworkData() + " ");
  str += ("checkpoint_id:" + checkpoint_id_.ToNetworkData());
  return str;
}

ID<worker_id_t> PrepareRewindCommand::worker_id() {
  return worker_id_;
}

ID<checkpoint_id_t> PrepareRewindCommand::checkpoint_id() {
  return checkpoint_id_;
}

bool PrepareRewindCommand::ReadFromProtobuf(const PrepareRewindPBuf& buf) {
  worker_id_.set_elem(buf.worker_id());
  checkpoint_id_.set_elem(buf.checkpoint_id());
  return true;
}

bool PrepareRewindCommand::WriteToProtobuf(PrepareRewindPBuf* buf) {
  buf->set_worker_id(worker_id().elem());
  buf->set_checkpoint_id(checkpoint_id().elem());
  return true;
}

