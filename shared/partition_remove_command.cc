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
  * Command sent from scheduler to worker, to remove a partition
  * definition. Currently not used by scheduler.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "shared/partition_remove_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

PartitionRemoveCommand::PartitionRemoveCommand() {
  name_ = PARTITION_REMOVE_NAME;
  type_ = PARTITION_REMOVE;
}

PartitionRemoveCommand::PartitionRemoveCommand(const ID<partition_id_t>& part)
  :id_(part) {
  name_ = PARTITION_REMOVE_NAME;
  type_ = PARTITION_REMOVE;
}

PartitionRemoveCommand::~PartitionRemoveCommand() {}

SchedulerCommand* PartitionRemoveCommand::Clone() {
  return new PartitionRemoveCommand();
}

bool PartitionRemoveCommand::Parse(const std::string& data) {
  PartitionRemovePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse PartitionRemovePBuf from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool PartitionRemoveCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_partition_remove()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse PartitionRemoveCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.partition_remove());
  }
}

std::string PartitionRemoveCommand::toString() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_PARTITION_REMOVE);
  PartitionRemovePBuf* cmd = buf.mutable_partition_remove();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string PartitionRemoveCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + id_.toString());
  return str;
}

ID<partition_id_t> PartitionRemoveCommand::id() {
  return id_;
}


bool PartitionRemoveCommand::ReadFromProtobuf(const PartitionRemovePBuf& buf) {
  id_.set_elem(buf.id());
  return true;
}

bool PartitionRemoveCommand::WriteToProtobuf(PartitionRemovePBuf* buf) {
  buf->set_id(id().elem());
  return true;
}
