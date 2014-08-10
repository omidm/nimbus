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
  * A define partition command tells the scheduler about a particular
  * geometric region. Later data objects can then refer to this geometric
  * region. Using partition identifiers allows the scheduler to very
  * easily determine if two defined datas cover the same region and more
  * precisely defines the geometric boundaries the application is using.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "shared/define_partition_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

DefinePartitionCommand::DefinePartitionCommand() {
  name_ = DEFINE_PARTITION_NAME;
  type_ = DEFINE_PARTITION;
}

DefinePartitionCommand::DefinePartitionCommand(const ID<partition_id_t>& part,
                                               const GeometricRegion& r)
  :id_(part), region_(r) {
  name_ = DEFINE_PARTITION_NAME;
  type_ = DEFINE_PARTITION;
}

DefinePartitionCommand::~DefinePartitionCommand() {}

SchedulerCommand* DefinePartitionCommand::Clone() {
  return new DefinePartitionCommand();
}

bool DefinePartitionCommand::Parse(const std::string& data) {
  DefinePartitionPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse DefinePartitionCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool DefinePartitionCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_define_partition()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse DefinePartitionCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.define_partition());
  }
}

std::string DefinePartitionCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_DEFINE_PARTITION);
  DefinePartitionPBuf* cmd = buf.mutable_define_partition();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string DefinePartitionCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + id_.ToNetworkData() + " ");
  str += ("region:" + region_.ToNetworkData() + " ");
  return str;
}

ID<partition_id_t> DefinePartitionCommand::partition_id() {
  return id_;
}

const GeometricRegion* DefinePartitionCommand::region() {
  return &region_;
}

bool DefinePartitionCommand::ReadFromProtobuf(const DefinePartitionPBuf& buf) {
  id_.set_elem(buf.partition_id());
  region_.FillInValues(&buf.region());
  return true;
}

bool DefinePartitionCommand::WriteToProtobuf(DefinePartitionPBuf* buf) {
  buf->set_partition_id(id_.elem());
  region_.FillInMessage(buf->mutable_region());
  return true;
}
