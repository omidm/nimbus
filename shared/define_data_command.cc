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
  * This command is sent from an application to the scheduler. It defines
  * a new logical data object (logical data ID). The application provides
  * the variable name (type), logical ID, its geometric region (defined
  * as an abstract partition ID) and what partitions neighbor it.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/define_data_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

DefineDataCommand::DefineDataCommand() {
  name_ = DEFINE_DATA_NAME;
  type_ = DEFINE_DATA;
}

DefineDataCommand::DefineDataCommand(const std::string& data_name,
                                     const ID<logical_data_id_t>& logical_data_id,
                                     const ID<partition_id_t>& partition_id,
                                     const IDSet<partition_id_t>& neighbor_partitions,
                                     const ID<job_id_t>& parent_job_id)
  : data_name_(data_name),
    logical_data_id_(logical_data_id),
    partition_id_(partition_id),
    neighbor_partitions_(neighbor_partitions),
    parent_job_id_(parent_job_id) {
  name_ = DEFINE_DATA_NAME;
  type_ = DEFINE_DATA;
}

DefineDataCommand::~DefineDataCommand() {
}

SchedulerCommand* DefineDataCommand::Clone() {
  return new DefineDataCommand();
}

bool DefineDataCommand::Parse(const std::string& params) {
  DefineDataPBuf buf;
  bool result = buf.ParseFromString(params);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse DefineDataCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool DefineDataCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_define_data()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse DefineDataCommand from DefineDataPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.define_data());
  }
}

std::string DefineDataCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_DEFINE_DATA);
  DefineDataPBuf* cmd = buf.mutable_define_data();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string DefineDataCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + data_name_ + " ");
  str += ("logical-id:" + logical_data_id_.ToNetworkData() + " ");
  str += ("partition-id:" + partition_id_.ToNetworkData() + " ");
  str += ("neighbor-partitions:" + neighbor_partitions_.ToNetworkData() + " ");
  str += ("parent-id:" + parent_job_id_.ToNetworkData() + " ");
  return str;
}

std::string DefineDataCommand::data_name() {
  return data_name_;
}

ID<job_id_t> DefineDataCommand::parent_job_id() {
  return parent_job_id_;
}

ID<logical_data_id_t> DefineDataCommand::logical_data_id() {
  return logical_data_id_;
}

ID<partition_id_t> DefineDataCommand::partition_id() {
  return partition_id_;
}

IDSet<partition_id_t> DefineDataCommand::neighbor_partitions() {
  return neighbor_partitions_;
}


bool DefineDataCommand::ReadFromProtobuf(const DefineDataPBuf& buf) {
  data_name_ = buf.name();
  logical_data_id_.set_elem(buf.logical_data_id());
  partition_id_.set_elem(buf.partition_id());
  neighbor_partitions_.ConvertFromRepeatedField(buf.neighbor_partitions().ids());
  parent_job_id_.set_elem(buf.parent_id());
  return true;
}

bool DefineDataCommand::WriteToProtobuf(DefineDataPBuf* buf) {
  buf->set_name(data_name());
  buf->set_logical_data_id(logical_data_id().elem());
  buf->set_partition_id(partition_id().elem());
  neighbor_partitions().ConvertToRepeatedField(buf->mutable_neighbor_partitions()->mutable_ids());
  buf->set_parent_id(parent_job_id().elem());
  return true;
}
