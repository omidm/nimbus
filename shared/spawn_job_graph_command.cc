/*
 * Copyright 2014 Stanford University.
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
  * A SpawnJobGraphCommand is a message sent from a worker to the
  * controller to spawn an instance of the job graph template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/spawn_job_graph_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SpawnJobGraphCommand::SpawnJobGraphCommand() {
  name_ = SPAWN_JOB_GRAPH_NAME;
  type_ = SPAWN_JOB_GRAPH;
}

SpawnJobGraphCommand::SpawnJobGraphCommand(const std::string& job_graph_name,
                                           const std::vector<job_id_t>& inner_job_ids,
                                           const std::vector<job_id_t>& outer_job_ids,
                                           const std::vector<Parameter>& parameters,
                                           const ID<job_id_t>& parent_job_id)
  : job_graph_name_(job_graph_name),
    inner_job_ids_(inner_job_ids),
    outer_job_ids_(inner_job_ids),
    parameters_(parameters),
    parent_job_id_(parent_job_id) {
  name_ = SPAWN_JOB_GRAPH_NAME;
  type_ = SPAWN_JOB_GRAPH;
}

SpawnJobGraphCommand::~SpawnJobGraphCommand() {
}

SchedulerCommand* SpawnJobGraphCommand::Clone() {
  return new SpawnJobGraphCommand();
}


bool SpawnJobGraphCommand::Parse(const std::string& data) {
  SubmitJobGraphPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SpawnJobGraphCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool SpawnJobGraphCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_submit_job_graph()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SpawnJobGraphCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.submit_job_graph());
  }
}

std::string SpawnJobGraphCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn job_graph field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_SPAWN_JOB_GRAPH);
  SubmitJobGraphPBuf* cmd = buf.mutable_submit_job_graph();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string SpawnJobGraphCommand::ToString() {
  // TODO(omidm) complete the implementation!
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_graph_name_ + ",");
  str += ("parent-id:" + parent_job_id_.ToNetworkData() + ",");
  return str;
}

std::string SpawnJobGraphCommand::job_graph_name() {
  return job_graph_name_;
}

std::vector<job_id_t> SpawnJobGraphCommand::inner_job_ids() {
  return inner_job_ids_;
}

std::vector<job_id_t> SpawnJobGraphCommand::outer_job_ids() {
  return outer_job_ids_;
}

std::vector<Parameter> SpawnJobGraphCommand::parameters() {
  return parameters_;
}

ID<job_id_t> SpawnJobGraphCommand::parent_job_id() {
  return parent_job_id_;
}


bool SpawnJobGraphCommand::ReadFromProtobuf(const SubmitJobGraphPBuf& buf) {
  job_graph_name_ = buf.job_graph_name();
  {
    inner_job_ids_.clear();
    typename google::protobuf::RepeatedField<job_id_t>::const_iterator it =
      buf.inner_job_ids().ids().begin();
    for (; it != buf.inner_job_ids().ids().end(); ++it) {
      inner_job_ids_.push_back(*it);
    }
  }
  {
    outer_job_ids_.clear();
    typename google::protobuf::RepeatedField<job_id_t>::const_iterator it =
      buf.outer_job_ids().ids().begin();
    for (; it != buf.outer_job_ids().ids().end(); ++it) {
      outer_job_ids_.push_back(*it);
    }
  }
  {
    parameters_.clear();
    typename google::protobuf::RepeatedPtrField<std::string >::const_iterator it =
      buf.parameters().params().begin();
    for (; it != buf.parameters().params().end(); ++it) {
      SerializedData d(*it);
      Parameter p;
      p.set_ser_data(d);
      parameters_.push_back(p);
    }
  }

  parent_job_id_.set_elem(buf.parent_job_id());

  return true;
}

bool SpawnJobGraphCommand::WriteToProtobuf(SubmitJobGraphPBuf* buf) {
  // TODO(omidm) complete the implementation!
  buf->set_job_graph_name(job_graph_name());
  buf->set_parent_job_id(parent_job_id().elem());

  return true;
}
