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
  * A SpawnTemplateCommand is a message sent from a worker to the
  * controller to spawn an instance of the job graph template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/spawn_template_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SpawnTemplateCommand::SpawnTemplateCommand() {
  name_ = SPAWN_TEMPLATE_NAME;
  type_ = SPAWN_TEMPLATE;
}

SpawnTemplateCommand::SpawnTemplateCommand(const std::string& job_graph_name,
                                           const std::vector<job_id_t>& inner_job_ids,
                                           const std::vector<job_id_t>& outer_job_ids,
                                           const std::vector<Parameter>& parameters,
                                           const ID<job_id_t>& parent_job_id)
  : job_graph_name_(job_graph_name),
    inner_job_ids_(inner_job_ids),
    outer_job_ids_(outer_job_ids),
    parameters_(parameters),
    parent_job_id_(parent_job_id) {
  name_ = SPAWN_TEMPLATE_NAME;
  type_ = SPAWN_TEMPLATE;
}

SpawnTemplateCommand::~SpawnTemplateCommand() {
}

SchedulerCommand* SpawnTemplateCommand::Clone() {
  return new SpawnTemplateCommand();
}


bool SpawnTemplateCommand::Parse(const std::string& data) {
  SpawnTemplatePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SpawnTemplateCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool SpawnTemplateCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_spawn_template()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SpawnTemplateCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.spawn_template());
  }
}

std::string SpawnTemplateCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn template field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_SPAWN_TEMPLATE);
  SpawnTemplatePBuf* cmd = buf.mutable_spawn_template();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string SpawnTemplateCommand::ToString() {
  // TODO(omidm) complete the implementation!
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_graph_name_ + ",");
  str += ("parent-id:" + parent_job_id_.ToNetworkData() + ",");
  return str;
}

std::string SpawnTemplateCommand::job_graph_name() {
  return job_graph_name_;
}

std::vector<job_id_t> SpawnTemplateCommand::inner_job_ids() {
  return inner_job_ids_;
}

std::vector<job_id_t> SpawnTemplateCommand::outer_job_ids() {
  return outer_job_ids_;
}

std::vector<Parameter> SpawnTemplateCommand::parameters() {
  return parameters_;
}

ID<job_id_t> SpawnTemplateCommand::parent_job_id() {
  return parent_job_id_;
}


bool SpawnTemplateCommand::ReadFromProtobuf(const SpawnTemplatePBuf& buf) {
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

bool SpawnTemplateCommand::WriteToProtobuf(SpawnTemplatePBuf* buf) {
  buf->set_job_graph_name(job_graph_name());

  {
    typename google::protobuf::RepeatedField<job_id_t> *b =
      buf->mutable_inner_job_ids()->mutable_ids();
    std::vector<job_id_t>::const_iterator it = inner_job_ids_.begin();
    for (; it != inner_job_ids_.end(); ++it) {
      b->Add(*it);
    }
  }
  {
    typename google::protobuf::RepeatedField<job_id_t> *b =
      buf->mutable_outer_job_ids()->mutable_ids();
    std::vector<job_id_t>::const_iterator it = outer_job_ids_.begin();
    for (; it != outer_job_ids_.end(); ++it) {
      b->Add(*it);
    }
  }
  {
    ParameterVector *b = buf->mutable_parameters();
    std::vector<Parameter>::iterator it = parameters_.begin();
    for (; it != parameters_.end(); ++it) {
      b->add_params(it->ser_data().ToNetworkData());
    }
  }

  buf->set_parent_job_id(parent_job_id().elem());

  return true;
}
