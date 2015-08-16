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
  * A SpawnCommandTemplateCommand is a message sent from a controller to the
  * worker to instantiate a command template at worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/spawn_command_template_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SpawnCommandTemplateCommand::SpawnCommandTemplateCommand() {
  name_ = SPAWN_COMMAND_TEMPLATE_NAME;
  type_ = SPAWN_COMMAND_TEMPLATE;
}

SpawnCommandTemplateCommand::SpawnCommandTemplateCommand(const std::string& command_template_name,
                                                   const std::vector<job_id_t>& inner_job_ids,
                                                   const std::vector<job_id_t>& outer_job_ids,
                                                   const std::vector<Parameter>& parameters,
                                                   const std::vector<physical_data_id_t>& phy_ids,
                                                   const template_id_t& template_generation_id)
  : command_template_name_(command_template_name),
    inner_job_ids_(inner_job_ids),
    outer_job_ids_(outer_job_ids),
    parameters_(parameters),
    phy_ids_(phy_ids),
    template_generation_id_(template_generation_id) {
  name_ = SPAWN_COMMAND_TEMPLATE_NAME;
  type_ = SPAWN_COMMAND_TEMPLATE;
}

SpawnCommandTemplateCommand::~SpawnCommandTemplateCommand() {
}

SchedulerCommand* SpawnCommandTemplateCommand::Clone() {
  return new SpawnCommandTemplateCommand();
}


bool SpawnCommandTemplateCommand::Parse(const std::string& data) {
  SpawnCommandTemplatePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SpawnCommandTemplateCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool SpawnCommandTemplateCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_spawn_command_template()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SpawnCommandTemplateCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.spawn_command_template());
  }
}

std::string SpawnCommandTemplateCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn job_graph field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_SPAWN_COMMAND_TEMPLATE);
  SpawnCommandTemplatePBuf* cmd = buf.mutable_spawn_command_template();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string SpawnCommandTemplateCommand::ToString() {
  // TODO(omidm) complete the implementation!
  std::string str;
  str += (name_ + " ");
  str += ("name:" + command_template_name_ + ", ...");
  return str;
}

std::string SpawnCommandTemplateCommand::command_template_name() {
  return command_template_name_;
}

std::vector<job_id_t> SpawnCommandTemplateCommand::inner_job_ids() {
  return inner_job_ids_;
}

std::vector<job_id_t> SpawnCommandTemplateCommand::outer_job_ids() {
  return outer_job_ids_;
}

std::vector<Parameter> SpawnCommandTemplateCommand::parameters() {
  return parameters_;
}

std::vector<job_id_t> SpawnCommandTemplateCommand::phy_ids() {
  return phy_ids_;
}

template_id_t SpawnCommandTemplateCommand::template_generation_id() {
  return template_generation_id_;
}

bool SpawnCommandTemplateCommand::ReadFromProtobuf(const SpawnCommandTemplatePBuf& buf) {
  command_template_name_ = buf.command_template_name();

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
  {
    phy_ids_.clear();
    typename google::protobuf::RepeatedField<physical_data_id_t>::const_iterator it =
      buf.phy_ids().ids().begin();
    for (; it != buf.phy_ids().ids().end(); ++it) {
      phy_ids_.push_back(*it);
    }
  }

  template_generation_id_ = buf.template_generation_id();

  return true;
}

bool SpawnCommandTemplateCommand::WriteToProtobuf(SpawnCommandTemplatePBuf* buf) {
  buf->set_command_template_name(command_template_name());

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
  {
    typename google::protobuf::RepeatedField<physical_data_id_t> *b =
      buf->mutable_phy_ids()->mutable_ids();
    std::vector<physical_data_id_t>::const_iterator it = phy_ids_.begin();
    for (; it != phy_ids_.end(); ++it) {
      b->Add(*it);
    }
  }

  buf->set_template_generation_id(template_generation_id_);

  return true;
}
