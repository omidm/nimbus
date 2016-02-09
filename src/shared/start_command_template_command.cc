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
  * A StartCommandTemplateCommand is a message sent from a controller to the
  * worker to mark the definition of command template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/start_command_template_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

StartCommandTemplateCommand::StartCommandTemplateCommand() {
  name_ = START_COMMAND_TEMPLATE_NAME;
  type_ = START_COMMAND_TEMPLATE;
}

StartCommandTemplateCommand::StartCommandTemplateCommand(const std::string& command_template_name,
                                                   const std::vector<job_id_t>& inner_job_ids,
                                                   const std::vector<job_id_t>& outer_job_ids,
                                                   const std::vector<physical_data_id_t>& phy_ids)
  : command_template_name_(command_template_name),
    inner_job_ids_(inner_job_ids),
    outer_job_ids_(outer_job_ids),
    phy_ids_(phy_ids) {
  name_ = START_COMMAND_TEMPLATE_NAME;
  type_ = START_COMMAND_TEMPLATE;
}

StartCommandTemplateCommand::~StartCommandTemplateCommand() {
}

SchedulerCommand* StartCommandTemplateCommand::Clone() {
  return new StartCommandTemplateCommand();
}


bool StartCommandTemplateCommand::Parse(const std::string& data) {
  StartCommandTemplatePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse StartCommandTemplateCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool StartCommandTemplateCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_start_command_template()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse StartCommandTemplateCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.start_command_template());
  }
}

std::string StartCommandTemplateCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn job_graph field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_START_COMMAND_TEMPLATE);
  StartCommandTemplatePBuf* cmd = buf.mutable_start_command_template();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string StartCommandTemplateCommand::ToString() {
  // TODO(omidm) complete the implementation!
  std::string str;
  str += (name_ + " ");
  str += ("name:" + command_template_name_ + ", ...");
  return str;
}

std::string StartCommandTemplateCommand::command_template_name() {
  return command_template_name_;
}

std::vector<job_id_t> StartCommandTemplateCommand::inner_job_ids() {
  return inner_job_ids_;
}

std::vector<job_id_t> StartCommandTemplateCommand::outer_job_ids() {
  return outer_job_ids_;
}

std::vector<job_id_t> StartCommandTemplateCommand::phy_ids() {
  return phy_ids_;
}

bool StartCommandTemplateCommand::ReadFromProtobuf(const StartCommandTemplatePBuf& buf) {
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
    phy_ids_.clear();
    typename google::protobuf::RepeatedField<physical_data_id_t>::const_iterator it =
      buf.phy_ids().ids().begin();
    for (; it != buf.phy_ids().ids().end(); ++it) {
      phy_ids_.push_back(*it);
    }
  }

  return true;
}

bool StartCommandTemplateCommand::WriteToProtobuf(StartCommandTemplatePBuf* buf) {
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
    typename google::protobuf::RepeatedField<physical_data_id_t> *b =
      buf->mutable_phy_ids()->mutable_ids();
    std::vector<physical_data_id_t>::const_iterator it = phy_ids_.begin();
    for (; it != phy_ids_.end(); ++it) {
      b->Add(*it);
    }
  }

  return true;
}
