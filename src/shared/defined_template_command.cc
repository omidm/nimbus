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
  * A DefinedTemplateCommand is a message sent from a worker to the
  * controller to mark the definition of a job graph template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/defined_template_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

DefinedTemplateCommand::DefinedTemplateCommand() {
  name_ = DEFINED_TEMPLATE_NAME;
  type_ = DEFINED_TEMPLATE;
}

DefinedTemplateCommand::DefinedTemplateCommand(const std::string& job_graph_name)
  : job_graph_name_(job_graph_name) {
  name_ = DEFINED_TEMPLATE_NAME;
  type_ = DEFINED_TEMPLATE;
}

DefinedTemplateCommand::~DefinedTemplateCommand() {
}

SchedulerCommand* DefinedTemplateCommand::Clone() {
  return new DefinedTemplateCommand();
}


bool DefinedTemplateCommand::Parse(const std::string& data) {
  DefinedTemplatePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse DefinedTemplateCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool DefinedTemplateCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_defined_template()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse DefinedTemplateCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.defined_template());
  }
}

std::string DefinedTemplateCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn job_graph field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_DEFINED_TEMPLATE);
  DefinedTemplatePBuf* cmd = buf.mutable_defined_template();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string DefinedTemplateCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_graph_name_);
  return str;
}

std::string DefinedTemplateCommand::job_graph_name() {
  return job_graph_name_;
}

bool DefinedTemplateCommand::ReadFromProtobuf(const DefinedTemplatePBuf& buf) {
  job_graph_name_ = buf.job_graph_name();
  return true;
}

bool DefinedTemplateCommand::WriteToProtobuf(DefinedTemplatePBuf* buf) {
  buf->set_job_graph_name(job_graph_name());
  return true;
}

