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
  * A EndCommandTemplateCommand is a message sent from a controller to the
  * worker to mark the end of command template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/end_command_template_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

EndCommandTemplateCommand::EndCommandTemplateCommand() {
  name_ = END_COMMAND_TEMPLATE_NAME;
  type_ = END_COMMAND_TEMPLATE;
}

EndCommandTemplateCommand::EndCommandTemplateCommand(const std::string& command_template_name)
  : command_template_name_(command_template_name) {
  name_ = END_COMMAND_TEMPLATE_NAME;
  type_ = END_COMMAND_TEMPLATE;
}

EndCommandTemplateCommand::~EndCommandTemplateCommand() {
}

SchedulerCommand* EndCommandTemplateCommand::Clone() {
  return new EndCommandTemplateCommand();
}


bool EndCommandTemplateCommand::Parse(const std::string& data) {
  EndCommandTemplatePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse EndCommandTemplateCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool EndCommandTemplateCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_end_command_template()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse EndCommandTemplateCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.end_command_template());
  }
}

std::string EndCommandTemplateCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn job_graph field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_END_COMMAND_TEMPLATE);
  EndCommandTemplatePBuf* cmd = buf.mutable_end_command_template();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string EndCommandTemplateCommand::ToString() {
  // TODO(omidm) complete the implementation!
  std::string str;
  str += (name_ + " ");
  str += ("name:" + command_template_name_ + ", ...");
  return str;
}

std::string EndCommandTemplateCommand::command_template_name() {
  return command_template_name_;
}

bool EndCommandTemplateCommand::ReadFromProtobuf(const EndCommandTemplatePBuf& buf) {
  command_template_name_ = buf.command_template_name();

  return true;
}

bool EndCommandTemplateCommand::WriteToProtobuf(EndCommandTemplatePBuf* buf) {
  buf->set_command_template_name(command_template_name());

  return true;
}
