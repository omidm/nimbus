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
  * A StartTemplateCommand is a message sent from a worker to the
  * controller to mark the definition of a job graph template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/start_template_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

StartTemplateCommand::StartTemplateCommand() {
  name_ = START_TEMPLATE_NAME;
  type_ = START_TEMPLATE;
}

StartTemplateCommand::StartTemplateCommand(const std::string& job_graph_name,
                                           const ID<job_id_t>& parent_job_id)
  : job_graph_name_(job_graph_name),
    parent_job_id_(parent_job_id) {
  name_ = START_TEMPLATE_NAME;
  type_ = START_TEMPLATE;
}

StartTemplateCommand::~StartTemplateCommand() {
}

SchedulerCommand* StartTemplateCommand::Clone() {
  return new StartTemplateCommand();
}


bool StartTemplateCommand::Parse(const std::string& data) {
  StartTemplatePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse StartTemplateCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool StartTemplateCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_start_template()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse StartTemplateCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.start_template());
  }
}

std::string StartTemplateCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn job_graph field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_START_TEMPLATE);
  StartTemplatePBuf* cmd = buf.mutable_start_template();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string StartTemplateCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_graph_name_ + ",");
  str += ("parent-id:" + parent_job_id_.ToNetworkData());
  return str;
}

std::string StartTemplateCommand::job_graph_name() {
  return job_graph_name_;
}

ID<job_id_t> StartTemplateCommand::parent_job_id() {
  return parent_job_id_;
}

bool StartTemplateCommand::ReadFromProtobuf(const StartTemplatePBuf& buf) {
  job_graph_name_ = buf.job_graph_name();
  parent_job_id_.set_elem(buf.parent_job_id());
  return true;
}

bool StartTemplateCommand::WriteToProtobuf(StartTemplatePBuf* buf) {
  buf->set_job_graph_name(job_graph_name());
  buf->set_parent_job_id(parent_job_id().elem());
  return true;
}
