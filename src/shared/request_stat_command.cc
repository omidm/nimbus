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
  * Controller requests the worker stats including blocked, running, and idle
  * time with a RequestStatCommand.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/request_stat_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

RequestStatCommand::RequestStatCommand() {
  name_ = REQUEST_STAT_NAME;
  type_ = REQUEST_STAT;
}

RequestStatCommand::RequestStatCommand(const counter_t& query_id)
: query_id_(query_id) {
  name_ = REQUEST_STAT_NAME;
  type_ = REQUEST_STAT;
}

RequestStatCommand::~RequestStatCommand() {
}

SchedulerCommand* RequestStatCommand::Clone() {
  return new RequestStatCommand();
}

bool RequestStatCommand::Parse(const std::string& data) {
  RequestStatPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RequestStatCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool RequestStatCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_request_stat()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RequestStatCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.request_stat());
  }
}

std::string RequestStatCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_REQUEST_STAT);
  RequestStatPBuf* cmd = buf.mutable_request_stat();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string RequestStatCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("query_id:" + boost::lexical_cast<std::string>(query_id_));
  return str;
}


counter_t RequestStatCommand::query_id() {
  return query_id_;
}

bool RequestStatCommand::ReadFromProtobuf(const RequestStatPBuf& buf) {
  query_id_ = buf.query_id();
  return true;
}

bool RequestStatCommand::WriteToProtobuf(RequestStatPBuf* buf) {
  buf->set_query_id(query_id_);
  return true;
}

