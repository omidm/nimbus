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
  * Controller ask worker to prints the stats including blocked, running,
  * idle time, and parent execution time with a PrintStatCommand.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/print_stat_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

PrintStatCommand::PrintStatCommand() {
  name_ = PRINT_STAT_NAME;
  type_ = PRINT_STAT;
}

PrintStatCommand::PrintStatCommand(const counter_t& query_id)
: query_id_(query_id) {
  name_ = PRINT_STAT_NAME;
  type_ = PRINT_STAT;
}

PrintStatCommand::~PrintStatCommand() {
}

SchedulerCommand* PrintStatCommand::Clone() {
  return new PrintStatCommand();
}

bool PrintStatCommand::Parse(const std::string& data) {
  PrintStatPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse PrintStatCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool PrintStatCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_print_stat()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse PrintStatCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.print_stat());
  }
}

std::string PrintStatCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_PRINT_STAT);
  PrintStatPBuf* cmd = buf.mutable_print_stat();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string PrintStatCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("query_id:" + boost::lexical_cast<std::string>(query_id_));
  return str;
}


counter_t PrintStatCommand::query_id() {
  return query_id_;
}

bool PrintStatCommand::ReadFromProtobuf(const PrintStatPBuf& buf) {
  query_id_ = buf.query_id();
  return true;
}

bool PrintStatCommand::WriteToProtobuf(PrintStatPBuf* buf) {
  buf->set_query_id(query_id_);
  return true;
}

