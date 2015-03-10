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
  * worker responds to controller query for its stats with a
  * RespondStatCommand, reporting idle, block and run time.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/respond_stat_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

RespondStatCommand::RespondStatCommand() {
  name_ = RESPOND_STAT_NAME;
  type_ = RESPOND_STAT;
}

RespondStatCommand::RespondStatCommand(const counter_t& query_id,
                                       const int64_t& run_time,
                                       const int64_t& block_time,
                                       const int64_t& idle_time)
  : query_id_(query_id),
    run_time_(run_time),
    block_time_(block_time),
    idle_time_(idle_time) {
  name_ = RESPOND_STAT_NAME;
  type_ = RESPOND_STAT;
}

RespondStatCommand::~RespondStatCommand() {
}

SchedulerCommand* RespondStatCommand::Clone() {
  return new RespondStatCommand();
}

bool RespondStatCommand::Parse(const std::string& data) {
  RespondStatPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RespondStatCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool RespondStatCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_respond_stat()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse RespondStatCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.respond_stat());
  }
}

std::string RespondStatCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the job done field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_RESPOND_STAT);
  RespondStatPBuf* jdbuf = buf.mutable_respond_stat();
  WriteToProtobuf(jdbuf);

  buf.SerializeToString(&result);

  return result;
}

std::string RespondStatCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("query_id:" + boost::lexical_cast<std::string>(query_id_) + " ");
  str += ("run_time: " + boost::lexical_cast<std::string>(run_time_) + " ");
  str += ("block_time: " + boost::lexical_cast<std::string>(block_time_) + " ");
  str += ("idle_time: " + boost::lexical_cast<std::string>(idle_time_));
  return str;
}

counter_t RespondStatCommand::query_id() {
  return query_id_;
}

int64_t RespondStatCommand::run_time() {
  return run_time_;
}

int64_t RespondStatCommand::block_time() {
  return block_time_;
}

int64_t RespondStatCommand::idle_time() {
  return idle_time_;
}

bool RespondStatCommand::ReadFromProtobuf(const RespondStatPBuf& buf) {
  query_id_ = buf.query_id();
  run_time_ = buf.run_time();
  block_time_ = buf.block_time();
  idle_time_ = buf.idle_time();
  return true;
}

bool RespondStatCommand::WriteToProtobuf(RespondStatPBuf* buf) {
  buf->set_query_id(query_id());
  buf->set_run_time(run_time());
  buf->set_block_time(block_time());
  buf->set_idle_time(idle_time());
  return true;
}

