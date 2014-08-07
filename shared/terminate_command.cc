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
  * Terminate command used to signal the scheduler that the application is
  * complete and there  would be no more spawned jobs.
  * Also used by the scheduler to terminate the workers.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/terminate_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

TerminateCommand::TerminateCommand() {
  name_ = TERMINATE_NAME;
  type_ = TERMINATE;
}

TerminateCommand::TerminateCommand(const ID<exit_status_t>& exit_status)
: exit_status_(exit_status) {
  name_ = TERMINATE_NAME;
  type_ = TERMINATE;
}

TerminateCommand::~TerminateCommand() {
}

SchedulerCommand* TerminateCommand::Clone() {
  return new TerminateCommand();
}

bool TerminateCommand::Parse(const std::string& data) {
  TerminatePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse TerminateCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool TerminateCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_terminate()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse TerminateCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.terminate());
  }
}

std::string TerminateCommand::toString() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_TERMINATE);
  TerminatePBuf* cmd = buf.mutable_terminate();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string TerminateCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("exit-status:" + exit_status_.toString());
  return str;
}


ID<exit_status_t> TerminateCommand::exit_status() {
  return exit_status_;
}

bool TerminateCommand::ReadFromProtobuf(const TerminatePBuf& buf) {
  exit_status_.set_elem(buf.exit_status());
  return true;
}

bool TerminateCommand::WriteToProtobuf(TerminatePBuf* buf) {
  buf->set_exit_status(exit_status().elem());
  return true;
}
