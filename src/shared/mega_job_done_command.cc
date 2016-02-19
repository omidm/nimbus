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
  * A remote copy operation between two workers has two jobs: the
  * send and receive. This is the receive half.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "src/shared/mega_job_done_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

MegaJobDoneCommand::MegaJobDoneCommand() {
  name_ = MEGA_JOB_DONE_NAME;
  type_ = MEGA_JOB_DONE;
}

MegaJobDoneCommand::MegaJobDoneCommand(const std::vector<job_id_t>& job_ids)
: job_ids_(job_ids) {
  name_ = MEGA_JOB_DONE_NAME;
  type_ = MEGA_JOB_DONE;
}

MegaJobDoneCommand::~MegaJobDoneCommand() {
}

SchedulerCommand* MegaJobDoneCommand::Clone() {
  return new MegaJobDoneCommand();
}

bool MegaJobDoneCommand::Parse(const std::string& data) {
  MegaJobDonePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse MegaJobDoneCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool MegaJobDoneCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_mega_job_done()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse MegaJobDoneCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.mega_job_done());
  }
}

std::string MegaJobDoneCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_MEGA_JOB_DONE);
  MegaJobDonePBuf* cmd = buf.mutable_mega_job_done();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string MegaJobDoneCommand::ToString() {
  // TODO(omidm) complete the implementation!
  std::string str;
  str += (name_ + ", ...");
  return str;
}

std::vector<job_id_t> MegaJobDoneCommand::job_ids() {
  return job_ids_;
}

const std::vector<job_id_t>* MegaJobDoneCommand::job_ids_p() {
  return &job_ids_;
}

bool MegaJobDoneCommand::ReadFromProtobuf(const MegaJobDonePBuf& buf) {
  {
    job_ids_.clear();
    typename google::protobuf::RepeatedField<job_id_t>::const_iterator it =
      buf.job_ids().ids().begin();
    for (; it != buf.job_ids().ids().end(); ++it) {
      job_ids_.push_back(*it);
    }
  }

  return true;
}

bool MegaJobDoneCommand::WriteToProtobuf(MegaJobDonePBuf* buf) {
  {
    typename google::protobuf::RepeatedField<job_id_t> *b =
      buf->mutable_job_ids()->mutable_ids();
    std::vector<job_id_t>::const_iterator it = job_ids_.begin();
    for (; it != job_ids_.end(); ++it) {
      b->Add(*it);
    }
  }

  return true;
}

