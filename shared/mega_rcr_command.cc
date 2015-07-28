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

#include "shared/mega_rcr_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

MegaRCRCommand::MegaRCRCommand() {
  name_ = MEGA_RCR_NAME;
  type_ = MEGA_RCR;
}

MegaRCRCommand::MegaRCRCommand(const ID<job_id_t>& job_id,
                               const std::vector<job_id_t>& receive_job_ids,
                               const std::vector<physical_data_id_t>& to_physical_data_ids)
: job_id_(job_id),
  receive_job_ids_(receive_job_ids),
  to_physical_data_ids_(to_physical_data_ids) {
  name_ = MEGA_RCR_NAME;
  type_ = MEGA_RCR;
}

MegaRCRCommand::~MegaRCRCommand() {
}

SchedulerCommand* MegaRCRCommand::Clone() {
  return new MegaRCRCommand();
}

bool MegaRCRCommand::Parse(const std::string& data) {
  MegaRCRPBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse MegaRCRCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool MegaRCRCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_mega_rcr()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse MegaRCRCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.mega_rcr());
  }
}

std::string MegaRCRCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_MEGA_RCR);
  MegaRCRPBuf* cmd = buf.mutable_mega_rcr();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}

std::string MegaRCRCommand::ToString() {
  // TODO(omidm) complete the implementation!
  std::string str;
  str += (name_ + ",");
  str += ("job-id:" + job_id_.ToNetworkData() + ", ...");
  return str;
}

ID<job_id_t> MegaRCRCommand::job_id() {
  return job_id_;
}

std::vector<job_id_t> MegaRCRCommand::receive_job_ids() {
  return receive_job_ids_;
}

const std::vector<job_id_t>* MegaRCRCommand::receive_job_ids_p() {
  return &receive_job_ids_;
}

std::vector<physical_data_id_t> MegaRCRCommand::to_physical_data_ids() {
  return to_physical_data_ids_;
}

const std::vector<physical_data_id_t>* MegaRCRCommand::to_physical_data_ids_p() {
  return &to_physical_data_ids_;
}

bool MegaRCRCommand::ReadFromProtobuf(const MegaRCRPBuf& buf) {
  job_id_.set_elem(buf.job_id());

  {
    receive_job_ids_.clear();
    typename google::protobuf::RepeatedField<job_id_t>::const_iterator it =
      buf.receive_job_ids().ids().begin();
    for (; it != buf.receive_job_ids().ids().end(); ++it) {
      receive_job_ids_.push_back(*it);
    }
  }

  {
    to_physical_data_ids_.clear();
    typename google::protobuf::RepeatedField<physical_data_id_t>::const_iterator it =
      buf.to_phy_ids().ids().begin();
    for (; it != buf.to_phy_ids().ids().end(); ++it) {
      to_physical_data_ids_.push_back(*it);
    }
  }

  return true;
}

bool MegaRCRCommand::WriteToProtobuf(MegaRCRPBuf* buf) {
  buf->set_job_id(job_id().elem());

  {
    typename google::protobuf::RepeatedField<job_id_t> *b =
      buf->mutable_receive_job_ids()->mutable_ids();
    std::vector<job_id_t>::const_iterator it = receive_job_ids_.begin();
    for (; it != receive_job_ids_.end(); ++it) {
      b->Add(*it);
    }
  }

  {
    typename google::protobuf::RepeatedField<physical_data_id_t> *b =
      buf->mutable_to_phy_ids()->mutable_ids();
    std::vector<physical_data_id_t>::const_iterator it = to_physical_data_ids_.begin();
    for (; it != to_physical_data_ids_.end(); ++it) {
      b->Add(*it);
    }
  }

  return true;
}

