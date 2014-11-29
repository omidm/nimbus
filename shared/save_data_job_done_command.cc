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
  * SaveDataJobDone command to signal completion of a save data job.
  * In addition to normal job done it also sends the controller a handle that
  * can be used to load the data later.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/save_data_job_done_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

SaveDataJobDoneCommand::SaveDataJobDoneCommand() {
  name_ = SAVE_DATA_JOB_DONE_NAME;
  type_ = SAVE_DATA_JOB_DONE;
  run_time_ = 0;
  wait_time_ = 0;
  max_alloc_ = 0;
}

SaveDataJobDoneCommand::SaveDataJobDoneCommand(const ID<job_id_t>& job_id,
                                               const double run_time,
                                               const double wait_time,
                                               const size_t max_alloc,
                                               const std::string& handle)
  : job_id_(job_id),
    run_time_(run_time),
    wait_time_(wait_time),
    max_alloc_(max_alloc),
    handle_(handle) {
  name_ = SAVE_DATA_JOB_DONE_NAME;
  type_ = SAVE_DATA_JOB_DONE;
}

SaveDataJobDoneCommand::~SaveDataJobDoneCommand() {
}

SchedulerCommand* SaveDataJobDoneCommand::Clone() {
  return new SaveDataJobDoneCommand();
}

bool SaveDataJobDoneCommand::Parse(const std::string& data) {
  SaveDataJobDonePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SaveDataJobDoneCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}

bool SaveDataJobDoneCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_save_data_job_done()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse SaveDataJobDoneCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.save_data_job_done());
  }
}

std::string SaveDataJobDoneCommand::ToNetworkData() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the job done field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_SAVE_DATA_JOB_DONE);
  SaveDataJobDonePBuf* jdbuf = buf.mutable_save_data_job_done();
  WriteToProtobuf(jdbuf);

  buf.SerializeToString(&result);

  return result;
}

std::string SaveDataJobDoneCommand::ToString() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + job_id_.ToNetworkData() + " ");
  str += ("run_time: " + boost::lexical_cast<std::string>(run_time_) + " ");
  str += ("wait_time: " + boost::lexical_cast<std::string>(wait_time_) + " ");
  str += ("max_alloc: " + boost::lexical_cast<std::string>(max_alloc_) + " ");
  str += ("handle: " + handle_);
  return str;
}

ID<job_id_t> SaveDataJobDoneCommand::job_id() {
  return job_id_;
}

double SaveDataJobDoneCommand::run_time() {
  return run_time_;
}

double SaveDataJobDoneCommand::wait_time() {
  return wait_time_;
}

size_t SaveDataJobDoneCommand::max_alloc() {
  return max_alloc_;
}

std::string SaveDataJobDoneCommand::handle() {
  return handle_;
}

bool SaveDataJobDoneCommand::ReadFromProtobuf(const SaveDataJobDonePBuf& buf) {
  job_id_.set_elem(buf.job_id());
  run_time_ = buf.run_time();
  wait_time_ = buf.wait_time();
  max_alloc_ = buf.max_alloc();
  handle_ = buf.handle();
  return true;
}

bool SaveDataJobDoneCommand::WriteToProtobuf(SaveDataJobDonePBuf* buf) {
  buf->set_job_id(job_id().elem());
  buf->set_run_time(run_time());
  buf->set_wait_time(wait_time());
  buf->set_max_alloc(max_alloc());
  buf->set_handle(handle());
  return true;
}

