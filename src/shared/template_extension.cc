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
  * This is TemplateExtension module to hold the metadata for migrated jobs
  * within the execution template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/shared/template_extension.h"
#include "src/shared/execution_template.h"
#include "src/worker/job.h"

using namespace nimbus; // NOLINT

TemplateExtension::TemplateExtension(
    const bool& migrate_out,
    boost::shared_ptr<ComputeJobCommand> compute_command,
    const std::vector<boost::shared_ptr<RemoteCopySendCommand> >& send_commands,
    const std::vector<boost::shared_ptr<RemoteCopyReceiveCommand> >& receive_commands)
  : migrate_out_(migrate_out),
    compute_command_(compute_command),
    send_commands_(send_commands),
    receive_commands_(receive_commands) {
      job_done_counter_ = 0;
}

TemplateExtension::~TemplateExtension() {
}

bool TemplateExtension::migrate_out() const {
  return migrate_out_;
}

boost::shared_ptr<ComputeJobCommand> TemplateExtension::compute_command() const {
  return compute_command_;
}

std::vector<boost::shared_ptr<RemoteCopySendCommand> >*
TemplateExtension::send_commands_p() {
  return &send_commands_;
}

std::vector<boost::shared_ptr<RemoteCopyReceiveCommand> >*
TemplateExtension::receive_commands_p() {
  return &receive_commands_;
}




size_t TemplateExtension::MarkJobDone(JobList *ready_jobs,
                                      ExecutionTemplate *et,
                                      State *state,
                                      bool append) {
  ++job_done_counter_;

  if (migrate_out_) {
    if (job_done_counter_ == (send_commands_.size() + receive_commands_.size())) {
      *state = RELEASE;
      return 0;
    } else {
      *state = BLOCK;
      return 0;
    }
  } else {
    *state = BLOCK;
    if (job_done_counter_ == receive_commands_.size()) {
      return LoadComputeJob(ready_jobs, et, append);
    } else if (job_done_counter_ == (receive_commands_.size() + 1)) {
      return LoadSendJobs(ready_jobs, et, append);
    } else {
      return 0;
    }
  }
}




size_t TemplateExtension::LoadComputeJob(JobList *ready_jobs,
                                         ExecutionTemplate *et,
                                         bool append) {
  job_id_t shadow_job_id = compute_command_->job_id().elem();
  if (!append) {
    ready_jobs->clear();
  }

  Job* j = et->application()->CloneJob(compute_command_->job_name());
  j->set_name("Compute:" + compute_command_->job_name());
  j->set_id(compute_command_->job_id());
  j->set_read_set(compute_command_->read_set());
  j->set_write_set(compute_command_->write_set());
  j->set_scratch_set(compute_command_->scratch_set());
  j->set_reduce_set(compute_command_->reduce_set());
  j->set_after_set(compute_command_->after_set());
  assert(compute_command_->before_set().size() == 0);
  j->set_sterile(compute_command_->sterile());
  j->set_region(compute_command_->region());
  j->set_future_job_id(compute_command_->future_job_id());

  j->set_parameters(compute_command_->params());
  j->set_shadow_job_id(shadow_job_id);
  j->set_execution_template(et);

  ready_jobs->push_back((j));

  return 1;
}


size_t TemplateExtension::LoadSendJobs(JobList *ready_jobs,
                                       ExecutionTemplate *et,
                                       bool append) {
  job_id_t shadow_job_id = compute_command_->job_id().elem();
  if (!append) {
    ready_jobs->clear();
  }

  size_t count = 0;
  std::vector<boost::shared_ptr<RemoteCopySendCommand> >::iterator iter =
    send_commands_.begin();
  for (; iter != send_commands_.end(); ++iter) {
    boost::shared_ptr<RemoteCopySendCommand> cm = *iter;
    et->data_exchanger()->AddContactInfo(cm->to_worker_id().elem(),
                                         cm->to_ip(),
                                         cm->to_port().elem());
    RemoteCopySendJob * j = new RemoteCopySendJob(et->data_exchanger(), et->application());
    j->set_name("RemoteCopySend");
    j->set_id(cm->job_id());
    j->set_receive_job_id(cm->receive_job_id());
    j->set_mega_rcr_job_id(cm->mega_rcr_job_id());
    j->set_to_worker_id(cm->to_worker_id());
    j->set_to_ip(cm->to_ip());
    j->set_to_port(cm->to_port());
    IDSet<physical_data_id_t> read_set;
    read_set.insert(cm->from_physical_data_id().elem());
    j->set_read_set(read_set);
    j->set_shadow_job_id(shadow_job_id);
    j->set_execution_template(et);
    j->set_template_generation_id(et->template_generation_id());

    ready_jobs->push_back((j));
    ++count;
  }

  return count;
}


size_t TemplateExtension::LoadReceiveJob(RemoteCopyReceiveJob **rcr,
                                         const WorkerDataExchanger::Event& e,
                                         ExecutionTemplate *et) {
  job_id_t shadow_job_id = compute_command_->job_id().elem();

  bool found = false;
  std::vector<boost::shared_ptr<RemoteCopyReceiveCommand> >::iterator iter =
    receive_commands_.begin();
  for (; iter != receive_commands_.end(); ++iter) {
    boost::shared_ptr<RemoteCopyReceiveCommand> cm = *iter;
    if (cm->job_id().elem() == e.receive_job_id_) {
      RemoteCopyReceiveJob * j = new RemoteCopyReceiveJob(et->application());
      j->set_name("RemoteCopyReceive");
      j->set_id(cm->job_id());
      IDSet<physical_data_id_t> write_set;
      write_set.insert(cm->to_physical_data_id().elem());
      j->set_write_set(write_set);
      j->set_shadow_job_id(shadow_job_id);
      j->set_execution_template(et);

      j->set_data_version(e.version_);
      j->set_serialized_data(e.ser_data_);

      *rcr = j;
      found = true;
      break;
    }
  }
  assert(found);

  return 1;
}



