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

#ifndef NIMBUS_SRC_SHARED_TEMPLATE_EXTENSION_H_
#define NIMBUS_SRC_SHARED_TEMPLATE_EXTENSION_H_

#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <string>
#include <vector>
#include <list>
#include "src/shared/nimbus_types.h"
#include "src/shared/compute_job_command.h"
#include "src/shared/remote_copy_send_command.h"
#include "src/shared/remote_copy_receive_command.h"
#include "src/shared/worker_data_exchanger.h"
#include "src/worker/job.h"
#include "src/shared/dbg.h"
#include "src/shared/log.h"

namespace nimbus {

class ExecutionTemplate;

class TemplateExtension {
  public:
    TemplateExtension(
    const bool& migrate_out,
    boost::shared_ptr<ComputeJobCommand> compute_command,
    const std::vector<boost::shared_ptr<RemoteCopySendCommand> >& send_commands,
    const std::vector<boost::shared_ptr<RemoteCopyReceiveCommand> >& receive_commands);

    ~TemplateExtension();

    bool migrate_out() const;
    boost::shared_ptr<ComputeJobCommand> compute_command() const;
    std::vector<boost::shared_ptr<RemoteCopySendCommand> >* send_commands_p();
    std::vector<boost::shared_ptr<RemoteCopyReceiveCommand> >* receive_commands_p();

    enum State {
      BLOCK,
      RELEASE
    };

    size_t MarkJobDone(JobList *ready_jobs,
                       ExecutionTemplate *et,
                       State *state,
                       bool append);

    size_t LoadComputeJob(JobList *ready_jobs,
                          ExecutionTemplate *et,
                          bool append);

    size_t LoadSendJobs(JobList *ready_jobs,
                        ExecutionTemplate *et,
                        bool append);

    size_t LoadReceiveJob(RemoteCopyReceiveJob **rcr,
                          const WorkerDataExchanger::Event& e,
                          ExecutionTemplate *et);


  private:
    bool migrate_out_;
    boost::shared_ptr<ComputeJobCommand> compute_command_;
    std::vector<boost::shared_ptr<RemoteCopySendCommand> > send_commands_;
    std::vector<boost::shared_ptr<RemoteCopyReceiveCommand> > receive_commands_;

    size_t job_done_counter_;
};

typedef boost::unordered_map<worker_id_t, std::vector<TemplateExtension> > ExtensionsMap;

}  // namespace nimbus
#endif  // NIMBUS_SRC_SHARED_TEMPLATE_EXTENSION_H_
