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
  * A SpawnJobGraphCommand is a message sent from a worker to the
  * controller to spawn an instance of the job graph template.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_SPAWN_JOB_GRAPH_COMMAND_H_
#define NIMBUS_SRC_SHARED_SPAWN_JOB_GRAPH_COMMAND_H_


#include <string>
#include <vector>
#include "src/shared/scheduler_command.h"
#include "src/shared/geometric_region.h"
#include "src/shared/protobuf_compiled/commands.pb.h"

namespace nimbus {
class SpawnJobGraphCommand : public SchedulerCommand {
  public:
    SpawnJobGraphCommand();

    SpawnJobGraphCommand(const std::string& job_graph_name,
                         const std::vector<job_id_t>& inner_job_ids,
                         const std::vector<job_id_t>& outer_job_ids,
                         const std::vector<Parameter>& parameters,
                         const ID<job_id_t>& parent_job_id);

    ~SpawnJobGraphCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& param_segment);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();
    std::string job_graph_name();
    std::vector<job_id_t> inner_job_ids();
    std::vector<job_id_t> outer_job_ids();
    std::vector<Parameter> parameters();
    ID<job_id_t> parent_job_id();

  private:
    std::string job_graph_name_;
    std::vector<job_id_t> inner_job_ids_;
    std::vector<job_id_t> outer_job_ids_;
    std::vector<Parameter> parameters_;
    ID<job_id_t> parent_job_id_;

    bool ReadFromProtobuf(const SubmitJobGraphPBuf& buf);
    bool WriteToProtobuf(SubmitJobGraphPBuf* buf);
};



}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_SPAWN_JOB_GRAPH_COMMAND_H_
