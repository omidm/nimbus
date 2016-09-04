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
  * The controller sends compute jobs to workers to invoke application
  * code.  They are translations of SpawnComputeJobs from workers to
  * the controller, binding logical objects to physical instances.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_COMPUTE_JOB_COMMAND_H_
#define NIMBUS_SRC_SHARED_COMPUTE_JOB_COMMAND_H_


#include <string>
#include "src/shared/scheduler_command.h"

namespace nimbus {
class ComputeJobCommand : public SchedulerCommand {
  public:
    ComputeJobCommand();
    ComputeJobCommand(const std::string& job_name,
                      const ID<job_id_t>& job_id,
                      const IDSet<physical_data_id_t>& read,
                      const IDSet<physical_data_id_t>& write,
                      const IDSet<physical_data_id_t>& scratch,
                      const IDSet<physical_data_id_t>& reduce,
                      const IDSet<job_id_t>& before,
                      const IDSet<job_id_t>& extra_dependency,
                      const IDSet<job_id_t>& after,
                      const ID<job_id_t>& future_job_id,
                      const bool& sterile,
                      const GeometricRegion& region,
                      const Parameter& params);
    ~ComputeJobCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& param_segment);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();
    std::string job_name();
    ID<job_id_t> job_id();
    IDSet<physical_data_id_t> read_set();
    IDSet<physical_data_id_t>* read_set_p();
    IDSet<physical_data_id_t> write_set();
    IDSet<physical_data_id_t>* write_set_p();
    IDSet<physical_data_id_t> scratch_set();
    IDSet<physical_data_id_t>* scratch_set_p();
    IDSet<physical_data_id_t> reduce_set();
    IDSet<physical_data_id_t>* reduce_set_p();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t>* before_set_p();
    IDSet<job_id_t> extra_dependency();
    IDSet<job_id_t>* extra_dependency_p();
    IDSet<job_id_t> after_set();
    IDSet<job_id_t>* after_set_p();
    ID<job_id_t> future_job_id();
    bool sterile();
    GeometricRegion region();
    Parameter params();

  private:
    std::string job_name_;
    ID<job_id_t> job_id_;
    IDSet<physical_data_id_t> read_set_;
    IDSet<physical_data_id_t> write_set_;
    IDSet<physical_data_id_t> scratch_set_;
    IDSet<physical_data_id_t> reduce_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> extra_dependency_;
    IDSet<job_id_t> after_set_;
    ID<job_id_t> future_job_id_;
    bool sterile_;
    GeometricRegion region_;
    Parameter params_;

  public:
    bool ReadFromProtobuf(const ExecuteComputeJobPBuf& buf);
    bool WriteToProtobuf(ExecuteComputeJobPBuf* buf);
};



}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_COMPUTE_JOB_COMMAND_H_
