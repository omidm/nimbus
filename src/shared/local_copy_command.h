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
  * A command sent from the controller to a worker to copy between two
  * physical objects on a worker. A common use is when data for
  * timestep t is computed with data from timestep t, or when one job
  * needs to read data from timestep t while another writes it for
  * timestep t+1. Copies between workers are different and called
  * remote copies.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_LOCAL_COPY_COMMAND_H_
#define NIMBUS_SRC_SHARED_LOCAL_COPY_COMMAND_H_


#include <string>
#include "src/shared/scheduler_command.h"

namespace nimbus {
class LocalCopyCommand : public SchedulerCommand {
  public:
    LocalCopyCommand();
    LocalCopyCommand(const ID<job_id_t>& job_id,
                     const ID<physical_data_id_t>& from_physical_data_id,
                     const ID<physical_data_id_t>& to_physical_data_id,
                     const IDSet<job_id_t>& before);
    ~LocalCopyCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& data);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();
    ID<job_id_t> job_id();
    ID<physical_data_id_t> from_physical_data_id();
    ID<physical_data_id_t> to_physical_data_id();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t>* before_set_p();

  private:
    ID<job_id_t> job_id_;
    ID<physical_data_id_t> from_physical_data_id_;
    ID<physical_data_id_t> to_physical_data_id_;
    IDSet<job_id_t> before_set_;

    bool ReadFromProtobuf(const LocalCopyPBuf& buf);
    bool WriteToProtobuf(LocalCopyPBuf* buf);
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_LOCAL_COPY_COMMAND_H_
