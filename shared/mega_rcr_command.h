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

#ifndef NIMBUS_SHARED_MEGA_RCR_COMMAND_H_
#define NIMBUS_SHARED_MEGA_RCR_COMMAND_H_


#include <string>
#include <vector>
#include "shared/scheduler_command.h"

namespace nimbus {
class MegaRCRCommand : public SchedulerCommand {
  public:
    MegaRCRCommand();
    MegaRCRCommand(const ID<job_id_t>& job_id,
                   const std::vector<job_id_t>& receive_job_ids,
                   const std::vector<physical_data_id_t>& to_physical_data_ids);

    ~MegaRCRCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& data);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();
    ID<job_id_t> job_id();
    std::vector<job_id_t> receive_job_ids();
    const std::vector<job_id_t>* receive_job_ids_p();
    std::vector<physical_data_id_t> to_physical_data_ids();
    const std::vector<physical_data_id_t>* to_physical_data_ids_p();

  private:
    ID<job_id_t> job_id_;
    std::vector<job_id_t> receive_job_ids_;
    std::vector<physical_data_id_t> to_physical_data_ids_;

    bool ReadFromProtobuf(const MegaRCRPBuf& buf);
    bool WriteToProtobuf(MegaRCRPBuf* buf);
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_MEGA_RCR_COMMAND_H_
