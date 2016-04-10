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
  * This command batches job done ecvents and sends the batched command to
  * controller, it is used to signal the end of execution template by the
  * worker to the controller.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_MEGA_JOB_DONE_COMMAND_H_
#define NIMBUS_SRC_SHARED_MEGA_JOB_DONE_COMMAND_H_


#include <string>
#include <vector>
#include "src/shared/scheduler_command.h"

namespace nimbus {
class MegaJobDoneCommand : public SchedulerCommand {
  public:
    MegaJobDoneCommand();
    MegaJobDoneCommand(const std::vector<job_id_t>& job_ids,
                       const bool mark_stat);

    ~MegaJobDoneCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& data);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();
    std::vector<job_id_t> job_ids();
    const std::vector<job_id_t>* job_ids_p();
    bool mark_stat();

  private:
    std::vector<job_id_t> job_ids_;
    bool mark_stat_;

    bool ReadFromProtobuf(const MegaJobDonePBuf& buf);
    bool WriteToProtobuf(MegaJobDonePBuf* buf);
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_MEGA_JOB_DONE_COMMAND_H_
