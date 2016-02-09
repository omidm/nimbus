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
  * A command sent from the controller to a worker to save a physical data to
  * non-volatile memory. This command is used to create distributed checkpoint
  * in the system to rewind back to in case of failure.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_SAVE_DATA_COMMAND_H_
#define NIMBUS_SRC_SHARED_SAVE_DATA_COMMAND_H_


#include <string>
#include "src/shared/scheduler_command.h"

namespace nimbus {
class SaveDataCommand : public SchedulerCommand {
  public:
    SaveDataCommand();
    SaveDataCommand(const ID<job_id_t>& job_id,
                    const ID<physical_data_id_t>& from_physical_data_id,
                    const ID<checkpoint_id_t>& checkpoint_id,
                    const IDSet<job_id_t>& before);
    ~SaveDataCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& data);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();
    ID<job_id_t> job_id();
    ID<physical_data_id_t> from_physical_data_id();
    ID<checkpoint_id_t> checkpoint_id();
    IDSet<job_id_t> before_set();

  private:
    ID<job_id_t> job_id_;
    ID<physical_data_id_t> from_physical_data_id_;
    ID<checkpoint_id_t> checkpoint_id_;
    IDSet<job_id_t> before_set_;

    bool ReadFromProtobuf(const SaveDataPBuf& buf);
    bool WriteToProtobuf(SaveDataPBuf* buf);
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_SAVE_DATA_COMMAND_H_
