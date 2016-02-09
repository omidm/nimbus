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
  * A worker sends a define partition command to the controller to
  * name a geometric region. Later data objects can then refer to this
  * geometric region. Using partition identifiers allows the scheduler
  * to very easily determine if two defined datas cover the same
  * region and more precisely defines the geometric boundaries the
  * application is using.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_DEFINE_PARTITION_COMMAND_H_
#define NIMBUS_SRC_SHARED_DEFINE_PARTITION_COMMAND_H_


#include <string>
#include "src/shared/scheduler_command.h"

namespace nimbus {
class DefinePartitionCommand : public SchedulerCommand {
  public:
    DefinePartitionCommand();
    DefinePartitionCommand(const ID<partition_id_t>& partition_id,
                           const GeometricRegion& r);
    ~DefinePartitionCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& data);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();

    ID<partition_id_t> partition_id();
    const GeometricRegion* region();

  private:
    ID<partition_id_t> id_;
    GeometricRegion region_;

    bool ReadFromProtobuf(const DefinePartitionPBuf& buf);
    bool WriteToProtobuf(DefinePartitionPBuf* buf);
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_DEFINE_PARTITION_COMMAND_H_
