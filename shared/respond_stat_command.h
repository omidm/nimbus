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
  * worker responds to controller query for its stats with a
  * RespondStatCommand, reporting idle, block and run time.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_RESPOND_STAT_COMMAND_H_
#define NIMBUS_SHARED_RESPOND_STAT_COMMAND_H_


#include <list>
#include <string>
#include "shared/scheduler_command.h"

namespace nimbus {
class RespondStatCommand : public SchedulerCommand {
  public:
    RespondStatCommand();
    RespondStatCommand(const counter_t& query_id,
                       const worker_id_t& worker_id,
                       const int64_t& run_time,
                       const int64_t& block_time,
                       const int64_t& idle_time);
    ~RespondStatCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& data);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string ToNetworkData();
    virtual std::string ToString();
    counter_t query_id();
    worker_id_t worker_id();
    int64_t run_time();
    int64_t block_time();
    int64_t idle_time();

    void set_final(bool flag);

  private:
    counter_t query_id_;
    worker_id_t worker_id_;
    int64_t run_time_;
    int64_t block_time_;
    int64_t idle_time_;

    bool ReadFromProtobuf(const RespondStatPBuf& buf);
    bool WriteToProtobuf(RespondStatPBuf* buf);
};

typedef std::list<RespondStatCommand*> RespondStatCommandList;

}  // namespace nimbus

#endif  // NIMBUS_SHARED_RESPOND_STAT_COMMAND_H_
