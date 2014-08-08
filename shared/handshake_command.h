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
  * The handshake command is sent by a worker to the controller when
  * it first starts. It allows the controller to register the worker
  * and assign it a unique ID.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_SHARED_HANDSHAKE_COMMAND_H_
#define NIMBUS_SHARED_HANDSHAKE_COMMAND_H_


#include <string>
#include "shared/scheduler_command.h"

namespace nimbus {
class HandshakeCommand : public SchedulerCommand {
  public:
    HandshakeCommand();
    explicit HandshakeCommand(const ID<worker_id_t>& worker_id,
                              const std::string& ip,
                              const ID<port_t>& port);
    ~HandshakeCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& data);
    virtual bool Parse(const SchedulerPBuf& buf);
    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<worker_id_t> worker_id();
    std::string ip();
    ID<port_t> port();

  private:
    ID<worker_id_t> worker_id_;
    std::string ip_;
    ID<port_t> port_;

    bool ReadFromProtobuf(const HandshakePBuf& buf);
    bool WriteToProtobuf(HandshakePBuf* buf);
};


}  // namespace nimbus

#endif  // NIMBUS_SHARED_HANDSHAKE_COMMAND_H_
