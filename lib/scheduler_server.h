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
  * Server side of the Nimbus scheduler protocol.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */


#ifndef NIMBUS_LIB_SCHEDULER_SERVER_H_
#define NIMBUS_LIB_SCHEDULER_SERVER_H_

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <string>
#include <map>
#include "lib/nimbus.h"
#include "lib/parser.h"
#include "lib/scheduler_worker.h"
#include "lib/scheduler_server_connection.h"
#include "lib/scheduler_command.h"

namespace nimbus {

typedef uint32_t ConnectionId;

using boost::asio::ip::tcp;

class SchedulerServer {
 public:
  explicit SchedulerServer(ConnectionId port_no);
  virtual ~SchedulerServer();

  virtual bool Initialize();
  virtual void Run();
  virtual bool ReceiveCommands(SchedulerCommandList* storage,
                               uint32_t maxCommands);
  virtual void SendCommand(SchedulerWorker* connection,
                           SchedulerCommand* command);
  virtual void SendCommands(SchedulerWorker* connection,
                            SchedulerCommandList* commands);

  SchedulerWorkerList* workers();

 private:
  ConnectionId connection_port_;
  boost::mutex command_mutex_;
  SchedulerCommandList received_commands_;
  boost::mutex worker_mutex_;
  SchedulerWorkerList workers_;

  boost::asio::io_service* io_service_;
  tcp::acceptor* acceptor_;

  void ListenForNewConnections();
  void HandleAccept(SchedulerServerConnection* connection,
                    const boost::system::error_code& error);
  void HandleRead(SchedulerWorker* worker,
                  const boost::system::error_code& error,
                  size_t bytes_transferred);
  void HandleWrite(SchedulerWorker* worker,
                   const boost::system::error_code& error,
                   size_t bytes_transferred);

  SchedulerWorker* AddWorker(SchedulerServerConnection* connection);
  void MarkWorkerDead(SchedulerWorker* worker);
};

}  // namespace nimbus

#endif  // NIMBUS_LIB_SCHEDULER_SERVER_H_
