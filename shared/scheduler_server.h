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


#ifndef NIMBUS_SHARED_SCHEDULER_SERVER_H_
#define NIMBUS_SHARED_SCHEDULER_SERVER_H_

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <string>
#include <map>
#include "shared/nimbus_types.h"
#include "shared/parser.h"
#include "shared/scheduler_command_include.h"
#include "shared/scheduler_server_connection.h"
#include "scheduler/scheduler_worker.h"

namespace nimbus {

typedef uint32_t ConnectionId;

using boost::asio::ip::tcp;

/** Networking services for a Nimbus scheduler to talk to workers.
 *  All receive and send calls are non-blocking. Currently
 *  maintains a list of workers asynchronously: future versions
 *  will provide callbacks to notify scheduler of updates to the
 *  worker list.
 */
class SchedulerServer {
 public:
  explicit SchedulerServer(port_t listening_port);
  virtual ~SchedulerServer();

  /** Start the server, does not return. A running server accepts
   *  connections, spools received messages into its message queue,
   *  and services send requests. */
  virtual void Run();

  /** Pull incoming commands from the received queue. Puts at most
   *  maxCommands into storage, returning true if it placed one or more.
   *  Returns false if no commands were placed in storage. */
  virtual bool ReceiveCommands(SchedulerCommandList* storage,
                               uint32_t maxCommands);

  /** Send command to destinationWorker. Returns immediately and
   *   processes the send asynchronously.*/

  virtual void SendCommand(SchedulerWorker* destinationWorker,
                           SchedulerCommand* command);

  /** Send commands to destinationWorker. Returns immediately and
   *  processes the send asynchronously. */
  virtual void SendCommands(SchedulerWorker* destinationWorker,
                            SchedulerCommandList* commands);


  /** Broadcast command to all workers. Returns immediately and processes
   *  the send asynchronously. */
  virtual void BroadcastCommand(SchedulerCommand* command);

  /** Broadcast commands to all workers. Returns immediately and processes
      the send asynchronously. */
  virtual void BroadcastCommands(SchedulerCommandList* commands);


  /** Returns the current list of workers. Note that SchedulerServer
   *  may modify this list or the workers on it. Accessing the list
   *  therefore involves locks for thread safety. */
  SchedulerWorkerList* workers();

  void set_worker_command_table(SchedulerCommand::PrototypeTable* cmt);

 private:
  port_t listening_port_;
  boost::mutex command_mutex_;
  SchedulerCommandList received_commands_;
  boost::mutex worker_mutex_;
  SchedulerWorkerList workers_;
  SchedulerCommand::PrototypeTable* worker_command_table_;

  boost::asio::io_service* io_service_;
  tcp::acceptor* acceptor_;

  /** Create server socket, set up networking and state. */
  virtual bool Initialize();

  /** Asynchronous call to accept incoming connections. */
  virtual void ListenForNewConnections();

  /** Callback for a new connection request. */
  virtual void HandleAccept(SchedulerServerConnection* connection,
                            const boost::system::error_code& error);

  /** Asynchronous callback to read data. */
  virtual void HandleRead(SchedulerWorker* worker,
                          const boost::system::error_code& error,
                          size_t bytes_transferred);

  /** Asynchronous calback for writing data. */
  virtual void HandleWrite(SchedulerWorker* worker,
                           const boost::system::error_code& error,
                           size_t bytes_transferred);

  /** Call to take a buffer from the network of size, parse them
      into commands and put them on the pending command queue.
      Return value is how many bytes were read/parsed, this can
      be less than size (for example, if the buffer end has an
      incomplete command. */
  virtual int EnqueueCommands(char* buffer, size_t size);

  /* Add a worker to the worker list. */
  virtual SchedulerWorker* AddWorker(SchedulerServerConnection* connection);

  /* Mark a worker dead. */
  virtual void MarkWorkerDead(SchedulerWorker* worker);
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_SCHEDULER_SERVER_H_
