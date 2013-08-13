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

#include "lib/scheduler_server.h"
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include <string>
#include <sstream>
#include <iostream>  // NOLINT
#include "lib/dbg.h"


using boost::asio::ip::tcp;

#define BUF_SIZE 102400
namespace nimbus {

SchedulerServer::SchedulerServer(ConnectionId port_no)
  : connection_port_(port_no) {}

SchedulerServer::~SchedulerServer() {
  {
    boost::mutex::scoped_lock lock(worker_mutex_);
    for (SchedulerWorkerList::iterator iter = workers_.begin();
         iter != workers_.end();
         ++iter)   {
      SchedulerWorker* worker = *iter;
      worker->MarkDead();
      delete worker;
    }
    workers_.clear();
  }
  delete acceptor_;
  delete io_service_;
}

bool SchedulerServer::Initialize() {
  io_service_ = new boost::asio::io_service();
  acceptor_ = new tcp::acceptor(*io_service_,
                                tcp::endpoint(tcp::v4(), connection_port_));
  return true;
}

bool SchedulerServer::ReceiveCommands(SchedulerCommandList* storage,
                                      uint32_t maxCommands) {
  boost::mutex::scoped_lock(command_mutex_);
  uint32_t pending = received_commands_.size();
  if (pending == 0) {
    return false;
  } else if (pending < maxCommands) {
    maxCommands = pending;
  }
  for (uint32_t i = 0; i < maxCommands; i++) {
    SchedulerCommand* command = received_commands_.front();
    received_commands_.pop_front();
    dbg(DBG_NET, "Copying command %s to user buffer.\n", command->toString().c_str());
    storage->push_back(command);
  }
  return true;
}



void SchedulerServer::SendCommand(SchedulerWorker* worker,
                                  SchedulerCommand* command) {
  SchedulerServerConnection* connection = worker->connection();
  std::string msg = command->toString() + ";";
  dbg(DBG_NET, "Sending command %s.\n", msg.c_str());
  boost::system::error_code ignored_error;
  // Why are we IGNORING ERRORS!??!?!?
  boost::asio::write(*(connection->socket()), boost::asio::buffer(msg),
                     boost::asio::transfer_all(), ignored_error);
}

void SchedulerServer::SendCommands(SchedulerWorker* worker,
                                    SchedulerCommandList* commands) {
  SchedulerCommandList::iterator iter = commands->begin();
  for (; iter != commands->end(); ++iter) {
    SchedulerCommand* command = *iter;
    SendCommand(worker, command);
  }
}

void SchedulerServer::ListenForNewConnections() {
  dbg(DBG_NET, "Scheduler server listening for new connections.\n");
  tcp::socket* socket = new tcp::socket(*io_service_);
  SchedulerServerConnection* server_connection =
    new SchedulerServerConnection(socket);
  acceptor_->async_accept(*server_connection->socket(),
                          boost::bind(&SchedulerServer::HandleAccept,
                                      this,
                                      server_connection,
                                      boost::asio::placeholders::error));
}

void SchedulerServer::HandleAccept(SchedulerServerConnection* connection,
                                   const boost::system::error_code& error) {
  if (!error) {
    dbg(DBG_NET, "Scheduler accepted new connection.\n");
    SchedulerWorker* worker =  AddWorker(connection);
    ListenForNewConnections();
    boost::asio::async_read(*(worker->connection()->socket()),
                            boost::asio::buffer(worker->read_buffer(),
                                                worker->read_buffer_length()),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerServer::HandleRead,
                                        this,
                                        worker,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  } else {
    delete connection;
  }
}

using boost::tokenizer;
using boost::char_separator;

  /** Reads commands from buffer with size bytes. Puts commands into
   * internal command list for later reads. Returns the number of
   * bytes read. */
int SchedulerServer::EnqueueCommands(char* buffer, size_t size) {
  buffer[size] = '\0';
  dbg(DBG_NET, "Read string %s from worker.\n", buffer);

  char* start_pointer = buffer;
  for (size_t i = 0; i < size; i++) {
    // When we find a semicolon, replace it with a string terminator \0.
    // Then when we pass start_pointer to the constructor it will terminate.
    if (buffer[i] == ';') {
      buffer[i] = '\0';
      std::string input(start_pointer);
      SchedulerCommand* command = new SchedulerCommand(input);
      {
        dbg(DBG_NET, "Adding command %s to queue.\n", command->toString().c_str());
        boost::mutex::scoped_lock lock(command_mutex_);
        received_commands_.push_back(command);
      }
      // Next string starts after the semicolon
      start_pointer = buffer + i + 1;
    }
  }
  // We've read this many bytes successfully into
  // commands
  dbg(DBG_NET, "Scheduler now has %i entries in command queue.\n", received_commands_.size());
  return start_pointer - buffer;
}

void SchedulerServer::HandleRead(SchedulerWorker* worker,
                                 const boost::system::error_code& error,
                                 size_t bytes_transferred) {
  if (error) {
    dbg(DBG_NET|DBG_ERROR,
        "Error %s receiving %i bytes from worker %i.\n",
        error.message().c_str(), bytes_transferred, worker->worker_id());
    return;
  }

  dbg(DBG_NET, "Scheduler received %i bytes from worker %i.\n",
      bytes_transferred, worker->worker_id());
  int real_length = bytes_transferred + worker->existing_bytes();
  int len = EnqueueCommands(worker->read_buffer(),
                           real_length);
  int remaining = (real_length - len);

  // This is for the case when the string buffer had an incomplete
  // command at its end. Copy the fragement of the command to the beginning
  // of the buffer, mark how many bytes are valid with existing_bytes.
  if (remaining > 0) {
    char* buffer = worker->read_buffer();
    memcpy(buffer, (buffer + len), remaining);
    char* read_start_ptr = buffer + remaining;
    int read_buffer_length = worker->read_buffer_length() - remaining;

    worker->set_existing_bytes(remaining);
    boost::asio::async_read(*(worker->connection()->socket()),
                            boost::asio::buffer(read_start_ptr,
                                                read_buffer_length),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerServer::HandleRead,
                                        this,
                                        worker,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  } else {
    worker->set_existing_bytes(0);
    boost::asio::async_read(*(worker->connection()->socket()),
                            boost::asio::buffer(worker->read_buffer(),
                                                worker->read_buffer_length()),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerServer::HandleRead,
                                        this,
                                        worker,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  }
}


void SchedulerServer::HandleWrite(SchedulerWorker* worker,
                                  const boost::system::error_code& error,
                                  size_t bytes_transferred) {}


void SchedulerServer::Run() {
  dbg(DBG_NET, "Running the scheduler networking server.\n");
  Initialize();
  ListenForNewConnections();
  io_service_->run();
}

SchedulerWorkerList* SchedulerServer::workers() {
  return &workers_;
}

SchedulerWorker* SchedulerServer::AddWorker(SchedulerServerConnection* connection) { //NOLINT
  boost::mutex::scoped_lock lock(worker_mutex_);
  static worker_id_t workerIdentifier = 0;
  workerIdentifier++;
  SchedulerWorker* worker = new SchedulerWorker(workerIdentifier,
                                                connection, NULL);
  workers_.push_back(worker);
  return worker;
}

void SchedulerServer::MarkWorkerDead(SchedulerWorker* worker) {
  worker->MarkDead();
}

}  // namespace nimbus
