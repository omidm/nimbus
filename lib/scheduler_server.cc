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
#include <string>
#include <sstream>
#include <iostream>  // NOLINT
#include "lib/dbg.h"


using boost::asio::ip::tcp;

SchedulerServer::SchedulerServer(ConnectionId port_no)
  : connection_port_(port_no) {}

SchedulerServer::~SchedulerServer() {
  boost::mutex::scoped_lock lock(map_mutex_);
  for (ConnectionMapIter iter = connections_.begin();
       iter != connections_.end();
       ++iter)   {
    SchedulerServerConnection* scon = iter->second;
    delete scon;
  }
  connections_.clear();
  delete acceptor_;
  delete io_service_;
}

bool SchedulerServer::initialize() {
  io_service_ = new boost::asio::io_service();
  acceptor_ = new tcp::acceptor(*io_service_,
                                tcp::endpoint(tcp::v4(), connection_port_));
  return true;
}

bool SchedulerServer::receiveCommands(SchedulerCommandList* storage,
                                      uint32_t maxCommands) {return true;}



void SchedulerServer::sendCommand(SchedulerServerConnection* connection,
                                  SchedulerCommand* command) {
  std::string msg = command->toString() + ";";
  boost::system::error_code ignored_error;
  boost::asio::write(*(connection->socket()), boost::asio::buffer(msg),
                     boost::asio::transfer_all(), ignored_error);
}

void SchedulerServer::sendCommands(SchedulerServerConnection* connection,
                                    SchedulerCommandList* commands) {
  SchedulerCommandList::iterator iter = commands->begin();
  for (; iter != commands->end(); ++iter) {
    SchedulerCommand* command = *iter;
    sendCommand(connection, command);
  }
}

uint32_t SchedulerServer::receiveMessages() {
    // Implementation currently just receives commands from the
  // first connection, so it can be tested. Should pull commands
  // from all connections.
  ConnectionMapIter iter = connections_.begin();
  if (iter == connections_.end()) {
    return false;  // No active connections
  }

  SchedulerServerConnection* connection = iter->second;
  boost::system::error_code ignored_error;
  int bytes_available = connection->socket()->available(ignored_error);

  boost::asio::streambuf::mutable_buffers_type bufs =
    connection->read_buffer()->prepare(bytes_available);
  std::size_t bytes_read = connection->socket()->receive(bufs);
  connection->read_buffer()->commit(bytes_read);

  std::string str(boost::asio::buffer_cast<char*>(bufs), bytes_read);
  command_id_t num = connection->command_num() + countOccurence(str, ";");
  connection->set_command_num(num);

  // Why is this a conditional?
  if (connection->command_num() > 0) {
    std::istream input(connection->read_buffer());
    std::string commandText;
    std::getline(input, commandText, ';');
    connection->set_command_num(connection->command_num() - 1);
    SchedulerCommand* command = new SchedulerCommand(commandText);
    command->set_worker_id(connection->worker_id());
    received_commands_.push_back(command);
    return true;
  } else {
    return false;
  }
}

void SchedulerServer::listenForNewConnections() {
  dbg(DBG_NET, "Scheduler server listening for new connections.\n");
  tcp::socket* socket = new tcp::socket(*io_service_);
  SchedulerServerConnection* server_connection =
    new SchedulerServerConnection(socket);
  acceptor_->async_accept(*server_connection->socket(),
                          boost::bind(&SchedulerServer::handleAccept,
                                      this,
                                      server_connection,
                                      boost::asio::placeholders::error));
}

void SchedulerServer::handleAccept(SchedulerServerConnection* connection,
                                   const boost::system::error_code& error) {
  if (!error) {
    dbg(DBG_NET, "Scheduler accepted new connection.\n");
    boost::mutex::scoped_lock lock(map_mutex_);
    connections_[connection->id()] = connection;
    listenForNewConnections();
  } else {
    delete connection;
  }
}

void SchedulerServer::run() {
  dbg(DBG_NET, "Running the scheduler networking server.\n");
  initialize();
  listenForNewConnections();
  // receiveMessages();
  io_service_->run();
}

// connection_subscription_thread_ = new boost::thread(
//    boost::bind(&SchedulerServer::listenForNewConnections, this));
//  }

ConnectionMap* SchedulerServer::connections() {
  return &connections_;
}

SchedulerServerConnection::SchedulerServerConnection(tcp::socket* sock)
  :socket_(sock) {
  static ConnectionId id_assigner = 0;
  id_ = ++id_assigner;
  read_buffer_ = new boost::asio::streambuf();
  command_num_ = 0;
}

SchedulerServerConnection::~SchedulerServerConnection() {
  // FIXME: not actually cleaning up listening thread.
}


boost::asio::streambuf* SchedulerServerConnection::read_buffer() {
  return read_buffer_;
}

tcp::socket* SchedulerServerConnection::socket() {
  return socket_;
}

int SchedulerServerConnection::command_num() {
  return command_num_;
}

void SchedulerServerConnection::set_command_num(int n) {
  command_num_ = n;
}

worker_id_t SchedulerServerConnection::worker_id() {
  return worker_id_;
}

void SchedulerServerConnection::set_worker_id(worker_id_t w) {
  worker_id_ = w;
}

ConnectionId SchedulerServerConnection::id() {
  return id_;
}

