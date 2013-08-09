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
  */

#include "lib/scheduler_server.h"
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <string>
#include <sstream>
#include <iostream>  // NOLINT

using boost::asio::ip::tcp;

SchedulerServer::SchedulerServer(ConnectionId port_no)
  : connection_port_(port_no) {
  io_service_ = new boost::asio::io_service();
}

SchedulerServer::~SchedulerServer() {
  boost::mutex::scoped_lock lock(map_mutex_);
  for (ConnectionMapIter iter = connections_.begin();
       iter != connections_.end();
       ++iter)   {
    SchedulerServerConnection* scon = iter->second;
    delete scon;
  }
  connections_.clear();
}


SchedulerCommand* SchedulerServer::receiveCommand(SchedulerServerConnection* conn) { // NOLINT
  // boost::asio::read_until(*(conn->socket), *(conn->read_buffer), ';');
  // std::streamsize size = conn->read_buffer->in_avail();
  // boost::array<char, 128> buf;
  // boost::asio::read(*(conn->socket), boost::asio::buffer(buf, b));

  boost::system::error_code ignored_error;
  int bytes_available =conn->socket()->available(ignored_error);

  boost::asio::streambuf::mutable_buffers_type bufs =
    conn->read_buffer()->prepare(bytes_available);
  std::size_t bytes_read = conn->socket()->receive(bufs);
  conn->read_buffer()->commit(bytes_read);

  std::string str(boost::asio::buffer_cast<char*>(bufs), bytes_read);
  conn->set_command_num(conn->command_num() + countOccurence(str, ";"));

  if (conn->command_num() > 0) {
    std::istream input(conn->read_buffer());
    std::string command;
    std::getline(input, command, ';');
    conn->set_command_num(conn->command_num() - 1);
    SchedulerCommand* com = new SchedulerCommand(command);
    return com;
  } else {
    return new SchedulerCommand("no-command");
  }
}


void SchedulerServer::sendCommand(SchedulerServerConnection* conn,
    SchedulerCommand* command) {
  std::string msg = command->toString() + ";";
  boost::system::error_code ignored_error;
  boost::asio::write(*(conn->socket()), boost::asio::buffer(msg),
      boost::asio::transfer_all(), ignored_error);
}


void SchedulerServer::listenForNewConnections() {
  tcp::acceptor acceptor(
      *io_service_, tcp::endpoint(tcp::v4(), connection_port_));

  while (true) {
    tcp::socket* socket = new tcp::socket(*io_service_);
    acceptor.accept(*socket);
    {
      std::cout << "\nCreating new connection\n";
      boost::mutex::scoped_lock lock(map_mutex_);
      SchedulerServerConnection* sc =
        new SchedulerServerConnection(socket);
      connections_[sc->id()] = sc;
    }
  }
}

void SchedulerServer::run() {
  connection_subscription_thread_ = new boost::thread(
      boost::bind(&SchedulerServer::listenForNewConnections, this));
}


SchedulerServerConnection::SchedulerServerConnection(tcp::socket* sock)
  :socket_(sock) {
  static ConnectionId id_assigner = 0;
  id_ = id_assigner++;
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

ConnectionId SchedulerServerConnection::id() {
  return id_;
}

