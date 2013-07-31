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
  : connection_port_no(port_no) {
  io_service = new boost::asio::io_service();
}

SchedulerServer::~SchedulerServer() {
  boost::mutex::scoped_lock lock(map_mutex);
  for (ConnectionMapIter iter = connections.begin();
       iter != connections.end();
       ++iter)   {
    SchedulerServerConnection* scon = iter->second;
    delete scon;
  }
  connections.clear();
}


SchedulerCommand* SchedulerServer::receiveCommand(SchedulerServerConnection* conn) { // NOLINT
  boost::asio::read_until(*(conn->socket), *(conn->read_buffer), ';');

  std::streamsize size = conn->read_buffer->in_avail();
  if (size > 0) {
    std::istream input(conn->read_buffer);
    std::string command;
    std::getline(input, command, ';');

    SchedulerCommand* com = new SchedulerCommand(command);

    // std::cout << "\nReceived " << command <<
    // " as " << com->toString() << "\n";

    return com;
  } else {
    return new SchedulerCommand("halt");
  }
}


void SchedulerServer::sendCommand(SchedulerServerConnection* conn,
    SchedulerCommand* command) {
  std::string msg = command->toString() + ";";
  boost::system::error_code ignored_error;
  boost::asio::write(*(conn->socket), boost::asio::buffer(msg),
      boost::asio::transfer_all(), ignored_error);
}


void SchedulerServer::listenForNewConnections() {
  tcp::acceptor acceptor(
      *io_service, tcp::endpoint(tcp::v4(), connection_port_no));

  while (true) {
    tcp::socket* socket = new tcp::socket(*io_service);
    acceptor.accept(*socket);
    {
      std::cout << "\nCreating new connection\n";
      boost::mutex::scoped_lock lock(map_mutex);
      SchedulerServerConnection* sc =
        new SchedulerServerConnection(socket);
      connections[sc->get_id()] = sc;
    }
  }
}

void SchedulerServer::run() {
  connection_subscription_thread = new boost::thread(
      boost::bind(&SchedulerServer::listenForNewConnections, this));
}


SchedulerServerConnection::SchedulerServerConnection(tcp::socket* sock)
  :socket(sock) {
  static ConnectionId id_assigner = 0;
  id = id_assigner++;
  read_buffer = new boost::asio::streambuf();
}

SchedulerServerConnection::~SchedulerServerConnection() {
  // FIXME: not actually cleaning up listening thread.
}



