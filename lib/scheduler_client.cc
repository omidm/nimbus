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
  * Client (worker) side interface of the Nimbus scheduler protocol.
  *
  * This interface is a simple message passing interface to the scheduler.
  * It's responsible for transforming commands between object and wire
  * formats. On the reception side, this means instantiating the proper
  * command subclass using a factory.

  * Calls to the network should be non-blocking. In the case of
  * sending, this means we need to have error return values for
  * whether the send was successful. It should also properly handle no
  * messages to receive through non-blocking read calls.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "lib/scheduler_client.h"
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <string>
#include <sstream>
#include <iostream>  // NOLINT

using boost::asio::ip::tcp;


SchedulerClient::SchedulerClient(ConnectionId port_no)
  : connection_port_no(port_no) {
    io_service = new boost::asio::io_service();
    socket = new tcp::socket(*io_service);
}

SchedulerClient::~SchedulerClient() {}

SchedulerCommand* SchedulerClient::receiveCommand() {
  boost::asio::streambuf response;
  boost::asio::read_until(*socket, response, ";");

  std::istream is(&response);
  std::string msg;
  std::getline(is, msg);

  std::cout << "\nReceived msg: " << msg << "\n";
  SchedulerCommand* com = new SchedulerCommand(msg);
  return com;
}

void SchedulerClient::sendCommand(SchedulerCommand* command) {
  std::string msg = command->toString();
  boost::system::error_code ignored_error;
  boost::asio::write(*socket, boost::asio::buffer(msg),
      boost::asio::transfer_all(), ignored_error);
}

void SchedulerClient::createNewConnections() {
  std::cout << "Opening connections." << std::endl;
  tcp::resolver resolver(*io_service);
  tcp::resolver::query query("127.0.0.1", boost::to_string(connection_port_no));
  tcp::resolver::iterator iterator = resolver.resolve(query);
  boost::system::error_code error = boost::asio::error::host_not_found;
  socket->connect(*iterator, error);
}

void SchedulerClient::run() {
  createNewConnections();
}


