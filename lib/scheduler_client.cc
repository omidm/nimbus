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

SchedulerClient::SchedulerClient(uint _connection_port_no, Worker* w)
: connection_port_no(_connection_port_no),
  worker(w) {
  io_service = new boost::asio::io_service();
  socket = new tcp::socket(*io_service);
}

SchedulerClient::~SchedulerClient() {}

/* 
 * This function never exits, and keeps listening to scheduler. Upon command
 * reception it will take the appropriate action. 
 */
void SchedulerClient::receiveCommand() {
  while (true) {
    boost::asio::streambuf response;
    boost::asio::read_until(*socket, response, ";");

    std::istream is(&response);
    std::string msg;
    std::getline(is, msg);

    std::cout << "\nReceived msg: " << msg << "\n";
    SchedulerCommand* com = new SchedulerCommand(msg);
    /*
     * will add the code to take the appropriate action based on com. - Omid
     */
  }
}

/* 
 * This function never exits, and sends the commands loaded in the command
 * transmission buffer to the scheduler.
 */
void SchedulerClient::sendCommand() {
  while (true) {
    /*
     * will add the code to sleep until the command buffer is loaded, then wake
     * up and send the commends. Forr now just the cin interface. - Omid
     */
    // std::string msg = command->toString();

    std::string msg;
    std::getline(std::cin, msg);
    boost::system::error_code ignored_error;
    boost::asio::write(*socket, boost::asio::buffer(msg),
        boost::asio::transfer_all(), ignored_error);
  }
}

void SchedulerClient::create_new_connections() {
  tcp::resolver resolver(*io_service);
  tcp::resolver::query query("127.0.0.1", boost::to_string(connection_port_no));
  tcp::resolver::iterator iterator = resolver.resolve(query);
  boost::system::error_code error = boost::asio::error::host_not_found;
  socket->connect(*iterator, error);
}

void SchedulerClient::run() {
  create_new_connections();
  sending_thread = new boost::thread(boost::bind(&SchedulerClient::sendCommand, this));  // NOLINT
  receiving_thread = new boost::thread(boost::bind(&SchedulerClient::receiveCommand, this)); // NOLINT

  // for now, have no other work to do, so just wait until the listening
  // thread terminates.
  sending_thread->join();
  receiving_thread->join();
}


