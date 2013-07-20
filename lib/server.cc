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

#include "server.h"
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <string>
#include <sstream>
#include <iostream>

using boost::asio::ip::tcp;

Server::Server(unsigned int _connection_port_no, Scheduler* sch)
: scheduler(sch),
  connection_subscription_thread(NULL),
  connection_port_no(_connection_port_no) {}

Server::~Server() {
  boost::mutex::scoped_lock lock(map_mutex);
  for (ConnectionMapIter iter = connections.begin();
       iter != connections.end();
       ++iter)   {
    ServerConn* scon = iter->second;
    delete scon;
  }
  connections.clear();
}

void Server::receive_msg(const std::string& msg, ServerConn* conn) {
  // FIXME: memory leak for killing and ending thread listening for
  // new connections.

  std::cout << "\nReceived msg " << msg << "\n";

  boost::mutex::scoped_lock lock(map_mutex);
  for (ConnectionMapIter iter = connections.begin();
       iter != connections.end();
       ++iter)   {
    if (iter->first != conn->get_id()) {
      iter->second->send_msg(msg);
    }
  }
}

void Server::listen_for_new_connections() {
  boost::asio::io_service io_service;
  tcp::acceptor acceptor(
      io_service, tcp::endpoint(tcp::v4(), connection_port_no));

  while (true) {
    tcp::socket* socket = new tcp::socket(io_service);
    acceptor.accept(*socket);
    {
      std::cout << "\nCreating new connection\n";
      boost::mutex::scoped_lock lock(map_mutex);
      ServerConn* sc = new ServerConn(this, socket);
      connections[sc->get_id()] = sc;
      sc->start_listening();
    }
  }
}

void Server::run() {
  connection_subscription_thread = new boost::thread(
      boost::bind(&Server::listen_for_new_connections, this));

  connection_subscription_thread->join();
}


ServerConn::ServerConn(Server* s, tcp::socket* sock): server(s), socket(sock) {
  static ConnectionId id_assigner = 0;
  id = id_assigner++;
}


void ServerConn::listen_for_msgs() {
  while (true)   {
    boost::asio::streambuf response;
    boost::asio::read_until(*socket, response, ";");

    std::istream is(&response);
    std::string msg;
    std::getline(is, msg);
    server->receive_msg(msg, this);
  }
}

void ServerConn::start_listening() {
  listening_thread = new boost::thread(
      boost::bind(&ServerConn::listen_for_msgs, this));
}

void ServerConn::send_msg(const std::string& msg) {
  boost::system::error_code ignored_error;
  boost::asio::write(*socket, boost::asio::buffer(msg),
      boost::asio::transfer_all(), ignored_error);
}

ServerConn::~ServerConn() {
  // FIXME: not actually cleaning up listening thread.
}



