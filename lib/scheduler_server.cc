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
                                      uint32_t maxCommands) {return true;}



void SchedulerServer::SendCommand(SchedulerWorker* worker,
                                  SchedulerCommand* command) {
  SchedulerServerConnection* connection = worker->connection();
  std::string msg = command->toString() + ";";
  boost::system::error_code ignored_error;
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

  std::list<std::string> stringList;
  int end = separateCommands(worker->buffer(), bytes_transferred,
                             stringList);

  std::list<std::string>::iterator listIterator = stringList.begin();
  for (; listIterator != stringList.end(); ++listIterator) {
    std::string sval = *listIterator;
  }
  // Go through worker buffer, parsing commands and putting
  // on commmand list

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
