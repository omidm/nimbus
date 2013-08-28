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
  * The interface that the workers use to exchange data among each other.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/worker_data_exchanger.h"


using boost::asio::ip::tcp;

#define BUF_SIZE 102400
namespace nimbus {

WorkerDataExchanger::WorkerDataExchanger(port_t port_no)
  : listening_port_(port_no) {
  io_service_ = new boost::asio::io_service();
  acceptor_ = new tcp::acceptor(*io_service_,
      tcp::endpoint(tcp::v4(), listening_port_));
}

WorkerDataExchanger::~WorkerDataExchanger() {
  {
    boost::mutex::scoped_lock lock(send_connection_mutex_);
    WorkerDataExchangerConnectionMap::iterator iter = send_connections_.begin();
    for (; iter != send_connections_.end(); iter++)   {
      WorkerDataExchangerConnection* connection = iter->second;
      delete connection;
    }
    send_connections_.clear();
  }
  {
    boost::mutex::scoped_lock lock(receive_connection_mutex_);
    WorkerDataExchangerConnectionList::iterator iter = receive_connections_.begin();
    for (; iter != receive_connections_.end(); iter++)   {
      WorkerDataExchangerConnection* connection = *iter;
      delete connection;
    }
    receive_connections_.clear();
  }
  delete acceptor_;
  delete io_service_;
}

void WorkerDataExchanger::Run() {
  dbg(DBG_NET, "Running the worker data exchanger.\n");
  ListenForNewConnections();
  io_service_->run();
}

void WorkerDataExchanger::ListenForNewConnections() {
  dbg(DBG_NET, "Worker data exchanger listening for new connections.\n");
  tcp::socket* socket = new tcp::socket(*io_service_);
  WorkerDataExchangerConnection* connection =
    new WorkerDataExchangerConnection(socket);
  acceptor_->async_accept(*connection->socket(),
      boost::bind(&WorkerDataExchanger::HandleAccept,
        this, connection, boost::asio::placeholders::error));
}

void WorkerDataExchanger::HandleAccept(WorkerDataExchangerConnection* connection,
                                   const boost::system::error_code& error) {
  if (!error) {
    dbg(DBG_NET, "Worker accepted new connection.\n");
    AddReceiveConnection(connection);
    ListenForNewConnections();
    boost::asio::async_read(*(connection->socket()),
        boost::asio::buffer(connection->read_buffer(),
          connection->read_buffer_max_length()),
        boost::asio::transfer_at_least(1),
        boost::bind(&WorkerDataExchanger::HandleRead,
          this,
          connection,
          boost::asio::placeholders::error,
          boost::asio::placeholders::bytes_transferred));
  } else {
    delete connection;
  }
}

void WorkerDataExchanger::AddReceiveConnection(WorkerDataExchangerConnection* connection) { //NOLINT
  boost::mutex::scoped_lock lock(receive_connection_mutex_);
  receive_connections_.push_back(connection);
}

void WorkerDataExchanger::AddSendConnection(worker_id_t worker_id,
      WorkerDataExchangerConnection* connection) {
  boost::mutex::scoped_lock lock(send_connection_mutex_);
  send_connections_[worker_id] = connection;
}

void WorkerDataExchanger::HandleRead(WorkerDataExchangerConnection* connection,
                                 const boost::system::error_code& error,
                                 size_t bytes_transferred) {
  if (error) {
    dbg(DBG_NET|DBG_ERROR, "Error %s.\n", error.message().c_str());
    return;
  }

  size_t read_data_len = 0;
  size_t read_header_len = 0;
  size_t bytes_available = bytes_transferred + connection->existing_bytes();

  if (connection->middle_of_header()) {
    read_header_len = ReadHeader(connection, connection->read_buffer(), bytes_available);
  }
  if (connection->middle_of_data()) {
    read_data_len = ReadData(connection, connection->read_buffer() + read_header_len,
        bytes_available - read_header_len);
  }
  size_t read_len = read_data_len + read_header_len;
  size_t remaining = bytes_available - read_len;
  char* buffer = connection->read_buffer();
  memcpy(buffer, (buffer + read_len), remaining);
  connection->set_existing_bytes(remaining);

  char* read_start_ptr = buffer + remaining;
  size_t read_buffer_length = connection->read_buffer_max_length() - remaining;

  boost::asio::async_read(*(connection->socket()),
      boost::asio::buffer(read_start_ptr,
        read_buffer_length),
      boost::asio::transfer_at_least(1),
      boost::bind(&WorkerDataExchanger::HandleRead,
        this,
        connection,
        boost::asio::placeholders::error,
        boost::asio::placeholders::bytes_transferred));
}

size_t WorkerDataExchanger::ReadHeader(WorkerDataExchangerConnection* connection,
      char* buffer, size_t size) {
  size_t i = 0;
  for (; i < size; i++) {
    if (buffer[i] == ';') {
      job_id_t job_id = 0;
      size_t data_length = 0;
      buffer[i] = '\0';
      std::string input(buffer);
      // TODO(omidm) parse the input.
      connection->set_middle_of_data(true);
      connection->set_middle_of_header(false);
      connection->set_job_id(job_id);
      connection->set_data_length(data_length);
      connection->AllocateData(data_length);
      return (i + 1);
    }
  }
  return 0;
}

size_t WorkerDataExchanger::ReadData(WorkerDataExchangerConnection* connection,
      char* buffer, size_t size) {
  size_t remaining = connection->remaining_data_length();
  if (size < remaining) {
    connection->AppendData(buffer, size);
    return size;
  } else {
    connection->AppendData(buffer, remaining);
    connection->set_middle_of_data(false);
    connection->set_middle_of_header(true);
    SerializedData* ser_data =
      new SerializedData(connection->data_ptr(), connection->data_length());
    data_map_[connection->job_id()] = ser_data;
    return remaining;
  }
}

/*

bool WorkerDataExchanger::ReceiveCommands(SchedulerCommandList* storage,
                                      uint32_t maxCommands) {
  boost::mutex::scoped_lock(command_mutex_);
  uint32_t pending = received_commands_.size();
  if (pending == 0) {
    return false;
  } else if (pending < maxCommands) {
    maxCommands = pending;
  }
  storage->clear();
  for (uint32_t i = 0; i < maxCommands; i++) {
    SchedulerCommand* command = received_commands_.front();
    received_commands_.pop_front();
    dbg(DBG_NET, "Copying command %s to user buffer.\n", command->toString().c_str());
    storage->push_back(command);
  }
  return true;
}



void WorkerDataExchanger::SendCommand(SchedulerWorker* worker,
                                  SchedulerCommand* command) {
  WorkerDataExchangerConnection* connection = worker->connection();
  std::string msg = command->toString() + ";";
  dbg(DBG_NET, "Sending command %s.\n", msg.c_str());
  boost::system::error_code ignored_error;
  // Why are we IGNORING ERRORS!??!?!?
  boost::asio::write(*(connection->socket()), boost::asio::buffer(msg),
                     boost::asio::transfer_all(), ignored_error);
}

void WorkerDataExchanger::SendCommands(SchedulerWorker* worker,
                                    SchedulerCommandList* commands) {
  SchedulerCommandList::iterator iter = commands->begin();
  for (; iter != commands->end(); ++iter) {
    SchedulerCommand* command = *iter;
    SendCommand(worker, command);
  }
}

using boost::tokenizer;
using boost::char_separator;

int WorkerDataExchanger::EnqueueCommands(char* buffer, size_t size) {
  buffer[size] = '\0';
  dbg(DBG_NET, "Read string %s from worker.\n", buffer);

  char* start_pointer = buffer;
  for (size_t i = 0; i < size; i++) {
    // When we find a semicolon, replace it with a string terminator \0.
    // Then when we pass start_pointer to the constructor it will terminate.
    if (buffer[i] == ';') {
      buffer[i] = '\0';
      std::string input(start_pointer);

      SchedulerCommand* command;
      if (SchedulerCommand::GenerateSchedulerCommandChild(
          input, worker_command_set_, command)) {
        dbg(DBG_NET, "Adding command %s to queue.\n", command->toString().c_str());
        boost::mutex::scoped_lock lock(command_mutex_);
        received_commands_.push_back(command);
      } else {
        dbg(DBG_NET, "Ignored unknown command: %s.\n", input.c_str());
      }
      // Next string starts after the semicolon
      start_pointer = buffer + i + 1;
    }
  }
  // We've read this many bytes successfully into
  // commands
  return start_pointer - buffer;
}


void WorkerDataExchanger::HandleWrite(SchedulerWorker* worker,
                                  const boost::system::error_code& error,
                                  size_t bytes_transferred) {}


SchedulerWorkerList* WorkerDataExchanger::workers() {
  return &workers_;
}

void WorkerDataExchanger::set_worker_command_set(CommandSet* cms) {
  worker_command_set_ = cms;
}

void WorkerDataExchanger::MarkWorkerDead(SchedulerWorker* worker) {
  worker->MarkDead();
}
*/

}  // namespace nimbus
