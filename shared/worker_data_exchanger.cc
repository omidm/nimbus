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

#include <sstream> // NOLINT
#include "shared/worker_data_exchanger.h"

#define _NIMBUS_NO_NETWORK_LOG

using boost::asio::ip::tcp;

namespace nimbus {

WorkerDataExchanger::WorkerDataExchanger(port_t port_no)
  : listening_port_(port_no) {
  io_service_ = new boost::asio::io_service();
  acceptor_ = new tcp::acceptor(*io_service_,
      tcp::endpoint(tcp::v4(), listening_port_));

  std::ostringstream ss;
  ss << port_no;
  log_.set_file_name("log_wdx-" + ss.str());
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

/*
void WorkerDataExchanger::AddSendConnection(worker_id_t worker_id,
      WorkerDataExchangerConnection* connection) {
  boost::mutex::scoped_lock lock(send_connection_mutex_);
  send_connections_[worker_id] = connection;
}
*/

void WorkerDataExchanger::HandleRead(WorkerDataExchangerConnection* connection,
                                 const boost::system::error_code& error,
                                 size_t bytes_transferred) {
  if (error) {
    dbg(DBG_NET|DBG_ERROR, "Error %s.\n", error.message().c_str());
    return;
  }

  size_t read_len = 0;
  size_t read_data_len, read_header_len, remaining;
  size_t bytes_available = bytes_transferred + connection->existing_bytes();

  while (true) {
    read_data_len = 0;
    read_header_len = 0;
    remaining = bytes_available - read_len;

    if (connection->middle_of_header()) {
      read_header_len = ReadHeader(connection,
          connection->read_buffer() + read_len, remaining);
      read_len += read_header_len;
      remaining = bytes_available - read_len;
    }

    if (connection->middle_of_data()) {
      read_data_len = ReadData(connection,
          connection->read_buffer() + read_len, remaining);
      read_len += read_data_len;
      remaining = bytes_available - read_len;
    }

    if (((read_data_len + read_header_len) == 0) || (remaining == 0))
      break;
  }

  // This is the case when the buffer is full with incomplete message.
  // Then it cannot make progress anymore.
  if (remaining == connection->read_buffer_max_length()) {
    std::string err_msg = "ERROR: ";
    err_msg += "worker data exchanger buffer is full with incomplete data. ";
    err_msg += "You need to increase the WORKER_DATA_BUFSIZE ";
    err_msg += "in worker_data_exchanger_connection.cc.\n";
    dbg(DBG_ERROR, "%s\n.", err_msg.c_str());
    exit(-1);
  }

  char* buffer = connection->read_buffer();
  memmove(buffer, (buffer + read_len), remaining);
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
      data_version_t version = 0;
      buffer[i] = '\0';
      std::string input(buffer);
      ParseWorkerDataHeader(input, job_id, data_length, version);
      connection->set_middle_of_data(true);
      connection->set_middle_of_header(false);
      connection->set_job_id(job_id);
      connection->set_data_length(data_length);
      connection->set_data_version(version);
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
    AddSerializedData(connection->job_id(), ser_data, connection->data_version());
    return remaining;
  }
}

void WorkerDataExchanger::AddSerializedData(job_id_t job_id,
    SerializedData* ser_data, data_version_t version) {
#ifndef _NIMBUS_NO_NETWORK_LOG
  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff), "R %10.9lf j: %5.0lu s: %5.0lu",
      Log::GetRawTime(), job_id, ser_data->size());
  log_.log_WriteToFile(std::string(buff));
#endif
  boost::mutex::scoped_lock lock(data_map_mutex_);
  assert(data_map_.find(job_id) == data_map_.end());
  data_map_[job_id] = std::make_pair(ser_data, version);
  receive_events.push_back(job_id);
}

bool WorkerDataExchanger::GetReceiveEvent(job_id_t* job_id) {
  boost::mutex::scoped_lock lock(data_map_mutex_);
  if (receive_events.empty()) {
    return false;
  } else {
    *job_id = receive_events.front();
    receive_events.pop_front();
    return true;
  }
}

void WorkerDataExchanger::RemoveSerializedData(job_id_t job_id) {
  boost::mutex::scoped_lock lock(data_map_mutex_);
  data_map_.erase(job_id);
}

bool WorkerDataExchanger::AddContactInfo(worker_id_t worker_id,
      std::string ip_address, port_t port_no) {
  boost::mutex::scoped_lock lock(address_book_mutex_);
  address_book_[worker_id] = std::make_pair(ip_address, port_no);
  return true;
}

bool WorkerDataExchanger::ReceiveSerializedData(job_id_t job_id,
      SerializedData** ser_data, data_version_t& version) {
  int available;
  {
    boost::mutex::scoped_lock lock(data_map_mutex_);
    available = data_map_.count(job_id);
  }

  if (!available) {
    return false;
  } else {
    *ser_data = data_map_[job_id].first;
    version = data_map_[job_id].second;
    RemoveSerializedData(job_id);
    return true;
  }
}

bool WorkerDataExchanger::SendSerializedData(job_id_t job_id,
      worker_id_t worker_id, SerializedData& ser_data, data_version_t version) {
  boost::shared_array<char> buf = ser_data.data_ptr();
  size_t size = ser_data.size();
  boost::mutex::scoped_lock lock1(send_connection_mutex_);
  boost::mutex::scoped_lock lock2(address_book_mutex_);
  if (send_connections_.count(worker_id) == 0) {
    if (address_book_.count(worker_id) == 0) {
      std::cout << "ERROR: could not find the address info of the worker id: " <<
        worker_id << std::endl;
      return false;
    } else {
      std::pair<std::string, port_t> add = address_book_[worker_id];
      CreateNewSendConnection(worker_id, add.first, add.second);
    }
  }
  WorkerDataExchangerConnection* connection = send_connections_[worker_id];
  boost::system::error_code ignored_error;

  std::string header;
  std::ostringstream ss_j;
  ss_j << job_id;
  header += (ss_j.str() + " ");
  std::ostringstream ss_s;
  ss_s << size;
  header += (ss_s.str() + " ");
  std::ostringstream ss_v;
  ss_v << version;
  header += (ss_v.str() + ";");

#ifndef _NIMBUS_NO_NETWORK_LOG
  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff), "S %10.9lf j: %5.0lu s: %5.0lu",
      Log::GetRawTime(), job_id, size);
  log_.log_WriteToFile(std::string(buff));
#endif

  boost::asio::write(*(connection->socket()), boost::asio::buffer(header),
      boost::asio::transfer_all(), ignored_error);
  boost::asio::write(*(connection->socket()),
      boost::asio::buffer(buf.get(), size),
      boost::asio::transfer_all(), ignored_error);
  return true;
}

bool WorkerDataExchanger::CreateNewSendConnection(worker_id_t worker_id,
      std::string ip_address, port_t port_no) {
  std::cout << "Opening new worker-worker connection." << std::endl;
  tcp::resolver resolver(*io_service_);
  tcp::resolver::query query(ip_address,
      boost::to_string(port_no));
  tcp::resolver::iterator iterator = resolver.resolve(query);
  boost::system::error_code error =
    boost::asio::error::host_not_found;
  tcp::socket* socket = new tcp::socket(*io_service_);
  WorkerDataExchangerConnection* connection =
    new WorkerDataExchangerConnection(socket);
  socket->connect(*iterator, error);
  send_connections_[worker_id] = connection;
  return true;
}

void WorkerDataExchanger::HandleWrite(WorkerDataExchangerConnection* connection,
                           const boost::system::error_code& error,
                           size_t bytes_transferred) {
}

WorkerDataExchangerConnectionMap* WorkerDataExchanger::send_connections() {
  return &send_connections_;
}

WorkerDataExchangerConnectionList* WorkerDataExchanger::receive_connections() {
  return &receive_connections_;
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
    dbg(DBG_NET, "Copying command %s to user buffer.\n", command->ToNetworkData().c_str());
    storage->push_back(command);
  }
  return true;
}



void WorkerDataExchanger::SendCommand(SchedulerWorker* worker,
                                  SchedulerCommand* command) {
  WorkerDataExchangerConnection* connection = worker->connection();
  std::string msg = command->ToNetworkData() + ";";
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
        dbg(DBG_NET, "Adding command %s to queue.\n", command->ToNetworkData().c_str());
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
