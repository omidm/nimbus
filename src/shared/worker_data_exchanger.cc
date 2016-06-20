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
#include "src/shared/worker_data_exchanger.h"
#include "src/shared/fast_log.h"

#define _NIMBUS_NO_NETWORK_LOG

#define MAX_THREAD_NUM 1
#define MAX_HEADER_SIZE 100

// Note: you may have to increase the OS limits first.
// Look at the nimbus/scripts/configure_tcp.sh for help.
#ifndef __MACH__
#define DATA_EXCHANGER_TCP_SEND_BUF_SIZE 10485760  // 10MB
#define DATA_EXCHANGER_TCP_RECEIVE_BUF_SIZE 10485760  // 10MB
#endif

using boost::asio::ip::tcp;

namespace nimbus {

WorkerDataExchanger::WorkerDataExchanger(port_t port_no)
  : listening_port_(port_no) {
  io_service_ = new boost::asio::io_service();
  acceptor_ = new tcp::acceptor(*io_service_,
      tcp::endpoint(tcp::v4(), listening_port_));
  assert(MAX_THREAD_NUM >= 1);

#ifndef _NIMBUS_NO_NETWORK_LOG
  std::ostringstream ss;
  ss << port_no;
  log_.set_file_name("log_wdx-" + ss.str());
#endif
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

void WorkerDataExchanger::set_receive_event_mutex(boost::recursive_mutex *mutex) {
  receive_event_mutex_ = mutex;
}

void WorkerDataExchanger::set_receive_event_cond(boost::condition_variable_any *cond) {
  receive_event_cond_ = cond;
}

void WorkerDataExchanger::ThreadEntry() {
  timer::InitializeTimers();
  io_service_->run();
}

void WorkerDataExchanger::Run() {
  dbg(DBG_NET, "Running the worker data exchanger.\n");
  ListenForNewConnections();
}

void WorkerDataExchanger::ListenForNewConnections() {
  dbg(DBG_NET, "Worker data exchanger listening for new connections.\n");
  tcp::socket* socket = new tcp::socket(*io_service_);
  WorkerDataExchangerConnection* connection =
    new WorkerDataExchangerConnection(socket);
  acceptor_->async_accept(*connection->socket(),
      boost::bind(&WorkerDataExchanger::HandleAccept,
        this, connection, boost::asio::placeholders::error));
  if (threads_.size() < MAX_THREAD_NUM) {
    boost::thread *t = new boost::thread(
        boost::bind(&WorkerDataExchanger::ThreadEntry, this));
    threads_.push_back(t);
  }
  assert(threads_.size() >= 1);
}

void WorkerDataExchanger::HandleAccept(WorkerDataExchangerConnection* connection,
                                       const boost::system::error_code& error) {
  if (!error) {
    dbg(DBG_NET, "Worker accepted new connection.\n");
    AddReceiveConnection(connection);
    // Set the tcp send and receive buf size.
    // Note: you may have to increase the OS limits first.
    // Look at the nimbus/scripts/configure_tcp.sh for help.
#ifndef __MACH__
    boost::asio::socket_base::send_buffer_size s_option(DATA_EXCHANGER_TCP_SEND_BUF_SIZE);
    boost::asio::socket_base::receive_buffer_size r_option(DATA_EXCHANGER_TCP_RECEIVE_BUF_SIZE);
    connection->socket()->set_option(s_option);
    connection->socket()->set_option(r_option);
#endif
    // Turn of Nagle algorithm.
    boost::asio::ip::tcp::no_delay nd_option(TCP_NODELAY_OPTION);
    connection->socket()->set_option(nd_option);
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
  if (remaining > 0) {
    timer::StartTimer(timer::kDataExchangerLock);
    memmove(buffer, (buffer + read_len), remaining);
    timer::StopTimer(timer::kDataExchangerLock);
  }
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
                                       char* buffer,
                                       size_t size) {
  size_t i = 0;
  for (; i < size; i++) {
    if (buffer[i] == ';') {
      job_id_t receive_job_id = 0;
      job_id_t mega_rcr_job_id = 0;
      size_t data_length = 0;
      data_version_t version = 0;
      template_id_t template_generation_id = 0;

      char *next;
      receive_job_id = strtoll(buffer, &next, 10);
      mega_rcr_job_id = strtoll(next, &next, 10);
      data_length = strtoll(next, &next, 10);
      version = strtoll(next, &next, 10);
      template_generation_id = strtoll(next, &next, 10);
      assert((*next) == ';');

      connection->set_middle_of_data(true);
      connection->set_middle_of_header(false);
      connection->set_receive_job_id(receive_job_id);
      connection->set_mega_rcr_job_id(mega_rcr_job_id);
      connection->set_data_length(data_length);
      connection->set_data_version(version);
      connection->set_template_generation_id(template_generation_id);
      connection->AllocateData(data_length);
      return (i + 1);
    }
  }
  return 0;
}

size_t WorkerDataExchanger::ReadData(WorkerDataExchangerConnection* connection,
                                     char* buffer,
                                     size_t size) {
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
    AddSerializedData(connection->receive_job_id(),
                      connection->mega_rcr_job_id(),
                      ser_data,
                      connection->data_version(),
                      connection->template_generation_id());
    return remaining;
  }
}

void WorkerDataExchanger::AddSerializedData(job_id_t receive_job_id,
                                            job_id_t mega_rcr_job_id,
                                            SerializedData* ser_data,
                                            data_version_t version,
                                            template_id_t template_generation_id) {
#ifndef _NIMBUS_NO_NETWORK_LOG
  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff), "R %10.9lf j: %5.0lu s: %5.0lu",
      Log::GetRawTime(), receive_job_id, ser_data->size());
  log_.log_WriteToFile(std::string(buff));
#endif
  timer::StartTimer(timer::kDataExchangerLock);
  boost::unique_lock<boost::recursive_mutex> lock(*receive_event_mutex_);
  timer::StopTimer(timer::kDataExchangerLock);
  event_list_.push_back(Event(receive_job_id,
                              mega_rcr_job_id,
                              version,
                              ser_data,
                              template_generation_id));
  receive_event_cond_->notify_all();
}

size_t WorkerDataExchanger::PullReceiveEvents(EventList *events,
                                             size_t max_num) {
  events->clear();
  timer::StartTimer(timer::kDataExchangerLock);
  boost::unique_lock<boost::recursive_mutex> lock(*receive_event_mutex_);
  timer::StopTimer(timer::kDataExchangerLock);
  size_t count = 0;
  EventList::iterator iter = event_list_.begin();
  for (; (iter != event_list_.end()) && (count < max_num);) {
    events->push_back(*iter);
    event_list_.erase(iter++);
    ++count;
  }

  return count;
}

bool WorkerDataExchanger::AddContactInfo(worker_id_t worker_id,
                                         std::string ip_address,
                                         port_t port_no) {
  boost::mutex::scoped_lock lock(address_book_mutex_);
  address_book_[worker_id] = std::make_pair(ip_address, port_no);
  return true;
}


bool WorkerDataExchanger::SendSerializedData(const job_id_t& receive_job_id,
                                             const job_id_t& mega_rcr_job_id,
                                             const worker_id_t& worker_id,
                                             const SerializedData& ser_data,
                                             const data_version_t& version,
                                             const template_id_t& template_generation_id) {
  char header[MAX_HEADER_SIZE];
  int header_size = snprintf(header, MAX_HEADER_SIZE, "%lu %lu %lu %lu %lu;",
                             receive_job_id,
                             mega_rcr_job_id,
                             ser_data.size(),
                             version,
                             template_generation_id);
  if (header_size >= MAX_HEADER_SIZE || header_size < 0) {
    dbg(DBG_ERROR, "ERROR: failed building the header %d.\n", header_size); // NOLINT
    assert(false);
  }
  ser_data.set_header(header, (size_t)(header_size));

  WorkerDataExchangerConnection* connection = NULL;
  {
    boost::mutex::scoped_lock lock1(send_connection_mutex_);
    boost::mutex::scoped_lock lock2(address_book_mutex_);
    WorkerDataExchangerConnectionMap::iterator iter =
      send_connections_.find(worker_id);
    if (iter == send_connections_.end()) {
      if (address_book_.count(worker_id) == 0) {
        dbg(DBG_ERROR, "ERROR: could not find the address info of the worker id: %lu.\n", worker_id); //NOLINT
        assert(false);
      } else {
        std::pair<std::string, port_t> add = address_book_[worker_id];
        CreateNewSendConnection(worker_id, add.first, add.second);
      }
      connection = send_connections_[worker_id];
    } else {
      connection = iter->second;
    }
  }


#ifndef _NIMBUS_NO_NETWORK_LOG
  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff), "S %10.9lf j: %5.0lu s: %5.0lu",
      Log::GetRawTime(), job_id, size);
  log_.log_WriteToFile(std::string(buff));
#endif

  // For asynchronous write
  {
    boost::mutex::scoped_lock lock(*(connection->mutex()));
    std::list<SerializedData>* q = connection->send_queue();
    if (q->size() == 0) {
      std::vector<boost::asio::const_buffer> buffers;
      // Note: header() method returns a reference so the hader buffer remains
      // valid as long as ser_data is valid. - omidm
      buffers.push_back(boost::asio::buffer(ser_data.header(), ser_data.header().size()));
      buffers.push_back(boost::asio::buffer(ser_data.data_ptr().get(), ser_data.size()));

      boost::asio::async_write(*(connection->socket()),
                               buffers,
                               boost::bind(&WorkerDataExchanger::HandleWrite,
                                           this,
                                           connection,
                                           boost::asio::placeholders::error,
                                           boost::asio::placeholders::bytes_transferred));
    }
    q->push_back(ser_data);
  }


  // For synchronous write
  // {
  //   boost::mutex::scoped_lock lock(*(connection->mutex()));
  //   boost::system::error_code ignored_error;
  //   boost::asio::write(*(connection->socket()),
  //       boost::asio::buffer(header, header_size),
  //       boost::asio::transfer_all(), ignored_error);
  //   boost::asio::write(*(connection->socket()),
  //       boost::asio::buffer(ser_data.data_ptr().get(), ser_data.size()),
  //     boost::asio::transfer_all(), ignored_error);
  // }

  return true;
}

void WorkerDataExchanger::HandleWrite(WorkerDataExchangerConnection* connection,
                                      const boost::system::error_code& error,
                                      size_t bytes_transferred) {
  if (error) {
    dbg(DBG_NET|DBG_ERROR, "Error %s.\n", error.message().c_str());
    return;
  }

  {
    boost::mutex::scoped_lock lock(*(connection->mutex()));
    std::list<SerializedData>* q = connection->send_queue();
    assert(q->size() > 0);
    q->pop_front();
    if (q->size() > 0) {
      SerializedData ser_data = q->front();

      std::vector<boost::asio::const_buffer> buffers;
      // Note: header() method returns a reference so the hader buffer remains
      // valid as long as ser_data is valid. - omidm
      buffers.push_back(boost::asio::buffer(ser_data.header(), ser_data.header().size()));
      buffers.push_back(boost::asio::buffer(ser_data.data_ptr().get(), ser_data.size()));

      boost::asio::async_write(*(connection->socket()),
                               buffers,
                               boost::bind(&WorkerDataExchanger::HandleWrite,
                                           this,
                                           connection,
                                           boost::asio::placeholders::error,
                                           boost::asio::placeholders::bytes_transferred));
    }
  }
}

bool WorkerDataExchanger::CreateNewSendConnection(worker_id_t worker_id,
                                                  std::string ip_address,
                                                  port_t port_no) {
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
  // Set the tcp send and receive buf size.
  // Note: you may have to increase the OS limits first.
  // Look at the nimbus/scripts/configure_tcp.sh for help.
#ifndef __MACH__
  boost::asio::socket_base::send_buffer_size s_option(DATA_EXCHANGER_TCP_SEND_BUF_SIZE);
  boost::asio::socket_base::receive_buffer_size r_option(DATA_EXCHANGER_TCP_RECEIVE_BUF_SIZE);
  connection->socket()->set_option(s_option);
  connection->socket()->set_option(r_option);
#endif
  // Turn of Nagle algorithm.
  boost::asio::ip::tcp::no_delay nd_option(TCP_NODELAY_OPTION);
  connection->socket()->set_option(nd_option);
  send_connections_[worker_id] = connection;
  return true;
}

WorkerDataExchangerConnectionMap* WorkerDataExchanger::send_connections() {
  return &send_connections_;
}

WorkerDataExchangerConnectionList* WorkerDataExchanger::receive_connections() {
  return &receive_connections_;
}

void WorkerDataExchanger::WriteTimeDriftToLog(double drift) {
  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff), "D %10.9lf", drift);
  log_.log_WriteToFile(std::string(buff));
}

}  // namespace nimbus
