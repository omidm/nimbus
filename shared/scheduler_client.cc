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

// Comment(quhang): I added the lock for sending commmand.
#include <pthread.h>
#include "shared/scheduler_client.h"

using boost::asio::ip::tcp;
using namespace nimbus; // NOLINT

SchedulerClient::SchedulerClient(std::string scheduler_ip,
                                 port_t scheduler_port)
  : scheduler_ip_(scheduler_ip),
    scheduler_port_(scheduler_port) {
  read_buffer_ = new boost::asio::streambuf();
  io_service_ = new boost::asio::io_service();
  socket_ = new tcp::socket(*io_service_);
  command_num_ = 0;
  pthread_mutex_init(&send_lock_, NULL);
}

SchedulerClient::~SchedulerClient() {
  pthread_mutex_destroy(&send_lock_);
}

SchedulerCommand* SchedulerClient::ReceiveCommand() {
  // boost::asio::read_until(*socket, *read_buffer, ';');
  // std::streamsize size = read_buffer->in_avail();
  SchedulerCommand* com = NULL;
  pthread_mutex_lock(&send_lock_);

  boost::system::error_code ignored_error;
  int bytes_available = socket_->available(ignored_error);

  boost::asio::streambuf::mutable_buffers_type bufs =
    read_buffer_->prepare(bytes_available);
  std::size_t bytes_read = socket_->receive(bufs);
  read_buffer_->commit(bytes_read);

  if (read_buffer_->size() >= sizeof(SchedulerCommand::length_field_t)) {
    SchedulerCommand::length_field_t len;
    char* ptr = reinterpret_cast<char*>(&len);
    std::istream input(read_buffer_);
    int header_len = sizeof(SchedulerCommand::length_field_t);
    input.read(ptr, header_len);
    len = ntohl(len);

    dbg(DBG_NET, "Reading a command of length %i.\n", len);
    // We have a complete command
    if (read_buffer_->size() >= (len - header_len)) {
      std::string command;
      command.resize(len - header_len);
      // input.seekg(header_len);
      // Then read the actual data
      input.read(&command[0], len - header_len);

      if (SchedulerCommand::GenerateSchedulerCommandChild(command,
                                                          scheduler_command_table_,
                                                          com)) {
        dbg(DBG_NET, "Scheduler client received command %s\n",
            com->toStringWTags().c_str());
      } else {
        com = NULL;
        dbg(DBG_NET, "Ignored unknown command: %s.\n", command.c_str());
      }
    }
  } else {
    // dbg(DBG_NET, "Read buffer has %i bytes.\n", read_buffer_->size());
  }
  pthread_mutex_unlock(&send_lock_);
  return com;
}

void SchedulerClient::SendCommand(SchedulerCommand* command) {
  // static double serialization_time = 0;
  // static double buffer_time = 0;
  // struct timespec start_time;
  // clock_gettime(CLOCK_REALTIME, &start_time);

  std::string data = command->toString();
  SchedulerCommand::length_field_t len;
  len = htonl((uint32_t)data.length() + sizeof(len));
  std::string msg;
  msg.append((const char*)&len, sizeof(len));
  msg.append(data.c_str(), data.length());

  dbg(DBG_NET, "Client sending command of length %i: %s\n", data.length(), command->toStringWTags().c_str()); // NOLINT

  // struct timespec t;
  // clock_gettime(CLOCK_REALTIME, &t);
  // serialization_time += difftime(t.tv_sec, start_time.tv_sec)
  //     + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
  // printf("Serialization command time %f\n", serialization_time);

  boost::system::error_code ignored_error;
  // clock_gettime(CLOCK_REALTIME, &start_time);

  pthread_mutex_lock(&send_lock_);
  boost::asio::write(*socket_, boost::asio::buffer(msg),
      boost::asio::transfer_all(), ignored_error);
  pthread_mutex_unlock(&send_lock_);
  // clock_gettime(CLOCK_REALTIME, &t);
  // buffer_time += difftime(t.tv_sec, start_time.tv_sec)
  //     + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
  // printf("Buffer command time %f\n", buffer_time);
}

void SchedulerClient::CreateNewConnections() {
  std::cout << "Opening connections." << std::endl;
  tcp::resolver resolver(*io_service_);
  tcp::resolver::query query(scheduler_ip_,
      boost::to_string(scheduler_port_));
  tcp::resolver::iterator iterator = resolver.resolve(query);
  boost::system::error_code error = boost::asio::error::host_not_found;
  socket_->connect(*iterator, error);
}

void SchedulerClient::Run() {
  CreateNewConnections();
  // io_service_->run();
}


void
SchedulerClient::set_scheduler_command_table(SchedulerCommand::PrototypeTable* cmt) {
  scheduler_command_table_ = cmt;
}
