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
#include <algorithm>
#include "shared/scheduler_client.h"

using boost::asio::ip::tcp;
using namespace nimbus; // NOLINT

#define CLIENT_BUFSIZE 4096000

SchedulerClient::SchedulerClient(std::string scheduler_ip,
                                 port_t scheduler_port)
  : scheduler_ip_(scheduler_ip),
    scheduler_port_(scheduler_port) {
  read_buffer_ = new boost::asio::streambuf();
  io_service_ = new boost::asio::io_service();
  socket_ = new tcp::socket(*io_service_);
  command_num_ = 0;
  pthread_mutex_init(&send_lock_, NULL);
  byte_array_ = new char[CLIENT_BUFSIZE];
  existing_bytes_ = 0;
  existing_offset_ = 0;
}

SchedulerClient::~SchedulerClient() {
  pthread_mutex_destroy(&send_lock_);
  delete read_buffer_;
}

SchedulerCommand* SchedulerClient::ReceiveCommand() {
  SchedulerCommand* com = NULL;
  pthread_mutex_lock(&send_lock_);

  if (existing_bytes_ >= sizeof(SchedulerCommand::length_field_t)) {
    // We have at least a header field -- see if we have a whole message
    int header_len = sizeof(SchedulerCommand::length_field_t);
    char* read_ptr = &byte_array_[existing_offset_];

    SchedulerCommand::length_field_t len;
    char* len_ptr = reinterpret_cast<char*>(&len);
    memcpy(len_ptr, read_ptr, header_len);
    len = (uint32_t)ntohl(len);

    // We have a whole message
    if (existing_bytes_ >= len) {
      std::string input(read_ptr + header_len, len - header_len);

      // We've read out a command worth of bytes; in the case that
      // there are no bytes remaining, restart at beginning of buffer
      existing_bytes_ -= len;
      if (existing_bytes_ == 0) {
        existing_offset_ = 0;
      } else {
        existing_offset_ += len;
      }
      if (SchedulerCommand::GenerateSchedulerCommandChild(input,
                                                          scheduler_command_table_,
                                                          com)) {
        dbg(DBG_NET, "Scheduler client received command %s\n",
            com->ToString().c_str());
      } else {
        com = NULL;
        dbg(DBG_NET, "Ignored unknown command: %s.\n", input.c_str());
      }
      pthread_mutex_unlock(&send_lock_);
      return com;
    } else {
      // We don't have a whole message: move the fragment to the beginning,
      // so that we don't run off the end of the buffer after lots of reads
      // with incomplete commands.
      memmove(byte_array_, byte_array_ + existing_offset_, existing_bytes_);
      existing_offset_ = 0;
    }
  }
  // If we reach here there wasn't enough data in the buffer for a message.
  // Try reading, then try calling ourselves again. The level of recursion
  // Is expected to be very small.

  boost::system::error_code ignored_error;
  uint32_t bytes_available = socket_->available(ignored_error);
  bytes_available = std::min(bytes_available, CLIENT_BUFSIZE - existing_offset_);
  boost::asio::streambuf::mutable_buffers_type bufs =
    read_buffer_->prepare(bytes_available);
  std::size_t bytes_read = socket_->receive(bufs);
  read_buffer_->commit(bytes_read);
  std::istream input(read_buffer_);
  assert(bytes_read == bytes_available);
  input.read(byte_array_ + existing_offset_ + existing_bytes_, bytes_available);
  existing_bytes_ += bytes_available;
  pthread_mutex_unlock(&send_lock_);

  if (bytes_available > 0) {
    return ReceiveCommand();
  } else {
    return NULL;
  }
}

void SchedulerClient::SendCommand(SchedulerCommand* command) {
  // static double serialization_time = 0;
  // static double buffer_time = 0;
  // struct timespec start_time;
  // clock_gettime(CLOCK_REALTIME, &start_time);

  std::string data = command->ToNetworkData();
  SchedulerCommand::length_field_t len;
  len = htonl((uint32_t)(data.length() + sizeof(len)));
  std::string msg;
  msg.append((const char*)&len, sizeof(len));
  msg.append(data.c_str(), data.length());

  dbg(DBG_NET, "Client sending command of length %i: %s\n", data.length(), command->ToString().c_str()); // NOLINT

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
