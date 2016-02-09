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

/**
  * Abstraction of a connection from a Nimbus scheduler to
  * a worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "src/shared/scheduler_server_connection.h"
#include "src/shared/dbg.h"

#define SERVER_CON_BUF_SIZE 4096000

using boost::asio::ip::tcp;

namespace nimbus {

SchedulerServerConnection::SchedulerServerConnection(tcp::socket* sock)
  :socket_(sock) {
  command_num_ = 0;
  existing_bytes_ = 0;
  mutex_ = new boost::mutex();
  read_buffer_ = new char[SERVER_CON_BUF_SIZE];
}

SchedulerServerConnection::~SchedulerServerConnection() {
  // FIXME: not actually cleaning up listening thread.
  delete mutex_;
  delete read_buffer_;
}

tcp::socket* SchedulerServerConnection::socket() {
  return socket_;
}

int SchedulerServerConnection::command_num() {
  return command_num_;
}

void SchedulerServerConnection::set_command_num(int n) {
  command_num_ = n;
}

char* SchedulerServerConnection::read_buffer() {
  return read_buffer_;
}

uint32_t SchedulerServerConnection::existing_bytes() {
  return existing_bytes_;
}

boost::mutex* SchedulerServerConnection::mutex() {
  return mutex_;
}

void SchedulerServerConnection::set_existing_bytes(uint32_t bytes) {
  existing_bytes_ = bytes;
}

uint32_t SchedulerServerConnection::max_read_buffer_length() {
  return SERVER_CON_BUF_SIZE;
}

}  // namespace nimbus
