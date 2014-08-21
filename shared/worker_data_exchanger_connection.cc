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
  * Abstraction of a connection between two Nimbus workers.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/worker_data_exchanger_connection.h"

#define WORKER_DATA_BUFSIZE 4096000
// #define WORKER_DATA_BUFSIZE 8

using boost::asio::ip::tcp;

namespace nimbus {

WorkerDataExchangerConnection::WorkerDataExchangerConnection(tcp::socket* sock)
  :socket_(sock) {
  read_buffer_ = new char[WORKER_DATA_BUFSIZE];
  existing_bytes_ = 0;
  middle_of_data_ = false;
  middle_of_header_ = true;
}

WorkerDataExchangerConnection::~WorkerDataExchangerConnection() {
  delete read_buffer_;
  delete socket_;
}

void WorkerDataExchangerConnection::AllocateData(size_t size) {
  data_ptr_ = boost::shared_array<char>(new char[size]);
  data_cursor_ = data_ptr_.get();
  remaining_data_length_ = size;
}

void WorkerDataExchangerConnection::AppendData(char* buffer, size_t size) {
  if (size > remaining_data_length_) {
    dbg(DBG_ERROR, "ERROR: Appending beyond the size of data, ignored.\n");
    exit(-1);
    return;
  }
  memcpy(data_cursor_, buffer, size);
  data_cursor_ += size;
  remaining_data_length_ -= size;
}

tcp::socket* WorkerDataExchangerConnection::socket() {
  return socket_;
}


boost::shared_array<char> WorkerDataExchangerConnection::data_ptr() {
  return data_ptr_;
}

char* WorkerDataExchangerConnection::read_buffer() {
  return read_buffer_;
}

size_t WorkerDataExchangerConnection::existing_bytes() {
  return existing_bytes_;
}

void WorkerDataExchangerConnection::set_existing_bytes(size_t len) {
  existing_bytes_ = len;
}

size_t WorkerDataExchangerConnection::data_length() {
  return data_length_;
}

data_version_t WorkerDataExchangerConnection::data_version() {
  return data_version_;
}

void WorkerDataExchangerConnection::set_data_length(size_t len) {
  data_length_ = len;
}

void WorkerDataExchangerConnection::set_data_version(data_version_t version) {
  data_version_ = version;
}

job_id_t WorkerDataExchangerConnection::job_id() {
  return job_id_;
}

void WorkerDataExchangerConnection::set_job_id(job_id_t job_id) {
  job_id_ = job_id;
}

bool WorkerDataExchangerConnection::middle_of_header() {
  return middle_of_header_;
}

void WorkerDataExchangerConnection::set_middle_of_header(bool flag) {
  middle_of_header_ = flag;
}

bool WorkerDataExchangerConnection::middle_of_data() {
  return middle_of_data_;
}

void WorkerDataExchangerConnection::set_middle_of_data(bool flag) {
  middle_of_data_ = flag;
}


size_t WorkerDataExchangerConnection::remaining_data_length() {
  return remaining_data_length_;
}

size_t WorkerDataExchangerConnection::read_buffer_max_length() {
  return WORKER_DATA_BUFSIZE;
}


}  // namespace nimbus

