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
  * Abstraction of a connection between two Nimbus workers.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#ifndef NIMBUS_SHARED_WORKER_DATA_EXCHANGER_CONNECTION_H_
#define NIMBUS_SHARED_WORKER_DATA_EXCHANGER_CONNECTION_H_

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <map>
#include <list>
#include <string>
#include "shared/dbg.h"
#include "shared/parser.h"
#include "shared/nimbus_types.h"

namespace nimbus {

using boost::asio::ip::tcp;

class WorkerDataExchangerConnection {
 public:
  explicit WorkerDataExchangerConnection(tcp::socket* sock);
  virtual ~WorkerDataExchangerConnection();

  void AllocateData(size_t size);

  void AppendData(char* buffer, size_t size);

  tcp::socket* socket();
  job_id_t job_id();
  char* data_ptr();
  char* read_buffer();
  size_t existing_bytes();
  size_t read_buffer_max_length();
  size_t data_length();
  size_t remaining_data_length();
  bool middle_of_data();
  bool middle_of_header();

  void set_job_id(job_id_t job_id);
  void set_data_length(size_t len);
  void set_existing_bytes(size_t len);
  void set_middle_of_data(bool flag);
  void set_middle_of_header(bool flag);

 private:
  tcp::socket* socket_;
  job_id_t job_id_;
  char* data_ptr_;
  char* data_cursor_;
  char* read_buffer_;
  size_t existing_bytes_;
  size_t data_length_;
  size_t remaining_data_length_;
  bool middle_of_data_;
  bool middle_of_header_;
};

typedef std::list<WorkerDataExchangerConnection*>
WorkerDataExchangerConnectionList;

typedef std::map<worker_id_t, WorkerDataExchangerConnection*>
WorkerDataExchangerConnectionMap;

}  // namespace nimbus

#endif  // NIMBUS_SHARED_WORKER_DATA_EXCHANGER_CONNECTION_H_
