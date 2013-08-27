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
#include <string>
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "shared/parser.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class WorkerDataExchangerConnection {
 public:
  explicit WorkerDataExchangerConnection(tcp::socket* sock);
  virtual ~WorkerDataExchangerConnection();

  boost::asio::streambuf* read_buffer();
  tcp::socket* socket();
  bool data_fully_received();
  void data_fully_received(bool flag);

 private:
  boost::asio::streambuf* read_buffer_;
  tcp::socket* socket_;
  bool data_fully_received_;
};

typedef std::map<worker_id_t, WorkerDataExchangerConnection*>
WorkerDataExchangerConnectionMap;

}  // namespace nimbus

#endif  // NIMBUS_SHARED_WORKER_DATA_EXCHANGER_CONNECTION_H_
