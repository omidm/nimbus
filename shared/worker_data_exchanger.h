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


#ifndef NIMBUS_SHARED_WORKER_DATA_EXCHANGER_H_
#define NIMBUS_SHARED_WORKER_DATA_EXCHANGER_H_

#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>  // NOLINT
#include <map>
#include <list>
#include <vector>
#include <sstream>
#include <string>
#include <utility>
#include "shared/dbg.h"
#include "shared/log.h"
#include "shared/parser.h"
#include "shared/nimbus_types.h"
#include "shared/serialized_data.h"
#include "shared/worker_data_exchanger_connection.h"

namespace nimbus {

using boost::asio::ip::tcp;

class WorkerDataExchanger {
 public:
  // All fields in Event should be primary types or you need to add copy
  // constructor.
  class Event {
    public:
      Event(const job_id_t& receive_job_id,
            const job_id_t& mega_rcr_job_id,
            const data_version_t& version,
            SerializedData *ser_data,
            const template_id_t& template_generation_id)
        : receive_job_id_(receive_job_id),
          mega_rcr_job_id_(mega_rcr_job_id),
          version_(version),
          ser_data_(ser_data),
          template_generation_id_(template_generation_id) {}
      ~Event() {}

      job_id_t receive_job_id_;
      job_id_t mega_rcr_job_id_;
      data_version_t version_;
      SerializedData *ser_data_;
      template_id_t template_generation_id_;
  };

  typedef std::list<Event> EventList;

  explicit WorkerDataExchanger(port_t port_no);
  virtual ~WorkerDataExchanger();

  virtual void Run();

  virtual bool AddContactInfo(worker_id_t worker_id,
                              std::string ip_address, port_t port_no);

  virtual size_t PullReceiveEvents(EventList *events,
                                      size_t max_num);

  virtual bool SendSerializedData(job_id_t receive_job_id,
                                  job_id_t mega_rcr_job_id,
                                  worker_id_t worker_id,
                                  SerializedData& ser_data,
                                  data_version_t version,
                                  template_id_t template_generation_id);

  WorkerDataExchangerConnectionMap* send_connections();

  WorkerDataExchangerConnectionList* receive_connections();

  void WriteTimeDriftToLog(double drift);

 private:
  typedef std::map<worker_id_t, std::pair<std::string, port_t> >AddressBook;

  Log log_;
  port_t listening_port_;
  tcp::acceptor* acceptor_;
  boost::asio::io_service* io_service_;

  AddressBook address_book_;
  boost::mutex address_book_mutex_;

  EventList event_list_;
  boost::mutex event_list_mutex_;

  WorkerDataExchangerConnectionMap send_connections_;
  boost::mutex send_connection_mutex_;

  WorkerDataExchangerConnectionList receive_connections_;
  boost::mutex receive_connection_mutex_;

  virtual void ListenForNewConnections();

  virtual void HandleAccept(WorkerDataExchangerConnection* connection,
                            const boost::system::error_code& error);

  virtual void HandleRead(WorkerDataExchangerConnection* connection,
                          const boost::system::error_code& error,
                          size_t bytes_transferred);

  virtual void HandleWrite(WorkerDataExchangerConnection* connection,
                           const boost::system::error_code& error,
                           size_t bytes_transferred);

  virtual size_t ReadData(WorkerDataExchangerConnection* connection,
                          char* buffer,
                          size_t size);

  virtual size_t ReadHeader(WorkerDataExchangerConnection* connection,
                            char* buffer,
                            size_t size);

  virtual bool CreateNewSendConnection(worker_id_t worker_id,
                                       std::string ip_address,
                                       port_t port_no);

  virtual void AddReceiveConnection(WorkerDataExchangerConnection* connection);

  virtual void AddSerializedData(job_id_t receive_job_id,
                                 job_id_t mega_rcr_job_id,
                                 SerializedData* ser_data,
                                 data_version_t version,
                                 template_id_t template_generation_id);
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_WORKER_DATA_EXCHANGER_H_
