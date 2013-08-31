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
  * Simple Nimbus Worker. It runs the commands it receives from the scheduler
  * without special discretion. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./simple_worker.h"

SimpleWorker::SimpleWorker(std::string scheduler_ip, port_t scheduler_port,
        port_t listening_port, Application * a)
: Worker(scheduler_ip, scheduler_port, listening_port, a) {
}

void SimpleWorker::WorkerCoreProcessor() {
  std::cout << "Simple Worker Core Processor" << std::endl;

  // Update the addressbook manually instead of query from scheduler;
  data_exchanger_->AddContactInfo(1, WORKER_IP_1, WORKER_PORT_1);
  data_exchanger_->AddContactInfo(2, WORKER_IP_2, WORKER_PORT_2);
  data_exchanger_->AddContactInfo(3, WORKER_IP_3, WORKER_PORT_3);
  data_exchanger_->AddContactInfo(4, WORKER_IP_4, WORKER_PORT_4);

  // So we have time to launch both workers, will remove when added
  // async_connect, and query from scheduler;
  // sleep(5);

  // Send some data;
  if (listening_port_ == WORKER_PORT_1) {
    std::string data_s = "hello hello hello to WORKER 2";
    SerializedData ser_data_s((char*)(data_s.c_str()), data_s.length()); // NOLINT
    data_exchanger_->SendSerializedData(201, 2, ser_data_s);
  }

  if (listening_port_ == WORKER_PORT_2) {
    std::string data_s = "Long Long Long Long Long Long Long Long to WORKER 1";
    SerializedData ser_data_s((char*)(data_s.c_str()), data_s.length()); // NOLINT
    data_exchanger_->SendSerializedData(101, 1, ser_data_s);
  }

  if (listening_port_ == WORKER_PORT_2) {
    std::string data_s = "Second Message to WORKER 1";
    SerializedData ser_data_s((char*)(data_s.c_str()), data_s.length()); // NOLINT
    data_exchanger_->SendSerializedData(102, 1, ser_data_s);
  }

  // Receive the data.
  if (listening_port_ == WORKER_PORT_1) {
    SerializedData ser_data_r;
    while (!data_exchanger_->ReceiveSerializedData(102, ser_data_r)) {
      // std::cout << "Waiting to receive data " << 1 << " ..." << std::endl;
      // sleep(1);
    }
    char * buf_r = new char[ser_data_r.size() + 1];
    memcpy(buf_r, ser_data_r.data_ptr(), ser_data_r.size());
    buf_r[ser_data_r.size()] = '\0';
    std::string data_r(buf_r);
    std::cout << "Received data. :) " << 102 << " : " << data_r << std::endl;
  }

  if (listening_port_ == WORKER_PORT_1) {
    SerializedData ser_data_r;
    while (!data_exchanger_->ReceiveSerializedData(101, ser_data_r)) {
      // std::cout << "Waiting to receive data " << 1 << " ..." << std::endl;
      // sleep(1);
    }
    char * buf_r = new char[ser_data_r.size() + 1];
    memcpy(buf_r, ser_data_r.data_ptr(), ser_data_r.size());
    buf_r[ser_data_r.size()] = '\0';
    std::string data_r(buf_r);
    std::cout << "Received data. :) " << 101 << " : " << data_r << std::endl;
  }

  if (listening_port_ == WORKER_PORT_2) {
    SerializedData ser_data_r;
    while (!data_exchanger_->ReceiveSerializedData(201, ser_data_r)) {
      // std::cout << "Waiting to receive data " << 1 << " ..." << std::endl;
      // sleep(1);
    }
    char * buf_r = new char[ser_data_r.size() + 1];
    memcpy(buf_r, ser_data_r.data_ptr(), ser_data_r.size());
    buf_r[ser_data_r.size()] = '\0';
    std::string data_r(buf_r);
    std::cout << "Received data. :) " << 201 << " : " << data_r << std::endl;
  }

  std::cout << "Just before the main loop of core processor." << std::endl;

  while (true) {
    // sleep(1);
    // std::cout << "While Loop Worker ID: " << id_ << std::endl;
    // Log::dbg_printLine("Worker running core loop.", INFO);
    // std::string str = "createjob name:main id:{0} read:{1,2} write:{1,2} ";
    // str += " before:{} after:{1,2,3} type:operation param:t=20,g=6";
    // SchedulerCommand cm(str);
    // std::cout << "Sending command: " << cm.toString() << std::endl;
    // client->sendCommand(&cm);
    SchedulerCommand* comm = client_->receiveCommand();
    if (comm != NULL) {
      std::cout << "Received command: " << comm->toString() << std::endl;
      ProcessSchedulerCommand(comm);
      delete comm;
    }
  }
}
