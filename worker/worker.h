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
  * Nimbus abstraction of an application. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_WORKER_WORKER_H_
#define NIMBUS_WORKER_WORKER_H_

#include <boost/thread.hpp>
#include <string>
#include <vector>
#include <map>
#include "worker/data.h"
#include "worker/job.h"
#include "worker/application.h"
#include "shared/scheduler_client.h"
#include "shared/scheduler_command.h"
#include "shared/cluster.h"
#include "shared/parser.h"
#include "shared/log.h"

namespace nimbus {

class Worker;
typedef std::map<int, Worker*> WorkerMap;

class Worker {
 public:
  Worker(unsigned int port,  Application* application);

  virtual void run();
  virtual void workerCoreProcessor() {}
  virtual void processSchedulerCommand(SchedulerCommand* command);

 protected:
  SchedulerClient* client_;

 private:
  int id_;
  Log log;
  Computer host_;
  unsigned int port_;
  DataMap data_map_;
  JobMap job_map_;
  Application* application_;
  CmSet scheduler_command_set_;

  virtual void setupSchedulerInterface();

  virtual void addJob(Job* job);
  virtual void deleteJob(int id);
  virtual void deleteJob(Job* job) {}
  virtual void addData(Data* data);
  virtual void deleteData(int id);
  virtual void deleteData(Data* data) {}
  virtual void loadSchedulerCommands();
};

}  // namespace nimbus

#endif  // NIMBUS_WORKER_WORKER_H_
