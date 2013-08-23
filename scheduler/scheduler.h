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
  * Nimbus scheduler. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_SCHEDULER_H_
#define NIMBUS_SCHEDULER_SCHEDULER_H_

#include <boost/thread.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "shared/nimbus.h"
#include "shared/scheduler_server.h"
#include "shared/cluster.h"
#include "shared/parser.h"

namespace nimbus {
class Scheduler {
  public:
    explicit Scheduler(unsigned int listening_port);
    virtual ~Scheduler() {}

    virtual void run();
    virtual void schedulerCoreProcessor() {}
    virtual void loadClusterMap(std::string) {}
    virtual void deleteWorker(Worker * worker) {}
    virtual void addWorker(Worker * worker) {}
    virtual Worker* getWorker(int workerId) {return worker_map_[workerId];}

  protected:
    SchedulerServer* server_;
    CmSet user_command_set_;
    CommandSet worker_command_set_;

  private:
    virtual void setupUserInterface();
    virtual void setupWorkerInterface();
    virtual void getUserCommand();
    virtual void loadUserCommands();
    virtual void loadWorkerCommands();

    boost::thread* user_interface_thread_;
    boost::thread* worker_thread_;
    Computer host_;
    uint16_t port_;
    uint64_t appId_;
    // AppMap app_map_;
    WorkerMap worker_map_;
    ClusterMap cluster_map_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_SCHEDULER_H_
