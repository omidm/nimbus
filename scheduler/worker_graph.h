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
  * Scheduler data structure for monitoring computational resources and
  * their capabilities.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */


#ifndef NIMBUS_SCHEDULER_WORKER_GRAPH_H_
#define NIMBUS_SCHEDULER_WORKER_GRAPH_H_

#include <map>
#include <string>
#include <vector>
#include "shared/nimbus_types.h"
#include "shared/cluster.h"
#include "shared/protobuf_compiled/workergraphmessage.pb.h"
#include "scheduler/scheduler_worker.h"

namespace nimbus {

  typedef std::vector<worker_id_t> WorkerIdVector;
  typedef std::vector<SchedulerWorker*> SchedulerWorkerVector;
  typedef std::map<worker_id_t, SchedulerWorker*> SchedulerWorkerTable;
  typedef std::map<worker_id_t, Computer*> ComputerTable;

  class WorkerGraph {
  public:
    WorkerGraph();
    virtual ~WorkerGraph();

    // The WorkerGraph does not take responsibility
    // for managing SchedulerWorker objects.
    bool AddWorker(SchedulerWorker* worker);
    bool RemoveWorker(SchedulerWorker* worker);
    int AllWorkers(WorkerIdVector* dest);

    // Assumes the WorkerGraphMessage* will be freed
    // elsewhere

    bool ProcessMessage(WorkerGraphMessage* message);
    const SchedulerWorker* WorkerById(worker_id_t id);
    const Computer* ComputerById(worker_id_t id);
    uint32_t InterComputerMbps(worker_id_t src, worker_id_t dest);
    uint32_t InterComputerMicroSec(worker_id_t src, worker_id_t dest);

  private:
    ClusterMap cluster_map_;
    SchedulerWorkerTable worker_table_;
    ComputerTable computer_table_;

    bool ProcessWorkerRegisterMessage(const WorkerRegisterMessage& msg);
    bool ProcessSwitchRegisterMessage(const SwitchRegisterMessage& msg);
    bool ProcessWorkerLinkMessage(const WorkerLinkMessage& msg);
    bool ProcessSwitchLinkMessage(const SwitchLinkMessage& msg);
    bool ProcessUpdateMessage(const UpdateMessage& msg);
  };

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_WORKER_GRAPH_H_
