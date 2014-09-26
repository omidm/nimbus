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
  * Worker monitor keeps track of the workers activity. It captures their idle,
  * blocked, and active time and detects the stragglers and fast workers for
  * load balancing strategy.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_WORKER_MONITOR_H_
#define NIMBUS_SCHEDULER_WORKER_MONITOR_H_

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <map>
#include <list>
#include <vector>
#include <string>
#include "shared/nimbus_types.h"
#include "scheduler/worker_profile.h"
#include "scheduler/region_map.h"

namespace nimbus {

  class WorkerMonitor {
  public:
    typedef boost::unordered_map<worker_id_t, WorkerProfile*> Map;
    typedef Map::iterator MapIter;

    WorkerMonitor();
    virtual ~WorkerMonitor();

    bool AddWorker(worker_id_t worker_id);

    bool RemoveWorker(worker_id_t worker_id);

    bool AddReadyJob(worker_id_t worker_id, job_id_t job_id);

    bool AddBlockedJob(worker_id_t worker_id, job_id_t job_id);

    bool NotifyJobDone(worker_id_t worker_id, job_id_t job_id);

    void ResetWorkerTimers();

    size_t GetStragglers(std::list<worker_id_t> *list, size_t percentile);

    size_t GetFastWorkers(std::list<worker_id_t> *list, size_t percentile);

    std::string PrintStats();

  private:
    WorkerMonitor(const WorkerMonitor& other) {}

    Map map_;
    size_t worker_num_;
  };

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_WORKER_MONITOR_H_
