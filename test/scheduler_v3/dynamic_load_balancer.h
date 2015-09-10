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
  * This is dynamic load balancer that has load ba;ancer and base class and
  * overrides some of its methods.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_TEST_SCHEDULER_V3_DYNAMIC_LOAD_BALANCER_H_
#define NIMBUS_TEST_SCHEDULER_V3_DYNAMIC_LOAD_BALANCER_H_

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <utility>
#include "shared/nimbus_types.h"
#include "scheduler/job_entry.h"
#include "scheduler/job_profile.h"
#include "scheduler/data_manager.h"
#include "scheduler/job_manager.h"
#include "scheduler/region_map.h"
#include "scheduler/straggler_map.h"
#include "scheduler/job_assigner.h"
#include "scheduler/load_balancer.h"
#include "scheduler/worker_monitor.h"
#include "shared/cluster.h"
#include "shared/id_maker.h"
#include "shared/scheduler_server.h"
#include "shared/geometric_region.h"
#include "shared/graph.h"

namespace nimbus {

  class DynamicLoadBalancer : public LoadBalancer {
  public:
    DynamicLoadBalancer();
    virtual ~DynamicLoadBalancer();

    virtual void Run();

    virtual void NotifyRegisteredWorker(SchedulerWorker *worker);

    virtual bool NotifyDownWorker(worker_id_t worker_id);

    virtual bool BalanceLoad(counter_t query_id);

    bool SetWorkerToAssignJob(JobEntry *job);

  private:
    typedef std::map<worker_id_t, SchedulerWorker*> WorkerMap;
    typedef std::pair<worker_id_t, worker_id_t>  Exchange;
    typedef std::list<Exchange>  BlackList;

    DynamicLoadBalancer(const DynamicLoadBalancer& other) {}

    Log log_;
    size_t worker_num_;
    bool init_phase_;
    RegionMap region_map_;
    WorkerMap worker_map_;
    boost::recursive_mutex mutex_;

    int64_t last_global_run_time_;
    Exchange last_exchange_;
    BlackList black_list_;

    void SplitDimensions(size_t worker_num);

    void BuildRegionMap();

    void UpdateRegionMap();

    int64_t GetGlobalRunTime(const StatQuery *st);

    bool InBlackList(worker_id_t w2, worker_id_t w1);

    bool CheckPossibleFluctuation(worker_id_t w2,
                                  worker_id_t w1,
                                  const StatQuery *st);
  };

}  // namespace nimbus

#endif  // NIMBUS_TEST_SCHEDULER_V3_DYNAMIC_LOAD_BALANCER_H_
