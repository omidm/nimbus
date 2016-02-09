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
  * StragglerMap is the data structure that keeps track of the load imbalance
  * among the workers. LoadBalancer uses this data structure to record the
  * information about the blamed workers that blocked a faster worker, and then
  * can query the data structure to fine the most imbalance workers
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_STRAGGLER_MAP_H_
#define NIMBUS_SRC_SCHEDULER_STRAGGLER_MAP_H_

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <map>
#include <list>
#include <vector>
#include "src/shared/nimbus_types.h"
#include "src/scheduler/job_entry.h"
#include "src/scheduler/job_profile.h"
#include "src/scheduler/data_manager.h"
#include "src/scheduler/job_manager.h"
#include "src/scheduler/region_map.h"
#include "src/shared/cluster.h"
#include "src/shared/geometric_region.h"
#include "src/shared/graph.h"

namespace nimbus {

  class StragglerMap {
  public:
    typedef boost::unordered_map<worker_id_t, size_t> Table;
    typedef Table::iterator TableIter;
    typedef boost::unordered_map<worker_id_t, Table*> Map;
    typedef Map::iterator MapIter;

    StragglerMap();
    virtual ~StragglerMap();

    void AddRecord(const worker_id_t& suffered,
                   const worker_id_t& blamed);

    void ClearRecords();

    bool GetMostImbalanceWorkers(worker_id_t *fast,
                                 worker_id_t *slow);

    size_t LookUp(const worker_id_t& suffered,
                  const worker_id_t& blamed);

  private:
    StragglerMap(const StragglerMap& other) {}

    Map map_;
  };

}  // namespace nimbus

#endif  // NIMBUS_SRC_SCHEDULER_STRAGGLER_MAP_H_
