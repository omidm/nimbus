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
  * This is the base class that serves the scheduler regarding job assignment
  * queries in the cluster. It tries to minimize the completion of simulation
  * by reducing the cost of communication by locality aware data placement and
  * mitigate the effect of stragglers in the system by adapting the job
  * assignment strategies to the dynamic changes of the cloud.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_LOAD_BALANCER_H_
#define NIMBUS_SCHEDULER_LOAD_BALANCER_H_

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <map>
#include <vector>
#include "shared/nimbus_types.h"
#include "scheduler/job_entry.h"
#include "scheduler/job_profile.h"
#include "shared/cluster.h"
#include "shared/geometric_region.h"
#include "shared/graph.h"

namespace nimbus {

  class LoadBalancer {
  public:
    typedef std::map<worker_id_t, GeometricRegion> RegionMap;
    typedef std::map<worker_id_t, SchedulerWorker*> WorkerMap;

    LoadBalancer();
    explicit LoadBalancer(ClusterMap* cluster_map);
    virtual ~LoadBalancer();

    void Run();

    ClusterMap* cluster_map();
    GeometricRegion global_region();

    void set_cluster_map(ClusterMap* cluster_map);
    void set_global_region(GeometricRegion global_region);

    bool GetWorkerToAssignJob(JobEntry *job, SchedulerWorker*& worker);

    void NotifyJobAssignment(const JobEntry *job, const SchedulerWorker* worker);

    void NotifyJobDone(const JobEntry *job);


  private:
    LoadBalancer(const LoadBalancer& other) {}

    ClusterMap* cluster_map_;
    GeometricRegion global_region_;
    Log log_;

    Graph<JobProfile, job_id_t> job_graph_;
    boost::mutex job_graph_mutex_;

    RegionMap region_map_;
    boost::mutex region_map_mutex_;

    WorkerMap worker_map_;
    boost::mutex worker_map_mutex_;


    bool updated_info_;
    boost::mutex updated_info_mutex_;
    boost::condition_variable update_info_cond_;

    void InitializeRegionMap();
    void UpdateRegionMap();
  };

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_LOAD_BALANCER_H_