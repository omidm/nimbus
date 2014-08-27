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
#include <list>
#include <vector>
#include <string>
#include "shared/nimbus_types.h"
#include "scheduler/job_entry.h"
#include "scheduler/job_profile.h"
#include "scheduler/data_manager.h"
#include "scheduler/job_manager.h"
#include "scheduler/region_map.h"
#include "scheduler/straggler_map.h"
#include "shared/cluster.h"
#include "shared/id_maker.h"
#include "shared/scheduler_server.h"
#include "shared/geometric_region.h"
#include "shared/graph.h"

namespace nimbus {

  class LoadBalancer {
  public:
    typedef std::map<worker_id_t, SchedulerWorker*> WorkerMap;
    typedef WorkerMap::iterator WorkerMapIter;
    typedef std::map<job_id_t, JobProfile*> JobHistory;

    LoadBalancer();
    explicit LoadBalancer(ClusterMap* cluster_map);
    virtual ~LoadBalancer();

    void Run();

    ClusterMap* cluster_map();

    void set_id_maker(IDMaker *id_maker);
    void set_server(SchedulerServer* server);
    void set_cluster_map(ClusterMap* cluster_map);
    void set_job_manager(JobManager *job_manager);
    void set_data_manager(DataManager *data_manager);

    void AssignJobs(const JobEntryList& list);

    bool GetWorkerToAssignJob(JobEntry *job, SchedulerWorker*& worker);

    void NotifyJobAssignment(const JobEntry *job, const SchedulerWorker* worker);

    void NotifyJobDone(const JobEntry *job);

    void NotifyRegisteredWorker(SchedulerWorker *worker);


  private:
    LoadBalancer(const LoadBalancer& other) {}

    Log log_;
    size_t worker_num_;
    GeometricRegion global_region_;
    int stamp_state_;
    boost::mutex stamp_mutex_;

    IDMaker *id_maker_;
    SchedulerServer* server_;
    ClusterMap* cluster_map_;
    JobManager *job_manager_;
    DataManager *data_manager_;

    JobHistory job_history_;
    boost::recursive_mutex job_history_mutex_;

    std::list<job_id_t> done_jobs_;

    RegionMap region_map_;
    boost::recursive_mutex region_map_mutex_;

    WorkerMap worker_map_;
    boost::recursive_mutex worker_map_mutex_;

    StragglerMap straggler_map_;
    boost::recursive_mutex straggler_map_mutex_;

    JobEntryList job_queue_;
    boost::recursive_mutex job_queue_mutex_;
    boost::condition_variable_any job_queue_cond_;
    size_t pending_assignment_;

    std::list<boost::thread*> job_assigner_threads_;

    bool update_;
    bool init_phase_;
    boost::recursive_mutex update_mutex_;
    boost::condition_variable_any update_cond_;

    std::map<worker_id_t, size_t> blame_map_;
    size_t blame_counter_;

    void Initialize();

    void InitializeRegionMap();

    void UpdateRegionMap();

    void JobAssignerThread();

    bool AssignJob(JobEntry* job);

     bool PrepareDataForJobAtWorker(JobEntry* job,
                                    SchedulerWorker* worker,
                                    logical_data_id_t l_id);

    bool AllocateLdoInstanceToJob(JobEntry* job,
                                  LogicalDataObject* ldo,
                                  PhysicalData pd);

    bool CreateDataAtWorker(SchedulerWorker* worker,
                            LogicalDataObject* ldo,
                            PhysicalData* created_data);

    bool RemoteCopyData(SchedulerWorker* from_worker,
                        SchedulerWorker* to_worker,
                        LogicalDataObject* ldo,
                        PhysicalData* from_data,
                        PhysicalData* to_data);

    bool LocalCopyData(SchedulerWorker* worker,
                       LogicalDataObject* ldo,
                       PhysicalData* created_data,
                       PhysicalData* to_data);

    bool GetFreeDataAtWorker(SchedulerWorker* worker,
                             LogicalDataObject* ldo,
                             PhysicalData* free_data);

    size_t GetObsoleteLdoInstancesAtWorker(SchedulerWorker* worker,
                                           LogicalDataObject* ldo,
                                           PhysicalDataVector* dest);




    bool SendComputeJobToWorker(SchedulerWorker* worker,
                                JobEntry* job);
  };

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_LOAD_BALANCER_H_
