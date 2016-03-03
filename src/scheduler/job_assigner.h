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
  * This is the base job assigner module. It provides methods required for
  * binding a job to a worker and submitting the compute and copy jobs to the
  * worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_JOB_ASSIGNER_H_
#define NIMBUS_SRC_SCHEDULER_JOB_ASSIGNER_H_

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <map>
#include <list>
#include <vector>
#include <string>
#include "src/shared/nimbus_types.h"
#include "src/scheduler/job_entry.h"
#include "src/scheduler/job_profile.h"
#include "src/scheduler/data_manager.h"
#include "src/scheduler/job_manager.h"
#include "src/scheduler/region_map.h"
#include "src/scheduler/straggler_map.h"
#include "src/scheduler/binding_template.h"
#include "src/shared/cluster.h"
#include "src/shared/id_maker.h"
#include "src/shared/scheduler_server.h"
#include "src/shared/geometric_region.h"
#include "src/shared/graph.h"

namespace nimbus {

  class Scheduler;
  class LoadBalancer;

  class JobAssigner {
  public:
    JobAssigner();
    virtual ~JobAssigner();

    void Run();

    virtual void set_id_maker(IDMaker *id_maker);
    virtual void set_server(SchedulerServer* server);
    virtual void set_job_manager(JobManager *job_manager);
    virtual void set_data_manager(DataManager *data_manager);
    virtual void set_load_balancer(LoadBalancer *load_balancer);
    virtual void set_thread_num(size_t thread_num);
    virtual void set_fault_tolerance_active(bool flag);
    virtual void set_scheduler(Scheduler *scheduler);
    virtual void set_ldo_map_p(const LdoMap* ldo_map_p);
    virtual void set_log_overhead(Log *log);

    virtual void AssignJobs(const JobEntryList& list);
    virtual void Reinitialize(checkpoint_id_t checkpoint_id);

  protected:
    IDMaker *id_maker_;
    SchedulerServer* server_;
    JobManager *job_manager_;
    DataManager *data_manager_;
    LoadBalancer *load_balancer_;
    size_t thread_num_;
    checkpoint_id_t checkpoint_id_;
    template_id_t template_generation_id_;
    bool fault_tolerance_active_;
    Scheduler *scheduler_;
    const LdoMap *ldo_map_p_;
    Log *log_overhead_;

    JobEntryList job_queue_;
    boost::recursive_mutex job_queue_mutex_;
    boost::condition_variable_any job_queue_cond_;
    size_t pending_assignment_;

    std::list<boost::thread*> job_assigner_threads_;

    virtual void Initialize();

    virtual void JobAssignerThread();

    virtual bool AssignJob(JobEntry* job);

    virtual bool AssignComplexJob(ComplexJobEntry* job);

    virtual size_t UpdateDataForCascading(BindingTemplate *bt,
                                          ComplexJobEntry *complex_job,
                                          const BindingTemplate::PatternEntry *pattern,
                                          data_version_t new_version_diff);

    virtual bool PrepareDataForJobAtWorker(JobEntry* job,
                                           SchedulerWorker* worker,
                                           logical_data_id_t l_id);

    virtual bool AllocateLdoInstanceToJob(JobEntry* job,
                                          LogicalDataObject* ldo,
                                          PhysicalData pd);

    virtual bool SaveJobContextForCheckpoint(JobEntry *job);

    virtual bool CreateDataAtWorker(SchedulerWorker* worker,
                                    LogicalDataObject* ldo,
                                    PhysicalData* created_data);

    virtual bool LocalCombineData(JobEntry* ref_job,
                                  const std::string& combiner,
                                  SchedulerWorker* worker,
                                  LogicalDataObject* ldo,
                                  PhysicalDataList* scratch_sources,
                                  PhysicalData* scratch_dest,
                                  BindingTemplate *bt);

    virtual job_id_t RemoteCopyData(SchedulerWorker* from_worker,
                                    SchedulerWorker* to_worker,
                                    LogicalDataObject* ldo,
                                    PhysicalData* from_data,
                                    PhysicalData* to_data,
                                    BindingTemplate *bt = NULL);

    virtual bool LocalCopyData(SchedulerWorker* worker,
                               LogicalDataObject* ldo,
                               PhysicalData* created_data,
                               PhysicalData* to_data,
                               BindingTemplate *bt = NULL);

    virtual bool SaveData(SchedulerWorker* worker,
                          LogicalDataObject* ldo,
                          PhysicalData* from_data,
                          checkpoint_id_t checkpoint_id);

    virtual bool LoadData(SchedulerWorker* worker,
                          LogicalDataObject* ldo,
                          data_version_t version,
                          PhysicalData* to_data,
                          std::string handle);

    virtual bool GetFreeDataAtWorker(SchedulerWorker* worker,
                                     LogicalDataObject* ldo,
                                     PhysicalData* free_data);

    virtual size_t GetObsoleteLdoInstancesAtWorker(SchedulerWorker* worker,
                                                   LogicalDataObject* ldo,
                                                   PhysicalDataList* dest);




    virtual bool SendComputeJobToWorker(SchedulerWorker* worker,
                                        JobEntry* job);

    virtual bool UpdateDataManagerByPatterns(
                    ComplexJobEntry* job,
                    const BindingTemplate *binding_template,
                    const std::vector<physical_data_id_t> *physical_ids);

    virtual bool QueryDataManagerForPatterns(
                    ComplexJobEntry* job,
                    const BindingTemplate *binding_template,
                    const std::vector<physical_data_id_t>*& physical_ids,
                    BindingTemplate::ExtraDependency *extra_dependency);

    class DataManagerQueryCache {
      public:
        DataManagerQueryCache();
        ~DataManagerQueryCache();

        enum State {
          INIT     = 0,
          VALID    = 1
        };

        void Invalidate();
        bool QueryIsHit(const std::string& record_name);
        bool Query(const std::string& record_name,
                   const std::vector<physical_data_id_t>*& phy_ids);
        void Learn(const std::string& record_name,
                   const std::vector<physical_data_id_t>* phy_ids);

      private:
        State state_;
        std::string record_name_;
        const std::vector<physical_data_id_t>* cached_phy_ids_;
        boost::mutex mutex_;

        bool SameQueries(const std::vector<physical_data_id_t>* phy_ids_1,
                         const std::vector<physical_data_id_t>* phy_ids_2);
    };

    DataManagerQueryCache dm_query_cache_;


    void FlushBatchUpdate();

    size_t update_batch_size_;
    std::vector<logical_data_id_t> batch_ldid_;
    std::vector<physical_data_id_t> batch_pdid_;
    std::vector<data_version_t> batch_base_version_;
    std::vector<data_version_t> batch_unit_diff_version_;
    std::vector<BindingTemplate::VersionType> batch_version_type_;
    job_id_t last_complex_job_id_;
    boost::mutex batch_update_mutex_;

  private:
    JobAssigner(const JobAssigner& other) {}
  };

}  // namespace nimbus

#endif  // NIMBUS_SRC_SCHEDULER_JOB_ASSIGNER_H_
