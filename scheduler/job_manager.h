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
  * Scheduler Job Manager object. This module serves the scheduler by providing
  * facilities about jobs ready to be maped, and their dependencies.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_JOB_MANAGER_H_
#define NIMBUS_SCHEDULER_JOB_MANAGER_H_

#include <boost/thread.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <map>
#include <set>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "shared/graph.h"
#include "scheduler/job_entry.h"
#include "scheduler/version_manager.h"
#include "scheduler/physical_data.h"
#include "scheduler/ldl_map.h"

namespace nimbus {
class JobManager {
  public:
    explicit JobManager();
    virtual ~JobManager();

    bool AddComputeJobEntry(
        const std::string& job_name,
        const job_id_t& job_id,
        const IDSet<logical_data_id_t>& read_set,
        const IDSet<logical_data_id_t>& write_set,
        const IDSet<job_id_t>& before_set,
        const IDSet<job_id_t>& after_set,
        const job_id_t& parent_job_id,
        const job_id_t& future_job_id,
        const bool& sterile,
        const Parameter& params);

    bool AddExplicitCopyJobEntry();

    bool AddKernelJobEntry();

    bool AddMainJobEntry(const job_id_t& job_id);

    bool AddCreateDataJobEntry(const job_id_t& job_id);

    bool AddLocalCopyJobEntry(const job_id_t& job_id);

    bool AddRemoteCopySendJobEntry(const job_id_t& job_id);

    bool AddRemoteCopyReceiveJobEntry(const job_id_t& job_id);

    bool AddFutureJobEntry(const job_id_t& job_id);

    bool AddJobEntryIncomingEdges(JobEntry *job);

    void ReceiveMetaBeforeSetDepthVersioningDependency(JobEntry* job);

    void PassMetaBeforeSetDepthVersioningDependency(JobEntry* job);

    bool GetJobEntry(job_id_t job_id, JobEntry*& job);

    bool RemoveJobEntry(JobEntry* job);

    bool RemoveJobEntry(job_id_t job_id);

    size_t GetJobsReadyToAssign(JobEntryList* list, size_t max_num);

    size_t RemoveObsoleteJobEntries();

    void CleanLdlMap();

    void NotifyJobDone(job_id_t job_id);

    void NotifyJobAssignment(JobEntry *job, const SchedulerWorker* worker);

    void DefineData(job_id_t job_id, logical_data_id_t ldid);

    size_t GetJobsNeedDataVersion(JobEntryList* list,
        VersionedLogicalData vld);

    bool AllJobsAreDone();

    void UpdateJobBeforeSet(JobEntry* job);

    void UpdateBeforeSet(IDSet<job_id_t>* before_set);

    bool CausingUnwantedSerialization(JobEntry* job,
        const logical_data_id_t& l_id, const PhysicalData& pd);

    void set_ldo_map_p(const std::map<logical_data_id_t, LogicalDataObject*>* ldo_map_p);

    Graph<JobEntry, job_id_t> *job_graph_p();

  private:
    Graph<JobEntry, job_id_t> job_graph_;
    VersionManager version_manager_;
    LdlMap ldl_map_;
    const std::map<logical_data_id_t, LogicalDataObject*>* ldo_map_p_;

    JobEntryMap jobs_done_;

    JobEntryMap jobs_need_version_;
    std::map<job_id_t, JobEntryList> pass_version_;

    JobEntryMap jobs_ready_to_assign_;
    JobEntryMap jobs_pending_to_assign_;

    Log log_version_;
    Log log_merge_;
    Log log_lookup_;
    Log log_sterile_;
    Log log_nonsterile_;
    size_t lookup_count_;


    size_t ResolveDataVersions();

    void PassDataVersionToJob(JobEntry *job,
                              const JobEntryList& from_jobs);

    bool LookUpVersion(JobEntry *job,
                       logical_data_id_t ldid,
                       data_version_t *version);

    bool JobVersionIsComplete(JobEntry *job);




    boost::mutex pass_version_mutex_;
    boost::condition_variable pass_version_put_cond_;
    boost::condition_variable pass_version_draw_cond_;
    size_t pass_version_in_progress_;
    void WaitToPassAllVersions();
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_JOB_MANAGER_H_
