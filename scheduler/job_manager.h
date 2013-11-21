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
#include <algorithm>
#include <map>
#include <set>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "scheduler/job_graph.h"
#include "scheduler/job_entry.h"

namespace nimbus {
class JobManager {
  public:
    explicit JobManager();
    virtual ~JobManager();

    bool AddJobEntry(const JobType& job_type,
        const std::string& job_name,
        const job_id_t& job_id,
        const IDSet<logical_data_id_t>& read_set,
        const IDSet<logical_data_id_t>& write_set,
        const IDSet<job_id_t>& before_set,
        const IDSet<job_id_t>& after_set,
        const job_id_t& parent_job_id,
        const Parameter& params);

    bool AddJobEntry(const JobType& job_type,
        const std::string& job_name,
        const job_id_t& job_id,
        const job_id_t& parent_job_id);


    bool GetJobEntry(job_id_t job_id, JobEntry*& job);

    bool RemoveJobEntry(JobEntry* job);

    bool RemoveJobEntry(job_id_t job_id);

    size_t GetJobsReadyToAssign(JobEntryList* list, size_t max_num);

    size_t RemoveObsoleteJobEntries();

    void JobDone(job_id_t job_id);

    void DefineData(job_id_t job_id, logical_data_id_t ldid);


  private:
    JobGraph job_graph_;

    bool ResolveJobDataVersions(JobEntry* job);
    size_t ResolveVersions();
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_JOB_MANAGER_H_
