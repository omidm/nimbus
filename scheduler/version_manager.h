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
  * Scheduler Version Manager. This module serves the job manager by keeping
  * track of the version numbers of each logical data object that are needed by
  * the jobs in the job graph.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_VERSION_MANAGER_H_
#define NIMBUS_SCHEDULER_VERSION_MANAGER_H_

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <utility>
#include <list>
#include <map>
#include <set>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "scheduler/job_entry.h"
#include "scheduler/complex_job_entry.h"
#include "scheduler/template_job_entry.h"
#include "scheduler/template_entry.h"
#include "scheduler/binding_template.h"
#include "scheduler/version_entry.h"
#include "shared/logical_data_object.h"

namespace nimbus {


typedef std::pair<logical_data_id_t, data_version_t> VersionedLogicalData;

class VersionManager {
  public:
    typedef std::pair<logical_data_id_t, data_version_t> VLD;
    typedef boost::unordered_map<logical_data_id_t, VersionEntry*> Index;
    typedef std::map<job_id_t, JobEntry*> ParentMap;
    typedef std::map<job_id_t, ComplexJobEntry*> ComplexJobMap;
    typedef std::map<job_id_t, counter_t> ChildCounter;

    VersionManager();
    virtual ~VersionManager();

    void set_snap_shot_rate(size_t rate);

    bool AddJobEntry(JobEntry *job);

    bool AddComplexJobEntry(ComplexJobEntry *complex_job);

    size_t GetJobsNeedDataVersion(
        JobEntryList* list, VLD vld);

    bool NotifyJobDone(JobEntry *job);

    bool RemoveJobEntry(JobEntry* job);

    bool ResolveJobDataVersions(JobEntry *job);

    bool ResolveEntireContextForJob(JobEntry *job);

    bool ResolveJobDataVersionsForPattern(JobEntry *job,
                  const BindingTemplate::PatternList* patterns);

    bool MemoizeVersionsForTemplate(JobEntry *job);

    bool DefineData(
        const logical_data_id_t ldid,
        const job_id_t& job_id,
        const job_depth_t& job_depth);

    bool CleanUp();

    void Reinitialize(const JobEntryList *list);

    void set_ldo_map_p(const LdoMap* ldo_map_p);

  private:
    Log log_;
    Index index_;
    bool snap_shot_pending_;
    size_t snap_shot_rate_;
    size_t non_sterile_counter_;
    IDSet<job_id_t> snap_shot_;
    ParentMap parent_map_;
    ComplexJobMap complex_jobs_;
    ChildCounter child_counter_;
    boost::recursive_mutex snap_shot_mutex_;
    const LdoMap* ldo_map_p_;

    bool LookUpVersion(JobEntry *job,
                       logical_data_id_t ldid,
                       data_version_t *version);

    bool InsertComplexJobInLdl(ComplexJobEntry *job);

    bool CreateCheckPoint(JobEntry *job);

    bool DetectNewJob(JobEntry *job);

    bool DetectNewComplexJob(ComplexJobEntry *xj);

    bool DetectVersionedJob(JobEntry *job);

    bool DetectDoneJob(JobEntry *job);

    bool GetSnapShot();

    bool InsertCheckPointLdlEntry(
        const logical_data_id_t ldid,
        const job_id_t& job_id,
        const data_version_t& version,
        const job_depth_t& job_depth);
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_VERSION_MANAGER_H_
