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
  * Scheduler Version Entry. It holds the meta data for each version of the
  * data containing logical data id, version number, and a pointer to a job
  * related to the version, and the type of relation meaning whether job needs
  * the version as input or dumps it as output.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_VERSION_ENTRY_H_
#define NIMBUS_SRC_SCHEDULER_VERSION_ENTRY_H_

#include <boost/thread.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <set>
#include <vector>
#include <list>
#include <algorithm>
#include "src/shared/nimbus_types.h"
#include "src/shared/dbg.h"
#include "src/scheduler/job_entry.h"
#include "src/scheduler/shadow_job_entry.h"
#include "src/scheduler/complex_job_entry.h"
#include "src/scheduler/logical_data_lineage.h"

namespace nimbus {

class VersionEntry {
  public:
    typedef boost::unordered_set<JobEntry*> Bucket;
    typedef Bucket::iterator BucketIter;
    typedef boost::unordered_map<data_version_t, Bucket*> Index;
    typedef Index::iterator IndexIter;

    explicit VersionEntry(logical_data_id_t ldid);
    VersionEntry(const VersionEntry& ve);
    virtual ~VersionEntry();

    VersionEntry& operator= (const VersionEntry& right);

    void InitializeLdl(
        const job_id_t& job_id,
        const job_depth_t& job_depth);

    bool AddJobEntryReader(JobEntry *job);

    bool AddJobEntryWriter(JobEntry *job);

    bool RemoveJobEntry(JobEntry *job);

    size_t GetJobsNeedVersion(
        JobEntryList* list, data_version_t version, bool append = false);

    bool LookUpVersion(JobEntry *job,
                       data_version_t *version);

    bool LookUpVersionByMetaBeforeSet(
                       boost::shared_ptr<MetaBeforeSet> mbs,
                       data_version_t *version);

    bool CleanLdl(const IDSet<job_id_t>& snap_shot_);

    void Reinitialize();

    bool InsertCheckPointLdlEntry(
        const job_id_t& job_id,
        const data_version_t& version,
        const job_depth_t& job_depth);

    bool is_empty();

  private:
    logical_data_id_t ldid_;
    Bucket seen_parent_jobs_;
    Bucket pending_reader_jobs_;
    Bucket pending_writer_jobs_;
    Index index_;
    LogicalDataLineage ldl_;
    boost::recursive_mutex mutex_;

    bool UpdateLdl();

    void ReinitializeLdl();

    void ClearIndex();
};

typedef std::list<VersionEntry> VersionEntryList;

}  // namespace nimbus

#endif  // NIMBUS_SRC_SCHEDULER_VERSION_ENTRY_H_

