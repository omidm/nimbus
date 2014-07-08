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
#include <utility>
#include <list>
#include <map>
#include <set>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "scheduler/job_entry.h"
#include "scheduler/version_entry.h"
#include "shared/logical_data_object.h"

namespace nimbus {


typedef std::pair<logical_data_id_t, data_version_t> VersionedLogicalData;

class VersionManager {
  public:
    typedef std::pair<logical_data_id_t, data_version_t> VLD;
    typedef boost::unordered_map<logical_data_id_t, VersionEntry*> Index;
    typedef Index::iterator IndexIter;

    VersionManager();
    virtual ~VersionManager();

    bool AddJobEntry(JobEntry *job);

    size_t GetJobsNeedDataVersion(
        JobEntryList* list, VLD vld);

    bool RemoveJobEntry(JobEntry* job);

    void set_ldo_map_p(const std::map<logical_data_id_t, LogicalDataObject*>* ldo_map_p);

  private:
    Index index_;
    const std::map<logical_data_id_t, LogicalDataObject*>* ldo_map_p_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_VERSION_MANAGER_H_
