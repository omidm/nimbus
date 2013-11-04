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
  * A mapping from a job type to the partitions that we want the job to run
  * over, with corresponding application-specific unique keys.
  */

#include <list>
#include <map>
#include <string>
#include "application_utils/job_partition_list.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"

namespace nimbus {

    JobPartitionList::JobPartitionList() {
        available_key_ = 0;
    }

    JobPartitionList::~JobPartitionList() {}

    void JobPartitionList::AddJob(
            std::string jobtype,
            GeometricRegion *region) {
        LocationID *job_loc_id = new LocationID();
        job_loc_id->region = region;
        job_loc_id->key = available_key_++;
        LocationIDList *jobtype_list;
        if (map_.find(jobtype) == map_.end()) {
            jobtype_list = new LocationIDList();
            map_[jobtype] = jobtype_list;
        } else {
            jobtype_list = map_[jobtype];
        }
        jobtype_list->push_back(job_loc_id);
    }

    JobKeySet* JobPartitionList::GetIntersectingJobs(
            std::string jobtype) {
        if (map_.find(jobtype) == map_.end()) {
            return NULL;
        } else {
            LocationIDList *jobtype_list = map_[jobtype];
            JobKeySet *job_keys = new JobKeySet();
            for (LocationIDList::iterator iter = jobtype_list->begin();
                    iter != jobtype_list->end(); ++iter) {
                job_keys->insert((*iter)->key);
            }
            return job_keys;
        }
    }

}  // namespace nimbus
