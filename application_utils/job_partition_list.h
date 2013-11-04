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

#ifndef NIMBUS_APPLICATION_UTILS_JOB_PARTITION_LIST_H_
#define NIMBUS_APPLICATION_UTILS_JOB_PARTITION_LIST_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"

namespace nimbus {

    typedef int app_job_key_t;

    struct LocationID {
        GeometricRegion *region;
        app_job_key_t key;
    };

    typedef std::list<LocationID*> LocationIDList;
    typedef std::set<app_job_key_t> JobKeySet;

    class JobPartitionList {
        public:
            JobPartitionList();
            virtual ~JobPartitionList();
            virtual void AddJob(std::string jobtype, GeometricRegion *region);
            virtual JobKeySet* GetIntersectingJobs(std::string jobtype);
        private:
            app_job_key_t available_key_;
            std::map<std::string, LocationIDList*> map_;
    };

}  // namespace nimbus

#endif  // NIMBUS_APPLICATION_UTILS_JOB_PARTITION_LIST_H_
