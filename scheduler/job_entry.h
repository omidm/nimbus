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
  * Job entry in the job table of the job manager. Each entry holds the
  * meta data of the compute and copy jobs received by the scheduler.   
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_JOB_ENTRY_H_
#define NIMBUS_SCHEDULER_JOB_ENTRY_H_

#include <vector>
#include <string>
#include <set>
#include <list>
#include <utility>
#include <map>
#include "worker/data.h"
#include "shared/idset.h"
#include "shared/parameter.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class JobEntry;
typedef std::map<job_id_t, JobEntry*> JobEntryMap;
typedef std::map<job_id_t, JobEntry*> JobEntryTable;
typedef std::list<JobEntry*> JobEntryList;
typedef std::vector<Data*> DataArray;

class JobEntry {
  public:
    typedef std::pair<logical_data_id_t, data_version_t> VersionedLogicalData;
    typedef std::map<logical_data_id_t, VersionedLogicalData> VersionTable;

    JobEntry();
    JobEntry(const JobType& job_type,
        const std::string& job_name,
        const job_id_t& job_id,
        const IDSet<logical_data_id_t>& read_set,
        const IDSet<logical_data_id_t>& write_set,
        const IDSet<job_id_t>& before_set,
        const IDSet<job_id_t>& after_set,
        const job_id_t& parent_job_id,
        const Parameter& params);
    virtual ~JobEntry();

    JobType job_type();
    std::string job_name();
    job_id_t job_id();
    IDSet<logical_data_id_t> read_set();
    IDSet<logical_data_id_t> write_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    job_id_t parent_job_id();
    Parameter params();
    VersionTable version_table();

  private:
    JobType job_type_;
    std::string job_name_;
    job_id_t job_id_;
    IDSet<physical_data_id_t> read_set_;
    IDSet<physical_data_id_t> write_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    job_id_t parent_job_id_;
    Parameter params_;

    VersionTable version_table_;

    bool versioned;
    bool assigned;
    bool done;
};

typedef std::map<job_id_t, JobEntry*> JobEntryTable;

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_JOB_ENTRY_H_


