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
  * Job profile in the job graph of load balancer.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_JOB_PROFILE_H_
#define NIMBUS_SCHEDULER_JOB_PROFILE_H_

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
#include "scheduler/version_table.h"
#include "scheduler/version_map.h"
#include "scheduler/ancestor_chain.h"
#include "scheduler/meta_before_set.h"
#include "scheduler/logical_data_lineage.h"

namespace nimbus {

class JobProfile;
typedef std::map<job_id_t, JobProfile*> JobProfileMap;
typedef std::map<job_id_t, JobProfile*> JobProfileTable;
typedef std::list<JobProfile*> JobProfileList;

class JobProfile {
  public:
    JobProfile();

    JobProfile(
        const JobType& job_type,
        const std::string& job_name,
        const job_id_t& job_id,
        const IDSet<job_id_t>& before_set,
        const job_id_t& parent_job_id,
        const bool& sterile);

    virtual ~JobProfile();

    JobType job_type();
    std::string job_name();
    job_id_t job_id();
    IDSet<job_id_t> before_set();
    job_id_t parent_job_id();
    bool sterile();
    bool done();
    IDSet<job_id_t>* before_set_p();
    const IDSet<job_id_t>* before_set_p() const;

    void set_job_type(JobType job_type);
    void set_job_name(std::string job_name);
    void set_job_id(job_id_t job_id);
    void set_before_set(IDSet<job_id_t> before_set);
    void set_parent_job_id(job_id_t parent_job_id);
    void set_sterile(bool flag);
    void set_done(bool flag);

    bool GetPhysicalReadSet(IDSet<physical_data_id_t>* set);
    bool GetPhysicalWriteSet(IDSet<physical_data_id_t>* set);

  private:
    JobType job_type_;
    std::string job_name_;
    job_id_t job_id_;
    IDSet<job_id_t> before_set_;
    job_id_t parent_job_id_;
    bool sterile_;
    bool done_;

    void Initialize();
};


}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_JOB_PROFILE_H_


