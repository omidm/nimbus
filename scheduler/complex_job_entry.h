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
  * Complex job entry in the job table of the job manager. This job contains a
  * group of compute or explicit copy jobs that are spawned with in a template.
  * The meta data calculated for the template including job dependencies and
  * versioning information is precomputed and is accessible by the complex job.
  * The idea is no avoid expanding the jobs within the template spawning. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_COMPLEX_JOB_ENTRY_H_
#define NIMBUS_SCHEDULER_COMPLEX_JOB_ENTRY_H_

#include <boost/unordered_map.hpp>
#include <vector>
#include <string>
#include <set>
#include <list>
#include <utility>
#include <map>
#include "shared/nimbus_types.h"
#include "shared/idset.h"
#include "shared/parameter.h"
#include "shared/geometric_region.h"
#include "scheduler/job_entry.h"
#include "scheduler/version_map.h"
#include "scheduler/data_manager.h"
#include "scheduler/meta_before_set.h"
#include "scheduler/logical_data_lineage.h"

namespace nimbus {

class TemplateEntry;

class ComplexJobEntry : public JobEntry {
  public:
    ComplexJobEntry();

    ComplexJobEntry(const job_id_t& job_id,
                    const job_id_t& parent_job_id,
                    TemplateEntry* template_entry,
                    const std::vector<job_id_t>& inner_job_ids,
                    const std::vector<job_id_t>& outer_job_ids,
                    const std::vector<Parameter>& parameters);

    virtual ~ComplexJobEntry();

    TemplateEntry* template_entry() const;
    std::vector<job_id_t> inner_job_ids() const;
    std::vector<job_id_t> outer_job_ids() const;
    std::vector<Parameter> parameters() const;
    const std::vector<job_id_t>* inner_job_ids_p() const;
    const std::vector<job_id_t>* outer_job_ids_p() const;
    const std::vector<Parameter>* parameters_p() const;


    size_t GetJobForAssignment(JobEntryList* list, size_t max_num);
    void MarkJobAssigned(job_id_t job_id);
    void MarkJobDone(job_id_t job_id);

  private:
    TemplateEntry* template_entry_;
    std::vector<job_id_t> inner_job_ids_;
    std::vector<job_id_t> outer_job_ids_;
    std::vector<Parameter> parameters_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_COMPLEX_JOB_ENTRY_H_

