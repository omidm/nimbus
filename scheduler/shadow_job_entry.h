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
  * Shadow job entry is a version of job entry that contanis only required
  * information to assign the job properly. It includes version maps and
  * read/write set. This data id=s read from the main job entry and shared
  * among multiple shadoes. This way the creation of shadow jobs is very cost
  * effective, cause only a couple od pointers are passed to each instance.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_SHADOW_JOB_ENTRY_H_
#define NIMBUS_SCHEDULER_SHADOW_JOB_ENTRY_H_

#include <boost/unordered_map.hpp>
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
#include "shared/geometric_region.h"
#include "scheduler/job_entry.h"
#include "scheduler/version_map.h"
#include "scheduler/scheduler_worker.h"
#include "scheduler/meta_before_set.h"
#include "scheduler/logical_data_lineage.h"

namespace nimbus {

class ShadowJobEntry : public JobEntry {
  public:
    ShadowJobEntry();

    explicit ShadowJobEntry(const job_id_t& job_id);

    ShadowJobEntry(const std::string& job_name,
                   const job_id_t& job_id,
                   const IDSet<logical_data_id_t>* read_set_p,
                   const IDSet<logical_data_id_t>* write_set_p,
                   const IDSet<logical_data_id_t>* union_set_p,
                   boost::shared_ptr<VersionMap> vmap_read,
                   boost::shared_ptr<VersionMap> vmap_write,
                   const job_id_t& parent_job_id,
                   const job_id_t& future_job_id,
                   const bool& sterile,
                   const GeometricRegion& region,
                   const Parameter& params);

    virtual ~ShadowJobEntry();

    virtual const IDSet<logical_data_id_t>* read_set_p() const;
    virtual const IDSet<logical_data_id_t>* write_set_p() const;
    virtual const IDSet<logical_data_id_t>* union_set_p() const;

    virtual void set_read_set_p(const IDSet<logical_data_id_t>* read_set_p);
    virtual void set_write_set_p(const IDSet<logical_data_id_t>* write_set_p);
    virtual void set_union_set_p(const IDSet<logical_data_id_t>* union_set_p);

  private:
    const IDSet<logical_data_id_t>* read_set_p_;
    const IDSet<logical_data_id_t>* write_set_p_;
    const IDSet<logical_data_id_t>* union_set_p_;
};


}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_SHADOW_JOB_ENTRY_H_


