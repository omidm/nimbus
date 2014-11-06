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
  * This is TemplateEntry module to hold and instantiate the templates.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_TEMPLATE_ENTRY_H_
#define NIMBUS_SCHEDULER_TEMPLATE_ENTRY_H_

#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "shared/log.h"
#include "scheduler/job_manager.h"

namespace nimbus {
class TemplateEntry {
  public:
    TemplateEntry();
    ~TemplateEntry();

    bool finalaized();

    bool Finalize();

    bool CleanPartiallyFilledTemplate();

    bool Instantiate(JobManager *job_manager,
                     const std::vector<job_id_t>& inner_job_ids,
                     const std::vector<job_id_t>& outer_job_ids,
                     const std::vector<Parameter>& parameters,
                     const job_id_t& parent_job_id);

    bool AddComputeJob(const std::string& job_name,
                       const job_id_t& job_id,
                       const IDSet<logical_data_id_t>& read_set,
                       const IDSet<logical_data_id_t>& write_set,
                       const IDSet<job_id_t>& before_set,
                       const IDSet<job_id_t>& after_set,
                       const job_id_t& parent_job_id,
                       const job_id_t& future_job_id,
                       const bool& sterile,
                       const GeometricRegion& region);

    bool AddExplicitCopyJob();

  private:
    typedef std::vector<job_id_t*> PtrList;
    typedef boost::unordered_set<job_id_t*> PtrSet;
    typedef boost::unordered_map<job_id_t, job_id_t*> PtrMap;

    class ComputeJobEntry {
      public:
        ComputeJobEntry(const std::string& job_name,
                        job_id_t* job_id_ptr,
                        const IDSet<logical_data_id_t>& read_set,
                        const IDSet<logical_data_id_t>& write_set,
                        const PtrSet& before_set_ptrs,
                        const PtrSet& after_set_ptrs,
                        job_id_t* future_job_id_ptr,
                        const bool& sterile,
                        const GeometricRegion& region)
          : job_name_(job_name),
            job_id_ptr_(job_id_ptr),
            read_set_(read_set),
            write_set_(write_set),
            before_set_ptrs_(before_set_ptrs),
            after_set_ptrs_(after_set_ptrs),
            future_job_id_ptr_(future_job_id_ptr),
            sterile_(sterile),
            region_(region) {}

        ~ComputeJobEntry() {}

        std::string job_name_;
        job_id_t* job_id_ptr_;
        IDSet<logical_data_id_t> read_set_;
        IDSet<logical_data_id_t> write_set_;
        PtrSet before_set_ptrs_;
        PtrSet after_set_ptrs_;
        job_id_t* future_job_id_ptr_;
        bool sterile_;
        GeometricRegion region_;
    };

    typedef std::vector<ComputeJobEntry> EntryList;

    bool finalized_;
    PtrList job_id_ptrs_;
    PtrMap job_id_ptrs_map_;
    EntryList entry_list_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_TEMPLATE_ENTRY_H_
