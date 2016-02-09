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
  * This is TemplateManager module. This module serves the controller by providing
  * facilities to detect, save, and instantiate templates in runtime. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_TEMPLATE_MANAGER_H_
#define NIMBUS_SRC_SCHEDULER_TEMPLATE_MANAGER_H_

#include <boost/thread.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include "src/shared/nimbus_types.h"
#include "src/shared/dbg.h"
#include "src/shared/log.h"
#include "src/scheduler/job_manager.h"
#include "src/scheduler/template_entry.h"

namespace nimbus {
class TemplateManager {
  public:
    typedef boost::unordered_map<std::string, TemplateEntry*> TemplateMap;
    TemplateManager();
    ~TemplateManager();

    void set_job_manager(JobManager* job_manager);
    void set_id_maker(IDMaker* id_maker);
    void set_worker_template_active(bool flag);
    void set_mega_rcr_job_active(bool flag);

    bool DetectNewTemplate(const std::string& template_name);

    bool FinalizeNewTemplate(const std::string& template_name);

    bool InstantiateTemplate(const std::string& template_name,
                             const std::vector<job_id_t>& inner_job_ids,
                             const std::vector<job_id_t>& outer_job_ids,
                             const std::vector<Parameter>& parameters,
                             const job_id_t& parent_job_id);

    bool GetComplexJobEntryForTemplate(ComplexJobEntry*& complex_job,
                                       const std::string& template_name,
                                       const job_id_t& parent_job_id,
                                       const std::vector<job_id_t>& inner_job_ids,
                                       const std::vector<job_id_t>& outer_job_ids,
                                       const std::vector<Parameter>& parameters);

    TemplateJobEntry* AddComputeJobToTemplate(const std::string& template_name,
                                              const std::string& job_name,
                                              const job_id_t& job_id,
                                              const IDSet<logical_data_id_t>& read_set,
                                              const IDSet<logical_data_id_t>& write_set,
                                              const IDSet<job_id_t>& before_set,
                                              const IDSet<job_id_t>& after_set,
                                              const job_id_t& parent_job_id,
                                              const job_id_t& future_job_id,
                                              const bool& sterile,
                                              const GeometricRegion& region);

    bool SetBaseVersionMapForTemplate(const std::string& template_name,
                                      boost::shared_ptr<VersionMap> vmap_base);

    bool AddExplicitCopyJobToTemplate();


  private:
    IDMaker *id_maker_;
    JobManager *job_manager_;
    bool worker_template_active_;
    bool mega_rcr_job_active_;
    TemplateMap template_map_;
    boost::mutex mutex_;
};

}  // namespace nimbus
#endif  // NIMBUS_SRC_SCHEDULER_TEMPLATE_MANAGER_H_
