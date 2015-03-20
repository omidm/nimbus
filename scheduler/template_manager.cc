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

#include "scheduler/template_manager.h"

using namespace nimbus; // NOLINT


TemplateManager::TemplateManager() {
  job_manager_ = NULL;
  worker_template_active_ = false;
}

TemplateManager::~TemplateManager() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  TemplateMap::iterator iter = template_map_.begin();
  for (; iter != template_map_.end(); ++iter) {
    delete iter->second;
  }
}

void TemplateManager::set_job_manager(JobManager* job_manager) {
  job_manager_ = job_manager;
}

void TemplateManager::set_id_maker(IDMaker *id_maker) {
  id_maker_ = id_maker;
}

void TemplateManager::set_worker_template_active(bool flag) {
  worker_template_active_ = flag;
}


bool TemplateManager::DetectNewTemplate(const std::string& template_name) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  TemplateMap::iterator iter = template_map_.find(template_name);
  if (iter == template_map_.end()) {
    template_map_[template_name] = new TemplateEntry(template_name,
                                                     worker_template_active_);
  } else {
    if (!iter->second->finalized()) {
      if (iter->second->CleanPartiallyFilledTemplate()) {
        dbg(DBG_WARN, "WARNING: cleaned up partially filled template.\n");
      } else {
        dbg(DBG_ERROR, "ERROR: could not clean partially filled template!\n");
        return false;
      }
    } else {
      dbg(DBG_ERROR, "ERROR: template already has been detected and finalized!\n");
      return false;
    }
  }

  return true;
}

bool TemplateManager::FinalizeNewTemplate(const std::string& template_name) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  TemplateMap::iterator iter = template_map_.find(template_name);
  if (iter != template_map_.end()) {
    return iter->second->Finalize();
  } else {
    dbg(DBG_ERROR, "ERROR: template has NOT been detected to get finalized!\n");
    return false;
  }

  return true;
}

bool TemplateManager::InstantiateTemplate(const std::string& template_name,
                                          const std::vector<job_id_t>& inner_job_ids,
                                          const std::vector<job_id_t>& outer_job_ids,
                                          const std::vector<Parameter>& parameters,
                                          const job_id_t& parent_job_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  TemplateMap::iterator iter = template_map_.find(template_name);
  if (iter != template_map_.end()) {
    if (job_manager_ == NULL) {
      dbg(DBG_ERROR, "ERROR: job manager is not initialized!\n");
      return false;
    }
    if (iter->second->finalized()) {
      return iter->second->Instantiate(job_manager_,
                                      inner_job_ids,
                                      outer_job_ids,
                                      parameters,
                                      parent_job_id);
    } else {
      dbg(DBG_ERROR, "ERROR: template has NOT been finalized to get instantiated!\n");
      return false;
    }
  } else {
    dbg(DBG_ERROR, "ERROR: template has NOT been detected to get instantiated!\n");
    return false;
  }

  return true;
}


bool TemplateManager::GetComplexJobEntryForTemplate(ComplexJobEntry*& complex_job,
                                                    const std::string& template_name,
                                                    const job_id_t& parent_job_id,
                                                    const std::vector<job_id_t>& inner_job_ids,
                                                    const std::vector<job_id_t>& outer_job_ids,
                                                    const std::vector<Parameter>& parameters) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  TemplateMap::iterator iter = template_map_.find(template_name);
  if (iter != template_map_.end()) {
    if (iter->second->finalized()) {
      std::vector<job_id_t> j;
      id_maker_->GetNewJobID(&j, 1);

      return iter->second->GetComplexJobEntry(complex_job,
                                              j[0],
                                              parent_job_id,
                                              inner_job_ids,
                                              outer_job_ids,
                                              parameters);
    } else {
      dbg(DBG_ERROR, "ERROR: template has NOT been finalized to get complex job!\n");
      return false;
    }
  } else {
    dbg(DBG_ERROR, "ERROR: template has NOT been detected to get complex job!\n");
    return false;
  }

  return true;
}

bool TemplateManager::SetBaseVersionMapForTemplate(const std::string& template_name,
                                                   boost::shared_ptr<VersionMap> vmap_base) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  TemplateMap::iterator iter = template_map_.find(template_name);
  if (iter != template_map_.end()) {
      iter->second->set_vmap_base(vmap_base);
  } else {
    dbg(DBG_ERROR, "ERROR: template has NOT been detected to get complex job!\n");
    return false;
  }

  return true;
}

TemplateJobEntry*
TemplateManager::AddComputeJobToTemplate(const std::string& template_name,
                                         const std::string& job_name,
                                         const job_id_t& job_id,
                                         const IDSet<logical_data_id_t>& read_set,
                                         const IDSet<logical_data_id_t>& write_set,
                                         const IDSet<job_id_t>& before_set,
                                         const IDSet<job_id_t>& after_set,
                                         const job_id_t& parent_job_id,
                                         const job_id_t& future_job_id,
                                         const bool& sterile,
                                         const GeometricRegion& region) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  TemplateMap::iterator iter = template_map_.find(template_name);
  if (iter != template_map_.end()) {
    if (job_manager_ == NULL) {
      dbg(DBG_ERROR, "ERROR: job manager is not initialized!\n");
      return false;
    }
    if (!iter->second->finalized()) {
      return iter->second->AddComputeJob(job_name,
                                         job_id,
                                         read_set,
                                         write_set,
                                         before_set,
                                         after_set,
                                         parent_job_id,
                                         future_job_id,
                                         sterile,
                                         region);
    } else {
      dbg(DBG_ERROR, "ERROR: template has been finalized and cannot add compute job!\n");
      return NULL;
    }
  } else {
    dbg(DBG_ERROR, "ERROR: template has NOT been detected to add compute job!\n");
    return NULL;
  }
}

bool TemplateManager::AddExplicitCopyJobToTemplate() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  dbg(DBG_ERROR, "ERROR: explicit copy jobs from application are not supported yet!.\n");
  assert(false);
  return false;
}


