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
  log_.set_file_name("log_template_manager");
}

TemplateManager::~TemplateManager() {
}

void TemplateManager::set_job_manager(JobManager* job_manager) {
  job_manager_ = job_manager;
}

bool TemplateManager::DetectNewTemplate(const std::string& template_name) {
  // TODO(omidm): Implement!
  return false;
}

bool TemplateManager::FinalizeNewTemplate(const std::string& template_name) {
  // TODO(omidm): Implement!
  return false;
}

bool TemplateManager::InstantiateTemplate(const std::string& template_name,
                                          const std::vector<job_id_t>& inner_job_ids,
                                          const std::vector<job_id_t>& outer_job_ids,
                                          const std::vector<Parameter>& parameters,
                                          const job_id_t& parent_job_id) {
  // TODO(omidm): Implement!
  return false;
}


bool TemplateManager::AddComputeJobToTemplate(const std::string& template_name,
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
  // TODO(omidm): Implement!
  return false;
}

bool TemplateManager::AddExplicitCopyJobToTemplate() {
  dbg(DBG_ERROR, "ERROR: explicit copy jobs from application are not supported yet!.\n");
  exit(-1);
  return false;
}




