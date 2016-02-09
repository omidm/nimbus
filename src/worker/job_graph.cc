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
  * Nimbus job graph base class. JobGraph graph defines a template of job
  * dependencies. Multiple instances of the template could be spawned only by
  * filling in the new job ids and specific parameters for each instance of the
  * job in the job graph.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <time.h>
#include "src/worker/job_graph.h"
#include "src/worker/application.h"

using namespace nimbus; // NOLINT

JobGraph::JobGraph(Application *application,
                   std::string job_graph_name,
                   size_t inner_job_num,
                   size_t outer_job_num) {
  application_ = application;
  job_graph_name_ = job_graph_name;
  inner_job_num_ = inner_job_num;
  outer_job_num_ = outer_job_num;
  InitializeIDLists();
}

void JobGraph::InitializeIDLists() {
  for (size_t i = 0; i < inner_job_num_; ++i) {
    inner_job_ids_.push_back(i);
  }
  for (size_t i = 0; i < outer_job_num_; ++i) {
    outer_job_ids_.push_back(i);
  }
}

JobGraph::~JobGraph() {
}

bool JobGraph::AddComputeJob(const std::string& name,
                             const job_id_t& id,
                             const IDSet<logical_data_id_t>& read,
                             const IDSet<logical_data_id_t>& write,
                             const IDSet<job_id_t>& before,
                             const IDSet<job_id_t>& after,
                             const bool& sterile,
                             const GeometricRegion& region,
                             const job_id_t& future_job_id) {
  application_->AddComputeJobToJobGraph(name,
                                        id,
                                        read,
                                        write,
                                        before,
                                        after,
                                        sterile,
                                        region,
                                        future_job_id,
                                        job_graph_name_);
  return true;
}

bool JobGraph::AddCopyJob(const job_id_t& id,
                          const logical_data_id_t& from_logical_id,
                          const logical_data_id_t& to_logical_id,
                          const IDSet<job_id_t>& before,
                          const IDSet<job_id_t>& after) {
  application_->AddCopyJobToJobGraph(id,
                                     from_logical_id,
                                     to_logical_id,
                                     before,
                                     after,
                                     job_graph_name_);
  return true;
}



job_id_t JobGraph::GetInnerJobID(size_t index) {
  assert(index < inner_job_num_);
  return inner_job_ids_[index];
}

job_id_t JobGraph::GetOuterJobID(size_t index) {
  assert(index < outer_job_num_);
  return outer_job_ids_[index];
}

bool JobGraph::GetPartition(partition_id_t id, GeometricRegion* r) const {
  return application_->GetPartition(id, r);
}

const LogicalDataObject* JobGraph::GetLogicalObject(logical_data_id_t id) const {
  return application_->GetLogicalObject(id);
}

int JobGraph::GetCoveredLogicalObjects(CLdoVector* result,
                                       const std::string& variable,
                                       const GeometricRegion* r) {
  return application_->GetCoveredLogicalObjects(result, variable, r);
}

int JobGraph::GetAdjacentLogicalObjects(CLdoVector* result,
                                        const std::string& variable,
                                        const GeometricRegion* r) {
  return application_->GetAdjacentLogicalObjects(result, variable, r);
}

int JobGraph::GetIntersectingLogicalObjects(CLdoVector* result,
                                            const std::string& variable,
                                            const GeometricRegion* r) {
  return application_->GetIntersectingLogicalObjects(result, variable, r);
}

Application* JobGraph::application() const {
  return application_;
}

size_t JobGraph::inner_job_num() const {
  return inner_job_num_;
}

size_t JobGraph::outer_job_num() const {
  return outer_job_num_;
}



