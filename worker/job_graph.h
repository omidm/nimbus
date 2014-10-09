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
  * Nimbus job graph base class. Job graph defines a template of job
  * dependencies. Multiple instances of the template could be spawned only by
  * filling in the new job ids and specific parameters for each instance of the
  * job in the job graph.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_WORKER_JOB_GRAPH_H_
#define NIMBUS_WORKER_JOB_GRAPH_H_

#include <vector>
#include <string>
#include <map>
#include <list>
#include "shared/nimbus_types.h"
#include "shared/id.h"
#include "shared/idset.h"
#include "shared/Parameter.h"
#include "shared/geometric_region.h"
#include "shared/logical_data_object.h"

namespace nimbus {

class Application;
class JobGraph;
typedef std::list<JobGraph*> JobGraphList;
typedef std::map<std::string, JobGraph*> JobGraphTable;

class JobGraph {
  public:
    JobGraph(Application *application,
             std::string job_graph_name,
             size_t inner_job_num,
             size_t outer_job_num);
    virtual ~JobGraph();

    virtual void Define() = 0;

    bool AddComputeJob(const std::string& name,
                       const job_id_t& id,
                       const IDSet<logical_data_id_t>& read,
                       const IDSet<logical_data_id_t>& write,
                       const IDSet<job_id_t>& before,
                       const IDSet<job_id_t>& after,
                       const bool& sterile = false,
                       const GeometricRegion& region = GeometricRegion(),
                       const job_id_t& future_job_id = 0);

    bool AddCopyJob(const job_id_t& id,
                    const logical_data_id_t& from_logical_id,
                    const logical_data_id_t& to_logical_id,
                    const IDSet<job_id_t>& before,
                    const IDSet<job_id_t>& after);

    job_id_t GetInnerJobID(size_t index);
    job_id_t GetOuterJobID(size_t index);

    bool GetPartition(partition_id_t id, GeometricRegion* r) const;

    const LogicalDataObject* GetLogicalObject(logical_data_id_t id) const;

    int GetCoveredLogicalObjects(CLdoVector* result,
                                 const std::string& variable,
                                 const GeometricRegion* r);

    int GetAdjacentLogicalObjects(CLdoVector* result,
                                  const std::string& variable,
                                  const GeometricRegion* r);

    int GetIntersectingLogicalObjects(CLdoVector* result,
                                      const std::string& variable,
                                      const GeometricRegion* r);

    Application* application() const;
    std::string job_graph_name() const;
    size_t inner_job_num() const;
    size_t outer_job_num() const;

  private:
    Application* application_;
    std::string job_graph_name_;
    size_t inner_job_num_;
    size_t outer_job_num_;

    std::vector<job_id_t> inner_job_ids_;
    std::vector<job_id_t> outer_job_ids_;

    void InitializeIDLists();
};

}  // namespace nimbus
#endif  // NIMBUS_WORKER_JOB_GRAPH_H_


