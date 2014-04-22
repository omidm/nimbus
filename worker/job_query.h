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
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_WORKER_JOB_QUERY_H_
#define NIMBUS_WORKER_JOB_QUERY_H_

#include <list>
#include <map>
#include <string>
#include "shared/nimbus.h"
#include "worker/job.h"

namespace nimbus {

class JobQuery {
 public:
  explicit JobQuery(Job* job);
  ~JobQuery();
  /*
  bool SpawnComputeJob(
      const std::string& name, const job_id_t& id,
      const IDSet<logical_data_id_t>& read,
      const IDSet<logical_data_id_t>& write,
      const IDSet<job_id_t>& before,
      const IDSet<job_id_t>& after,
      const Parameter& params,
      const bool sterile = false) {
    return StageJob(name, id, read, write, before, after, params, sterile,
                    false);
  }
  bool SpawnLastJob(
      const std::string& name, const job_id_t& id,
      const IDSet<logical_data_id_t>& read,
      const IDSet<logical_data_id_t>& write,
      const IDSet<job_id_t>& before,
      const IDSet<job_id_t>& after,
      const Parameter& params,
      const bool sterile = false) {
    return StageJob(name, id, read, write, before, after, params, sterile,
                    true);
  }
  */
  bool StageJob(
      const std::string& name, const job_id_t& id,
      const IDSet<logical_data_id_t>& read,
      const IDSet<logical_data_id_t>& write,
      const Parameter& params,
      const bool sterile,
      const bool barrier = false);
  bool StageJob(
      const std::string& name, const job_id_t& id,
      const IDSet<logical_data_id_t>& read,
      const IDSet<logical_data_id_t>& write,
      const IDSet<job_id_t>& before,
      const IDSet<job_id_t>& after,
      const Parameter& params,
      const bool sterile = false,
      const bool barrier = false);
  bool CommitStagedJobs();
  bool CommitJob(const job_id_t& id);
  void GenerateDotFigure(const std::string& file_name);

 private:
  void Eliminate(IDSet<job_id_t>* before);
  bool RESOLVE_WAW;
  bool RESOLVE_WAR;
  bool DISABLE;

  Job* job_;
  struct OutstandingAccessors {
    OutstandingAccessors() : has_outstanding_writer(false) {}
    bool has_outstanding_writer;
    job_id_t outstanding_writer;
    std::list<job_id_t> outstanding_reader_list;
  };
  typedef std::map<logical_data_id_t, OutstandingAccessors>
      OutstandingAccessorsMap;
  OutstandingAccessorsMap outstanding_accessors_map_;

  struct JobEntry {
    std::string name;
    job_id_t id;
    IDSet<logical_data_id_t> read, write;
    IDSet<job_id_t> before;
    Parameter params;
    bool sterile;
  };
  std::list<JobEntry> staged_jobs_;

  struct ShortJobEntry {
    std::string name;
    job_id_t id;
    IDSet<job_id_t> before;
  };
  std::list<ShortJobEntry> query_results_;
};

}  // namespace nimbus

#endif  // NIMBUS_WORKER_JOB_QUERY_H_
