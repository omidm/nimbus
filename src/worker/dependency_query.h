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
  * This is the DependencyQuery module to load before set for jobs. It is
  * written based on the code Hang Qu had provided as JobQuery module. It does
  * not have all the optimizations from the initial code and implements a
  * simpler algorithm. The reason is that we need to split different
  * functionalities, to make the code more modular and all optimizations were
  * interleaved in a complex way.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  * Modified by: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_WORKER_DEPENDENCY_QUERY_H_
#define NIMBUS_SRC_WORKER_DEPENDENCY_QUERY_H_

#include <time.h>
#include <boost/unordered_map.hpp>
#include <list>
#include <string>
#include <vector>
#include "src/shared/nimbus_types.h"
#include "src/shared/geometric_region.h"

namespace nimbus {

class DependencyQuery {
 public:
  DependencyQuery();
  ~DependencyQuery();

  // Read and write set will be cleared, due to performance consideration.
  bool StageJobAndLoadBeforeSet(IDSet<job_id_t> *before_set,
                                const std::string& name,
                                const job_id_t& id,
                                const IDSet<logical_data_id_t>& read,
                                const IDSet<logical_data_id_t>& write,
                                const IDSet<logical_data_id_t>& scratch,
                                const IDSet<logical_data_id_t>& reduce,
                                const bool barrier = false);

  bool MarkEndOfStage();

  void PrintTimeProfile();

  void GenerateDotFigure(const std::string& file_name);





 private:
  void Prune(IDSet<job_id_t>* before);

  typedef int64_t StageId;
  typedef int64_t IndexId;

  struct Job {
    job_id_t id;
    std::string name;
    IDSet<logical_data_id_t> read;
    IDSet<logical_data_id_t> write;
    IDSet<logical_data_id_t> scratch;
    IDSet<logical_data_id_t> reduce;
    IDSet<job_id_t> before;
    StageId stage_id;
  };

  struct LastBarrier {
    LastBarrier() : valid(false) {}
    bool valid;
    job_id_t id;
    StageId stage_id;
  };

  typedef std::vector<Job> JobList;
  typedef boost::unordered_map<logical_data_id_t, job_id_t> LastWriterMap;
  typedef std::vector<job_id_t> PartialWriters;
  typedef boost::unordered_map<logical_data_id_t, PartialWriters> PartialWritersMap;

  StageId stage_id_counter_;
  LastBarrier last_barrier_;
  JobList jobs_;
  JobList staged_jobs_;
  LastWriterMap last_writers_;
  PartialWritersMap *last_stage_partial_writers_;
  PartialWritersMap *current_stage_partial_writers_;
  boost::unordered_map<job_id_t, IndexId> job_id_to_index_;

  int64_t total_job_;
  double pruning_time_;
  double insertion_time_;
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_WORKER_DEPENDENCY_QUERY_H_
