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
  * Stores an index of all the job entries, including job id,
  * read set, write set of each job, and offers query service for jobs
  * that might cause IO conflict with given read set and write set.
  *
  * Read set and write set are specificed by "VariableRegionSet" for now.
  * May be changed in the future.
  *
  * Test Tag: N/A.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef JOB_INDEX_H_
#define JOB_INDEX_H_

#include <list>
#include "shared/idset.h"
#include "shared/nimbus_types.h"
#include "test/experimental_quh/variable_region_set.h"

namespace nimbus {

typedef IDSet<job_id_t> JobIdSet;

class JobIndex {
 public:
  JobIndex();
  virtual ~JobIndex();

  virtual int QueryJobEntry(
      const VariableRegionSet& read_set,
      const VariableRegionSet& write_set,
      JobIdSet* result);
  virtual bool AddJobEntry(
      const job_id_t job_id,
      const VariableRegionSet& read_set,
      const VariableRegionSet& write_set);
  virtual bool DeleteJobEntry(const job_id_t job_id);

  virtual int AllJobEntries(JobIdSet* result);

 private:
  struct JobEntry {
    job_id_t job_id;
    VariableRegionSet read_set;
    VariableRegionSet write_set;
  };
  typedef std::list<JobEntry*> JobEntryTable;
  JobEntryTable _job_table;
};

}  // namespace nimbus

#endif  // JOB_INDEX_H_
