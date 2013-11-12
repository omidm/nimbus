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
  * Nimbus job abstraction from scheduler point of view. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_SCHEDULER_JOB_H_
#define NIMBUS_SCHEDULER_SCHEDULER_JOB_H_

#include <vector>
#include <string>
#include <set>
#include <list>
#include <map>
#include "worker/data.h"
#include "shared/idset.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class SchedulerJob;
typedef std::map<int, SchedulerJob*> SchedulerJobMap;
typedef std::list<SchedulerJob*> SchedulerJobList;
typedef std::vector<Data*> DataArray;

class SchedulerJob {
 public:
  SchedulerJob();
  SchedulerJob(job_id_t id, app_id_t app_id, JobType type);
  virtual ~SchedulerJob();

  uint64_t id();
  void set_id(job_id_t id);


 private:
  job_id_t id_;
  app_id_t app_id_;
  JobType type_;
  IDSet<physical_data_id_t> read_set_;
  IDSet<physical_data_id_t> write_set_;
  IDSet<job_id_t> before_set_;
  IDSet<job_id_t> after_set_;
  std::string parameters_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_SCHEDULER_JOB_H_


