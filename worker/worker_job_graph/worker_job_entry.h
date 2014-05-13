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

#ifndef NIMBUS_WORKER_WORKER_JOB_GRAPH_WORKER_JOB_ENTRY_H_
#define NIMBUS_WORKER_WORKER_JOB_GRAPH_WORKER_JOB_ENTRY_H_

#include <limits>
#include "shared/nimbus_types.h"

namespace nimbus {

class Job;

class WorkerJobEntry {
 public:
  enum State {
    INIT,
    CONTROL,  // For internal usage.
    PENDING,  // Job not received yet.
    PENDING_DATA_RECEIVED,  // Job not received yet, but data received.
    BLOCKED,  // Blocked by other jobs/IO event.
    READY,  // Ready to run.
    FINISH  // Finished, but has not been deleted.
  };

  WorkerJobEntry();
  WorkerJobEntry(const job_id_t job_id, Job* job, State state = PENDING);
  virtual ~WorkerJobEntry() {}

  job_id_t get_job_id() {
    return job_id_;
  }
  Job* get_job() {
    return job_;
  }
  State get_state() {
    return state_;
  }
  void set_job_id(job_id_t job_id) {
    job_id_ = job_id;
  }
  void set_job(Job* job) {
    job_ = job;
  }
  void set_state(State state) {
    state_ = state;
  }

 private:
  job_id_t job_id_;
  Job* job_;
  State state_;
};

}  // namespace nimbus
#endif  // NIMBUS_WORKER_WORKER_JOB_GRAPH_WORKER_JOB_ENTRY_H_
