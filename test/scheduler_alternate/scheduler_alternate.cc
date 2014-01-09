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
  * This scheduler is written to alternate spawned jobs randomly among
  * registered workers of water simulation. It is guaranteed that after
  * convergence for each frame the write frame job whose name is defined as a
  * macro with tag WRITE_FRAME_NAME will be executed over the same worker
  * (first registered worker).This way, all the output frames are local to one
  * worker for sanity checks. The random seed can be changed by changing the
  * SEED_ macro. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./scheduler_alternate.h"

#define SEED_ 123
#define WRITE_FRAME_NAME "write_frame"

SchedulerAlternate::SchedulerAlternate(unsigned int p)
: Scheduler(p) {
  seed_ = SEED_;
}

bool SchedulerAlternate::GetWorkerToAssignJob(JobEntry* job, SchedulerWorker*& worker) {
  size_t worker_num = server_->worker_num();

  if (worker_num < 1) {
    dbg(DBG_SCHED, "ERROR: there is no worker in scheduler worker list for job assignment");
    return false;
  }

  worker_id_t w_id;
  if (job->job_name() == WRITE_FRAME_NAME) {
    w_id = 1;
  } else {
    w_id  = (rand_r(&seed_) % worker_num) + 1;
  }

  return server_->GetSchedulerWorkerById(worker, w_id);
}
