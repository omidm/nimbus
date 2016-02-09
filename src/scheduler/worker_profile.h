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
  * Worker profile used in load balancer.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_WORKER_PROFILE_H_
#define NIMBUS_SRC_SCHEDULER_WORKER_PROFILE_H_

#include <boost/thread.hpp>
#include <vector>
#include <string>
#include <set>
#include <list>
#include <utility>
#include <map>
#include "src/shared/idset.h"
#include "src/shared/chronometer.h"
#include "src/shared/nimbus_types.h"

namespace nimbus {

class WorkerProfile;
typedef std::map<job_id_t, WorkerProfile*> WorkerProfileMap;
typedef std::map<job_id_t, WorkerProfile*> WorkerProfileTable;
typedef std::list<WorkerProfile*> WorkerProfileList;

class WorkerProfile {
  public:
    explicit WorkerProfile(worker_id_t worker_id);
    virtual ~WorkerProfile();

    worker_id_t worker_id();
    double idle_time();
    double active_time();
    double blocked_time();
    IDSet<job_id_t>* ready_jobs();
    IDSet<job_id_t>* blocked_jobs();

    void ResetTimers();

    void AddReadyJob(job_id_t job_id);

    void AddBlockedJob(job_id_t job_id);

    void NotifyJobDone(job_id_t job_id);

    std::string PrintStats();

  private:
    worker_id_t worker_id_;
    IDSet<job_id_t> *ready_jobs_;
    IDSet<job_id_t> *blocked_jobs_;

    Chronometer idle_timer_;
    Chronometer active_timer_;
    Chronometer blocked_timer_;

    boost::recursive_mutex mutex_;
};


}  // namespace nimbus
#endif  // NIMBUS_SRC_SCHEDULER_WORKER_PROFILE_H_


