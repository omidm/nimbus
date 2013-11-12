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
  * Simple Nimbus scheduler that is supposed to run the application over
  * multiple workers. This scheduler is not still smart enough to make the
  * dynamic decisions. It assigned each part of the simulation to a worker.
  * This test is intended to check the data exchange functionality and how the
  * system correctly enforces the correct program frlow.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_TEST_SCHEDULER_MULTI_SCHEDULER_MULTI_H_
#define NIMBUS_TEST_SCHEDULER_MULTI_SCHEDULER_MULTI_H_

#define DEBUG_MODE

#include <boost/thread.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "shared/scheduler_server.h"
#include "shared/cluster.h"
#include "shared/parser.h"
#include "scheduler/scheduler.h"

#define WORKER_NUM 2

class SimpleScheduler : public Scheduler {
  public:
    explicit SimpleScheduler(unsigned int listening_port);

    virtual  void SchedulerCoreProcessor();
    virtual void ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm);
    virtual void ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm);
    virtual void ProcessDefineDataCommand(DefineDataCommand* cm);
    virtual void ProcessJobDoneCommand(JobDoneCommand* cm);

  private:
    SchedulerCommandList pending_compute_jobs_;
    SchedulerCommandList pending_copy_jobs_;
    std::map<job_id_t, physical_data_id_t> job_data_map_;
    std::map<physical_data_id_t, std::pair<worker_id_t, bool> > create_data_;
};

#endif  // NIMBUS_TEST_SCHEDULER_MULTI_SCHEDULER_MULTI_H_
