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
  * This scheduler is written to alternate the iterations for each frame of
  * water simulation between two workers. It is guaranteed that after
  * convergence for each frame the write_frame job will be executed over the
  * same worker, so that all the output frames are local to one worker for
  * sanity checks.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_TEST_SCHEDULER_ALTERNATE_COARSE_SCHEDULER_ALTERNATE_COARSE_H_
#define NIMBUS_TEST_SCHEDULER_ALTERNATE_COARSE_SCHEDULER_ALTERNATE_COARSE_H_

#define DEBUG_MODE

#include <boost/thread.hpp>
#include <stdlib.h>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "shared/scheduler_server.h"
#include "shared/cluster.h"
#include "shared/parser.h"
#include "scheduler/scheduler.h"

class SchedulerAlternateCoarse : public Scheduler {
  public:
    explicit SchedulerAlternateCoarse(unsigned int listening_port);

    virtual bool GetWorkerToAssignJob(JobEntry* job, SchedulerWorker*& worker);

  private:
    unsigned int seed_;
};

#endif  // NIMBUS_TEST_SCHEDULER_ALTERNATE_COARSE_SCHEDULER_ALTERNATE_COARSE_H_
