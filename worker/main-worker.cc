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
  * A Nimbus worker. 
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <cstdio>

#include "application/fluid-simulation/public/fluid-simulation-app.h"
#include "worker/job-worker-factory.h"
#include "worker/job-worker-interface.h"

// [TODO] Ugly code... Include class in main file. Reorgnize after discussion.
class InitJob : public JobWorkerInterface {
 public:
  InitJob() : argc_(0), argv_(NULL) {}
  InitJob(int argc, char **argv) : argc_(argc), argv_(argv) {}
  virtual void Run() {
    // [TODO] Run could be a non-blocking function, adding the job to the local
    // worker manager.
    main_job(argc_, argv_);
  }
 private:
  int argc_;
  char **argv_;
};
class ExecuteJob : public JobWorkerInterface {
 public:
  ExecuteJob() {}
  virtual void Run() {
    run_job();
  }
};

int main(int argc, char **argv) {
  JobWorkerFactory job_worker_factory;
  // [TODO] Parameter should be passed through run function as a parameter blob.
  InitJob init_job_proto(argc, argv);
  ExecuteJob execute_job_proto;
  job_worker_factory.RegisterJob(&init_job_proto, 0);
  job_worker_factory.RegisterJob(&execute_job_proto, 1);

  // [TODO} Maybe scheduler should give me an initialization interface,
  // and an execution interface.

  // This is not exactly what we expected finally, job_0 should run job_1.
  JobWorkerInterface *job_0 = job_worker_factory.New(0);
  JobWorkerInterface *job_1 = job_worker_factory.New(1);
  job_0->Run();
  job_1->Run();
  return 0;
}
