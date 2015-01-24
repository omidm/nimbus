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

#include <unistd.h>
#include <fstream>  // NOLINT
#include <sstream>
#include <string>

#include "worker/worker.h"
#include "worker/worker_manager.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_monitor.h"

namespace nimbus {

WorkerThreadMonitor::WorkerThreadMonitor(WorkerManager* worker_manager)
    : WorkerThread(worker_manager) {
}

WorkerThreadMonitor::~WorkerThreadMonitor() {
}

void WorkerThreadMonitor::Run() {
  std::ofstream output("worker_state.log");
  int64_t dispatched_computation_job_count_last = 0;

  int64_t dispatched_computation_job_count;
  int64_t ready_job_queue_length;

  output << "dispatched_computation_job_count "
      << "working_computation_thread_num "
      << "ready_job_queue_length "
      << std::endl;
  // int count = 0;
  while (true) {
    // count = (count + 1) % 10000;
    // if (count == 0) {
    //   AppDataManager* app_data_manager =
    //       worker_manager_->worker_->application_->app_data_manager();
    //   if (app_data_manager) {
    //     std::stringstream s;
    //     app_data_manager->PrintProfile(&s);
    //     // TODO(quhang): temporary usage.
    //     // printf("\nApplication data manager profile\n%s\n", s.str().c_str());
    //   }
    // }
    usleep(10000);

    dispatched_computation_job_count =
        worker_manager_->dispatched_computation_job_count_;

    ready_job_queue_length = worker_manager_->ready_jobs_count_;

    output << dispatched_computation_job_count
              - dispatched_computation_job_count_last
           << " "
           << ready_job_queue_length
           << std::endl;
    output.flush();

    dispatched_computation_job_count_last =
        dispatched_computation_job_count;
  }
}

}  // namespace nimbus
