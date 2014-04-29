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

#include <string>

#include "worker/worker.h"
#include "worker/worker_manager.h"
#include "worker/worker_thread.h"
#include "worker/worker_thread_finish.h"
#include "worker/util_dumping.h"

namespace nimbus {

WorkerThreadFinish::WorkerThreadFinish(
    WorkerManager* worker_manager, PhysicalDataMap* data_map)
    : WorkerThread(worker_manager) {
  data_map_ = data_map;
  assert(data_map_ != NULL);
}

WorkerThreadFinish::~WorkerThreadFinish() {
}

void WorkerThreadFinish::Run() {
  std::list<Job*> finish_job_list;
  while (true) {
    bool success_flag = worker_manager_->PullFinishJobs(this, &finish_job_list);
    assert(success_flag);
    assert(worker_manager_ != NULL);
    for (std::list<Job*>::iterator index = finish_job_list.begin();
         index != finish_job_list.end();
         index++) {
      ProcessJob(*index);
      delete *index;
    }
    finish_job_list.clear();
  }
}

void WorkerThreadFinish::ProcessJob(Job* job) {
  IDSet<physical_data_id_t>::IDSetIter iter;

  IDSet<physical_data_id_t> read = job->read_set();
  for (iter = read.begin(); iter != read.end(); iter++) {
    data_map_->ReleaseAccess(*iter, job->id().elem(), PhysicalDataMap::READ);
  }

  IDSet<physical_data_id_t> write = job->write_set();
  for (iter = write.begin(); iter != write.end(); iter++) {
    data_map_->ReleaseAccess(*iter, job->id().elem(), PhysicalDataMap::WRITE);
  }
  /*
  DataArray daw;
  IDSet<physical_data_id_t> write = job->write_set();
  IDSet<physical_data_id_t>::IDSetIter iter;
  for (iter = write.begin(); iter != write.end(); iter++) {
    daw.push_back((*data_map_)[*iter]);
  }
  // TODO(quhang) no sychronization effort here.
  DumpDataHashInformation(job, daw, data_hash_log_, "hash_out");

  if ((dynamic_cast<CreateDataJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<LocalCopyJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<RemoteCopySendJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<RemoteCopyReceiveJob*>(job) == NULL)) { // NOLINT
    for (iter = write.begin(); iter != write.end(); iter++) {
      // TODO(quhang) assume no sychronization effort for now.
      Data *d = (*data_map_)[*iter];
      data_version_t version = d->version();
      ++version;
      d->set_version(version);
    }
  }
  */

  Parameter params;
  JobDoneCommand cm(job->id(), job->after_set(), params, job->run_time(), job->wait_time());
  worker_manager_->SendCommand(&cm);
}

}  // namespace nimbus
