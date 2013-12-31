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
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./scheduler_v1.h"
#define WORKER_NUM 2

SchedulerV1::SchedulerV1(unsigned int p)
: Scheduler(p) {
}

bool SchedulerV1::GetWorkerToAssignJob(JobEntry* job, SchedulerWorker*& worker) {
  // Assumption is that partition Ids start from 0, and incrementally go up.
  size_t worker_num = server_->worker_num();
  size_t chunk = (data_manager_->max_defined_partition() + 1) / worker_num;
  std::vector<int> workers_rank(worker_num, 0);

  IDSet<logical_data_id_t> union_set = job->union_set();
  IDSet<logical_data_id_t>::IDSetIter iter;
  for (iter = union_set.begin(); iter != union_set.end(); ++iter) {
    const LogicalDataObject* ldo;
    ldo = data_manager_->FindLogicalObject(*iter);
    size_t poll = std::min((size_t)(ldo->partition()) / chunk, worker_num - 1);
    workers_rank[poll] = workers_rank[poll] + 1;
  }

  // find the worker that wins the poll.
  worker_id_t w_id = 1;
  int count = workers_rank[0];
  for (size_t i = 1; i < worker_num; ++i) {
    if (count < workers_rank[i]) {
      count = workers_rank[i];
      w_id = i + 1;
    }
  }

  return server_->GetSchedulerWorkerById(worker, w_id);
}

/*
void SchedulerV1::SchedulerCoreProcessor() {
  while (true) {
    std::cout << "OMID" << std::endl;
    GeometricRegion r;
    PhysicalData p1(0, 1);
    PhysicalData p2(0, 2);
    data_manager_->AddPartition(0, r);
    data_manager_->AddLogicalObject(0, "omid", 0);
    LogicalDataObject* ldo =
      const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(0));
    data_manager_->AddPhysicalInstance(ldo, p1);
    data_manager_->AddPhysicalInstance(ldo, p1);
    data_manager_->AddPhysicalInstance(ldo, p1);
    data_manager_->AddPhysicalInstance(ldo, p1);
    data_manager_->AddPhysicalInstance(ldo, p2);
    data_manager_->AddPhysicalInstance(ldo, p2);
    data_manager_->AddPhysicalInstance(ldo, p2);
    data_manager_->AddPhysicalInstance(ldo, p2);
    PhysicalDataVector pv;
    data_manager_->AllInstances(ldo, &pv);
    std::cout << pv.size() << std::endl;
    data_manager_->InstancesByWorker(ldo, 1, &pv);
    std::cout << pv.size() << std::endl;
    data_manager_->InstancesByWorker(ldo, 2, &pv);
    std::cout << pv.size() << std::endl;
    data_manager_->InstancesByWorker(ldo, 3, &pv);
    std::cout << pv.size() << std::endl;

    sleep(10);
  }
}
*/

