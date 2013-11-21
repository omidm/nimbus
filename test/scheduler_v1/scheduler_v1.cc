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

void SchedulerV1::SchedulerCoreProcessor() {
  while (true) {
    std::cout << "OMID" << std::endl;
    GeometricRegion r;
    PhysicalData p(0, 1);
    data_manager_->AddPartition(0, r);
    data_manager_->AddLogicalObject(0, "omid", 0);
    LogicalDataObject* ldo =
      const_cast<LogicalDataObject*>(data_manager_->FindLogicalObject(0));
    data_manager_->AddPhysicalInstance(ldo, p);
    data_manager_->AddPhysicalInstance(ldo, p);
    data_manager_->AddPhysicalInstance(ldo, p);
    data_manager_->AddPhysicalInstance(ldo, p);
    data_manager_->AddPhysicalInstance(ldo, p);
    data_manager_->AddPhysicalInstance(ldo, p);
    data_manager_->AddPhysicalInstance(ldo, p);
    data_manager_->AddPhysicalInstance(ldo, p);
    PhysicalDataVector pv;
    data_manager_->AllInstances(ldo, &pv);





    sleep(10);
  }
}

