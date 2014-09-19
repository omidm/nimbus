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
  * Class for producing unique job and data id.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "shared/id_maker.h"

using namespace nimbus; // NOLINT

IDMaker::IDMaker() {
  initialized_ = false;
  pthread_mutex_init(&lock_, NULL);
}

IDMaker::~IDMaker() {
  pthread_mutex_destroy(&lock_);
}

void IDMaker::Initialize(worker_id_t worker_id) {
  pthread_mutex_lock(&lock_);
  worker_id_ = worker_id;
  initialized_ = true;
  first_job_id_ = JOB_ID_BATCH * worker_id + NIMBUS_KERNEL_JOB_ID;
  last_job_id_ = first_job_id_;
  first_physical_data_id_ = PHYSICAL_DATA_ID_BATCH * worker_id;
  last_physical_data_id_ = first_physical_data_id_;
  first_logical_data_id_ = LOGICAL_DATA_ID_BATCH * worker_id;
  last_logical_data_id_ = first_logical_data_id_;
  pthread_mutex_unlock(&lock_);
}

bool IDMaker::GetNewJobID(std::vector<job_id_t>* result, size_t req_num) {
  pthread_mutex_lock(&lock_);
  result->clear();
  if (!initialized_) {
    std::cout << "ERROR: IDMaker has not been initialized yet." << std::endl;
    pthread_mutex_unlock(&lock_);
    return false;
  }
  if ((last_job_id_ + req_num) >= (first_job_id_ + JOB_ID_BATCH)) {
    std::cout << "ERROR: IDMaker ran out of job IDs." << std::endl;
    pthread_mutex_unlock(&lock_);
    return false;
  }
  for (size_t i = 0; i < req_num; i++)
    result->push_back(++last_job_id_);
  pthread_mutex_unlock(&lock_);
  return true;
}

bool IDMaker::GetNewPhysicalDataID(std::vector<physical_data_id_t>* result, size_t req_num) {
  pthread_mutex_lock(&lock_);
  result->clear();
  if (!initialized_) {
    std::cout << "ERROR: IDMaker has not been initialized yet." << std::endl;
    pthread_mutex_unlock(&lock_);
    return false;
  }
  if ((last_physical_data_id_ + req_num) >= (first_physical_data_id_ + PHYSICAL_DATA_ID_BATCH)) {
    std::cout << "ERROR: IDMaker ran out of physical data IDs." << std::endl;
    pthread_mutex_unlock(&lock_);
    return false;
  }
  for (size_t i = 0; i < req_num; i++)
    result->push_back(++last_physical_data_id_);
  pthread_mutex_unlock(&lock_);
  return true;
}

bool IDMaker::GetNewLogicalDataID(std::vector<logical_data_id_t>* result, size_t req_num) {
  pthread_mutex_lock(&lock_);
  result->clear();
  if (!initialized_) {
    std::cout << "ERROR: IDMaker has not been initialized yet." << std::endl;
    pthread_mutex_unlock(&lock_);
    return false;
  }
  if ((last_logical_data_id_ + req_num) >= (first_logical_data_id_ + LOGICAL_DATA_ID_BATCH)) {
    std::cout << "ERROR: IDMaker ran out of logical data IDs." << std::endl;
    pthread_mutex_unlock(&lock_);
    return false;
  }
  for (size_t i = 0; i < req_num; i++)
    result->push_back(++last_logical_data_id_);
  pthread_mutex_unlock(&lock_);
  return true;
}

bool IDMaker::SchedulerProducedJobID(job_id_t job_id) {
  static job_id_t first = JOB_ID_BATCH * NIMBUS_SCHEDULER_ID + NIMBUS_KERNEL_JOB_ID;
  static job_id_t last  = JOB_ID_BATCH + first;
  return ((job_id > first) && (job_id < last));
}


