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

#include <algorithm>
#include "worker/physical_data_map.h"

namespace nimbus {

PhysicalDataMap::PhysicalDataMap() {
  pthread_mutex_init(&lock_, NULL);
}

PhysicalDataMap::~PhysicalDataMap() {
  pthread_mutex_destroy(&lock_);
}

Data* PhysicalDataMap::AcquireAccess(
    physical_data_id_t physical_data_id,
    job_id_t job_id,
    AccessPattern access_pattern) {
  pthread_mutex_lock(&lock_);
  assert(internal_map_.find(physical_data_id) != internal_map_.end());
  AccessState& access_state = internal_map_[physical_data_id];
  Data* result = NULL;
  switch (access_pattern) {
    case READ:
      assert(access_state.initialized);
      if (access_state.flag_write) {
        assert(!access_state.flag_read_and_write);
        assert(access_state.write_job == job_id);
        access_state.flag_read_and_write = true;
      } else {
        access_state.read_jobs.push_back(job_id);
      }
      result = access_state.data;
      break;
    case WRITE:
      assert(access_state.initialized);
      assert(!access_state.flag_write);
      if (access_state.read_jobs.empty()) {
        access_state.flag_write = true;
        access_state.write_job = job_id;
      } else {
        assert(*access_state.read_jobs.begin() == job_id);
        assert((++access_state.read_jobs.begin()) ==
               access_state.read_jobs.end());
        access_state.flag_write = true;
        access_state.flag_read_and_write = true;
        access_state.read_jobs.clear();
        access_state.write_job = job_id;
      }
      result = access_state.data;
      break;
    case INIT:
      assert(!access_state.initialized);
      assert(!access_state.flag_write);
      assert(access_state.read_jobs.empty());
      access_state.initialized = true;
      access_state.flag_write = true;
      access_state.write_job = job_id;
      result = access_state.data;
      break;
    default:
      assert(false);
  }
  pthread_mutex_unlock(&lock_);
  return result;
}

bool PhysicalDataMap::ReleaseAccess(
    physical_data_id_t physical_data_id,
    job_id_t job_id,
    AccessPattern access_pattern) {
  pthread_mutex_lock(&lock_);
  assert(internal_map_.find(physical_data_id) != internal_map_.end());
  AccessState& access_state = internal_map_[physical_data_id];
  switch (access_pattern) {
    case READ: {
      assert(access_state.initialized);
      if (access_state.flag_write) {
        assert(access_state.flag_read_and_write);
        assert(access_state.write_job == job_id);
        access_state.flag_read_and_write = false;
      } else {
        std::list<job_id_t>::iterator iterator =
            std::find(access_state.read_jobs.begin(),
                      access_state.read_jobs.end(),
                      job_id);
        assert(iterator != access_state.read_jobs.end());
        access_state.read_jobs.erase(iterator);
      }
      break;
    }
    case WRITE:
      assert(access_state.initialized);
      assert(access_state.flag_write);
      assert(access_state.write_job == job_id);
      assert(access_state.read_jobs.empty());
      if (access_state.flag_read_and_write) {
        access_state.flag_read_and_write = false;
        access_state.flag_write = false;
        access_state.read_jobs.push_back(job_id);
      } else {
        access_state.flag_write = false;
      }
      break;
    case INIT:
      assert(false);
      assert(access_state.initialized);
      assert(access_state.flag_write);
      assert(access_state.write_job == job_id);
      assert(access_state.read_jobs.empty());
      access_state.flag_write = false;
      break;
    default:
      assert(false);
  }
  pthread_mutex_unlock(&lock_);
  return true;
}

bool PhysicalDataMap::AddMapping(
    physical_data_id_t physical_data_id,
    Data* data) {
  pthread_mutex_lock(&lock_);
  assert(data != NULL);
  assert(internal_map_.find(physical_data_id) == internal_map_.end());
  internal_map_[physical_data_id].data = data;
  pthread_mutex_unlock(&lock_);
  return true;
}

bool PhysicalDataMap::RemoveMapping(
    physical_data_id_t physical_data_id) {
  pthread_mutex_lock(&lock_);
  InternalMap::iterator index = internal_map_.find(physical_data_id);
  assert(index != internal_map_.end());
  internal_map_.erase(index);
  pthread_mutex_unlock(&lock_);
  return true;
}

}  // namespace nimbus
