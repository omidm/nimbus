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
#include <cstdio>
#include <ctime>
#include "worker/data.h"
#include "worker/physical_data_map.h"

namespace nimbus {

PhysicalDataMap::PhysicalDataMap() {
  sum_ = 0;
  physical_data_log = fopen("physical_data.txt", "w");
}

void PhysicalDataMap::PrintTimeStamp(const char* format, ...) {
  va_list args;
  va_start(args, format);
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(physical_data_log, "%f ", time_sum);
  vfprintf(physical_data_log, format, args);
  fflush(physical_data_log);
  va_end(args);
}

Data* PhysicalDataMap::AcquireAccess(
    physical_data_id_t physical_data_id,
    job_id_t job_id,
    AccessPattern access_pattern) {
  outstanding_used_data_[job_id].insert(physical_data_id);
  assert(internal_map_.find(physical_data_id) != internal_map_.end());
  return internal_map_[physical_data_id].first;
}

bool PhysicalDataMap::ReleaseAccess(
    job_id_t job_id) {
  PhysicalDataIdSet& data_set= outstanding_used_data_[job_id];
  /*
  for (PhysicalDataIdSet::iterator i_physical_data_id = data_set.begin();
       i_physical_data_id != data_set.end();
       ++i_physical_data_id) {
    physical_data_id_t physical_data_id= *i_physical_data_id;
    size_t temp_size = internal_map_[physical_data_id].first->memory_size();
    if (temp_size != internal_map_[physical_data_id].second) {
      sum_ = sum_ + temp_size - internal_map_[physical_data_id].second;
      internal_map_[physical_data_id].second = temp_size;
      PrintTimeStamp("%s %"PRIu64" %zu\n",
                     internal_map_[physical_data_id].first->name().c_str(),
                     physical_data_id,
                     temp_size);
      PrintTimeStamp("%zu\n", sum_);
    }
  }
  */
  outstanding_used_data_.erase(job_id);
  return true;
}

bool PhysicalDataMap::AddMapping(
    physical_data_id_t physical_data_id,
    Data* data) {
  assert(data != NULL);
  assert(internal_map_.find(physical_data_id) == internal_map_.end());
  internal_map_[physical_data_id].first = data;
  /*
  internal_map_[physical_data_id].second = data->memory_size();
  PrintTimeStamp("%s %"PRIu64" %zu\n",
                 data->name().c_str(),
                 physical_data_id,
                 data->memory_size());
                 */
  return true;
}

bool PhysicalDataMap::RemoveMapping(
    physical_data_id_t physical_data_id) {
  // TODO(quhang): this function should not be used.
  assert(false);
  InternalMap::iterator index = internal_map_.find(physical_data_id);
  assert(index != internal_map_.end());
  internal_map_.erase(index);
  return true;
}

}  // namespace nimbus
