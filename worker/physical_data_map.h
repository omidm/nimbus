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

#ifndef NIMBUS_WORKER_PHYSICAL_DATA_MAP_H_
#define NIMBUS_WORKER_PHYSICAL_DATA_MAP_H_

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <pthread.h>
#include <cassert>
#include <list>
#include <utility>

#include "shared/nimbus_types.h"

namespace nimbus {
class Data;
class PhysicalDataMap {
 public:
  enum AccessPattern {
    READ,
    WRITE,
    INIT
  };
  explicit PhysicalDataMap();
  virtual ~PhysicalDataMap() {}

  Data* AcquireAccess(
      physical_data_id_t physical_data_id,
      job_id_t job_id,
      AccessPattern access_pattern);
  bool ReleaseAccess(job_id_t job_id);

  bool AddMapping(
      physical_data_id_t physical_data_id,
      Data* data);
  bool RemoveMapping(
      physical_data_id_t physical_data_id);

 private:
  size_t sum_;
  typedef boost::unordered_set<physical_data_id_t> PhysicalDataIdSet;
  boost::unordered_map<job_id_t, PhysicalDataIdSet> outstanding_used_data_;
  typedef boost::unordered_map<physical_data_id_t, std::pair<Data*, size_t> >
      InternalMap;
  InternalMap internal_map_;
};
}  // namespace nimbus

#endif  // NIMBUS_WORKER_PHYSICAL_DATA_MAP_H_
