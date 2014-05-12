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

#include <pthread.h>
#include <cassert>
#include <list>
#include <map>

#include "shared/nimbus_types.h"

#ifndef MUTE_DATA_ACCESS_CHECK
#define MUTE_DATA_ACCESS_CHECK
#endif

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
  virtual ~PhysicalDataMap();

  Data* AcquireAccess(
      physical_data_id_t physical_data_id,
      job_id_t job_id,
      AccessPattern access_pattern);
  bool ReleaseAccess(
      physical_data_id_t physical_data_id,
      job_id_t job_id,
      AccessPattern access_pattern);

  bool AddMapping(
      physical_data_id_t physical_data_id,
      Data* data);
  bool RemoveMapping(
      physical_data_id_t physical_data_id);

 private:
  pthread_mutex_t lock_;
  // job_id_t
  struct AccessState {
    Data* data;
    bool initialized;
    bool flag_write;
    bool flag_read_and_write;
    job_id_t write_job;
    std::list<job_id_t> read_jobs;
    AccessState() {
      data = NULL;
      initialized = false;
      flag_write = false;
      flag_read_and_write = false;
    }
  };
  typedef std::map<physical_data_id_t, AccessState> InternalMap;
  InternalMap internal_map_;
};
}  // namespace nimbus

#endif  // NIMBUS_WORKER_PHYSICAL_DATA_MAP_H_
