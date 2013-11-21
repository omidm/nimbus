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
  * Class representing a physical instance of a data object, stored at
  * a particular worker with a certain version number.
  */

#ifndef NIMBUS_SCHEDULER_PHYSICAL_DATA_H_
#define NIMBUS_SCHEDULER_PHYSICAL_DATA_H_

#include <vector>
#include "shared/nimbus_types.h"

namespace nimbus {

  class PhysicalData {
  public:
    PhysicalData(physical_data_id_t id, worker_id_t worker);
    PhysicalData(physical_data_id_t id, worker_id_t worker, data_version_t version);
    virtual ~PhysicalData();

    physical_data_id_t id();
    worker_id_t worker();
    data_version_t version();
    job_id_t last_job_read();
    job_id_t last_job_write();

    void set_version(data_version_t v);
    void set_last_job_read(job_id_t id);
    void set_last_job_write(job_id_t id);

  private:
    physical_data_id_t id_;
    worker_id_t worker_;
    data_version_t version_;
    job_id_t last_job_read_;
    job_id_t last_job_write_;
  };

  typedef std::vector<PhysicalData> PhysicalDataVector;
}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_PHYSICAL_DATA_H_
