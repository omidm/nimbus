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

#ifndef NIMBUS_SHARED_ID_MAKER_H_
#define NIMBUS_SHARED_ID_MAKER_H_

#include <iostream> // NOLINT
#include <vector>
#include <string>
#include "shared/nimbus_types.h"

namespace nimbus {


class IDMaker {
  public:
    IDMaker();
    ~IDMaker();

    void Initialize(worker_id_t worker_id);
    bool GetNewJobID(std::vector<job_id_t>* result, size_t req_num);
    bool GetNewPhysicalDataID(std::vector<physical_data_id_t>* result, size_t req_num);
    bool GetNewLogicalDataID(std::vector<physical_data_id_t>* result, size_t req_num);

  private:
    bool initialized_;
    worker_id_t worker_id_;
    job_id_t first_job_id_;
    job_id_t last_job_id_;
    physical_data_id_t first_physical_data_id_;
    physical_data_id_t last_physical_data_id_;
    logical_data_id_t first_logical_data_id_;
    logical_data_id_t last_logical_data_id_;
    /*
    static const job_id_t JOB_ID_BATCH  = (1 << (sizeof(job_id_t)*3));
    static const physical_data_id_t PHYSICAL_DATA_ID_BATCH  = (1 << (sizeof(physical_data_id_t)*3));
    static const logical_data_id_t LOGICAL_DATA_ID_BATCH  = (1 << (sizeof(logical_data_id_t)*3));
    */
    static const job_id_t JOB_ID_BATCH  = 100000;
    static const physical_data_id_t PHYSICAL_DATA_ID_BATCH  = 100000;
    static const logical_data_id_t LOGICAL_DATA_ID_BATCH  = 100000;
};



}  // namespace nimbus


#endif  // NIMBUS_SHARED_ID_MAKER_H_

