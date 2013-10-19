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
 * This file enumerates the type of jobs, based on location --
 * (corner/ side/ interior).
 * It also defines some functions useful for jobs in the application code, to
 * handle nimbus data objects.
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_JOB_UTILS_H_
#define NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_JOB_UTILS_H_

#include "app_config.h"
#include "data_fwd_decl.h"
#include "physbam_include.h"
#include "shared/nimbus.h"

namespace water_app_job {

    // job locations in a simulation :
    // TODO LATER - need to add jobs on other possible regions
    enum JobRegion {
        kJobAll,
        kJobUpperLeft,
        kJobUpperRight,
        kJobBottomRight,
        kJobBottomLeft,
        kJobNum
    };

    typedef ::water_app_data::FaceArray<TV> FaceArray;
    typedef std::vector<FaceArray * > FaceArrayList;

    struct JobData {
        WaterDriver<TV> *driver;
        NonAdvData<TV, T> *sim_data;
        // this list is interpreted, based on the type of the job
        FaceArrayList boundary_vels;
        // this list is interpreted, based on the type of the job
        FaceArrayList central_vels;
        JobData();
    };

    class SimJob : public Job {
        private:
            JobRegion region_;
        public:
            SimJob() : region_(kJobAll) {}
            JobRegion region() {
                return region_;
            }
            void set_region(JobRegion region) {
                region_ = region;
            }
            void CollectData(const ::nimbus::DataArray& da, JobData& job_data);
    };

} // namespace water_app_job

#endif // NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_JOB_UTILS_H_
