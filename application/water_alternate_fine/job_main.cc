/* Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 vd* are met:
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
 * This file contains the "main" job that Nimbus launches after loading an
 * application. All subsequent jobs are spawned from here.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/data_app.h"
#include "application/water_alternate_fine/job_initialize.h"
#include "application/water_alternate_fine/job_main.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <vector>

namespace application {

    JobMain::JobMain(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobMain::Clone() {
        return new JobMain(application());
    }

    void JobMain::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing main job\n");

        // Partition setup
        nimbus::ID<partition_id_t> partition_id1(0);
        nimbus::ID<partition_id_t> partition_id2(1);
        nimbus::ID<partition_id_t> partition_id3(2);
        nimbus::ID<partition_id_t> partition_id4(3);
        nimbus::Parameter part_params;
        DefinePartition(partition_id1, kDomainFaceVel, part_params);
        DefinePartition(partition_id2, kDomainPhi, part_params);
        DefinePartition(partition_id3, kDomainPressure, part_params);
        DefinePartition(partition_id4, kDomainParticles, part_params);
        nimbus::IDSet<partition_id_t> neighbor_partitions;

        // Data setup
        int data_num = 8;
        std::vector<logical_data_id_t> data_ids;
        GetNewLogicalDataID(&data_ids, data_num);

        // Face arrays
        nimbus::Parameter fa_params;
        fa_params.set_ser_data(SerializedData(""));
        DefineData(APP_FACE_VEL,
                   data_ids[0],
                   partition_id1.elem(),
                   neighbor_partitions,
                   fa_params);

        // Scalar arrays
        nimbus::Parameter sa_params;
        sa_params.set_ser_data(SerializedData(""));
        DefineData(APP_PHI,
                   data_ids[1],
                   partition_id2.elem(),
                   neighbor_partitions,
                   sa_params);
        DefineData(APP_PRESSURE,
                   data_ids[2],
                   partition_id3.elem(),
                   neighbor_partitions,
                   sa_params);

        // Particles
        nimbus::Parameter particle_params;
        part_params.set_ser_data(SerializedData(""));
        DefineData(APP_POS_PARTICLES,
                   data_ids[3],
                   partition_id4.elem(),
                   neighbor_partitions,
                   particle_params);
        DefineData(APP_NEG_PARTICLES,
                   data_ids[4],
                   partition_id4.elem(),
                   neighbor_partitions,
                   particle_params);
        DefineData(APP_POS_REM_PARTICLES,
                   data_ids[5],
                   partition_id4.elem(),
                   neighbor_partitions,
                   particle_params);
        DefineData(APP_NEG_REM_PARTICLES,
                   data_ids[6],
                   partition_id4.elem(),
                   neighbor_partitions,
                   particle_params);
        DefineData(APP_LAST_UNIQUE_PARTICLE_ID,
                   data_ids[7],
                   partition_id4.elem(),
                   neighbor_partitions,
                   particle_params);

        // Job setup
        int job_num = 1;
        std::vector<nimbus::job_id_t> job_ids;
        GetNewJobID(&job_ids, job_num);
        nimbus::IDSet<nimbus::logical_data_id_t> read, write;
        nimbus::IDSet<nimbus::job_id_t> before, after;

        // Init job
        read.insert(data_ids[0]);
        read.insert(data_ids[1]);
        read.insert(data_ids[2]);
        read.insert(data_ids[3]);
        read.insert(data_ids[4]);
        read.insert(data_ids[5]);
        read.insert(data_ids[6]);
        read.insert(data_ids[7]);
        write.insert(data_ids[0]);
        write.insert(data_ids[1]);
        write.insert(data_ids[2]);
        write.insert(data_ids[3]);
        write.insert(data_ids[4]);
        write.insert(data_ids[5]);
        write.insert(data_ids[6]);
        write.insert(data_ids[7]);
        nimbus::Parameter init_params;
        init_params.set_ser_data(SerializedData(""));
        dbg(APP_LOG, "Spawning initialize\n");
        SpawnComputeJob(INITIALIZE,
                        job_ids[0],
                        read, write,
                        before, after,
                        init_params);

        dbg(APP_LOG, "Completed executing main job\n");
    }

} // namespace application
