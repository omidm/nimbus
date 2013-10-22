/* Copyright 2013 Stanford University.
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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "advection.h"
#include "app_config.h"
#include "app_utils.h"
#include "assert.h"
#include "data_face_arrays.h"
#include "data_utils.h"
#include "physbam_include.h"
#include "proto_files/adv_vel_par.pb.h"
#include "shared/nimbus.h"
#include "stdio.h"
#include "stdlib.h"
#include <vector>
#include "water_app.h"
#include "water_data_driver.h"
#include "water_driver.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

/* Assumption: types_list contains a maximum of one kWaterDriver and one
 * kFaceArray. */
#define GetJobData()                                                        \
    WaterApp *water_app = (WaterApp *)application();                        \
    ::application::JobData job_data;                                        \
    CollectData(da, job_data);                                              \
    if (!water_app)                                                         \
        asm volatile("" : : : "memory");
    // asm is just a barrier to allow code to compile with unused variables    

namespace application {
    std::string kDataRegionNames[kDataNum];
    TV_INT kDataRegionSizes[kDataNum];
    std::string kJobRegionNames[kJobNum];
} // namespace application

WaterApp::WaterApp() {
};

void WaterApp::Load() {
    printf("Worker beginning to load application\n");

    /* Initialize application constants. */
    TV_INT main_size(kMainSize, kMainSize);
    TV_INT main_all_size(kMainAllSize, kMainAllSize);
    ::application::GetDataRegionNames(::application::kDataRegionNames);
    ::application::GetDataRegionSizes(
            ::application::kDataRegionSizes,
            main_size,
            main_all_size,
            kGhostSize);
    ::application::GetJobRegionNames(::application::kJobRegionNames);

    LOG::Initialize_Logging(false, false, 1<<30, true, 1);

    /* Declare data types. */
    RegisterData("water_driver", new WaterDriver<TV>( STREAM_TYPE(T()) ) );
    RegisterData("sim_data", new NonAdvData<TV, T>(kMainAllSize));
    /* Declare velocity types. */
    for (int i = 0; i < ::application::kDataNum; i++) {
        std::string ntype_name = "velocities_"+
            ::application::kDataRegionNames[i];
        RegisterData(
                ntype_name,
                new FaceArray(
                    ::application::kDataRegionSizes[i],
                    (::application::DataRegion)i) );
    }

    /* Declare job types. */
    RegisterJob("main", new Main(this));
    RegisterJob("init", new Init(this));
    RegisterJob("loop", new Loop(this));
    RegisterJob("uptoadvect", new UptoAdvect(this));
    for (int i = 0; i < ::application::kJobNum; i++) {
        std::string ntype_name = "advection_"+
            ::application::kJobRegionNames[i];
        RegisterJob(
                ntype_name,
                new Advect(this, (::application::JobRegion)i));
    }
    RegisterJob("afteradvect", new AfterAdvect(this));
    RegisterJob("writeframe", new WriteFrame(this));

    printf("Finished creating job and data definitions\n");
    printf("Finished loading application\n");

    /* Application data initialization -- initialization of partition specific
     * data happens later. */
    set_advection_scalar(new ::PhysBAM::ADVECTION_SEMI_LAGRANGIAN_UNIFORM
            < ::PhysBAM::GRID<TV>,T>());

    set_boundary(new ::PhysBAM::BOUNDARY_UNIFORM< ::PhysBAM::GRID<TV>, T>());
    ::PhysBAM::VECTOR< ::PhysBAM::VECTOR<bool, 2>, TV::dimension>
        domain_boundary;
    for(int i=1;i<=TV::dimension;i++) {
        domain_boundary(i)(1)=true;
        domain_boundary(i)(2)=true;
    }
    domain_boundary(2)(2)=false;
    ::PhysBAM::VECTOR< ::PhysBAM::VECTOR<bool,2>,TV::dimension>
        domain_open_boundaries = ::PhysBAM::VECTOR_UTILITIES
        ::Complement(domain_boundary);
    boundary()->Set_Constant_Extrapolation(domain_open_boundaries);
}

Main::Main(Application *app) {
    set_application(app);
};

Job* Main::Clone() {
    printf("Cloning main job\n");
    return new Main(application());
};

void Main::Execute(std::string params, const DataArray& da) {
    printf("Begin main\n");

    /* Define data as necessary for jobs here. */
    std::string par_data;
    IDSet<partition_t> neighbor_partitions;
    partition_t partition_id = 0;
    // get number of data regions, and corresponding ids and type names for
    // advection
    int data_num = 2;
    std::vector<data_id_t> adv_data_ids[kAdvJobTypesNum];
    std::vector<std::string> adv_data_types[kAdvJobTypesNum];
    for (int i = 0; i < kAdvJobTypesNum; i++) {
        ::application::GetJobDataTypes(kAdvJobTypes[i], adv_data_types[i]);
        int num = adv_data_types[i].size();
        for (int j = 0; j < num; j++) {
            adv_data_ids[i].push_back((data_id_t)(data_num+j));
        }
        data_num += num;
    }
    std::vector<data_id_t> d;
    GetNewDataID(&d, data_num);
    // water driver: non adv data
    DefineData(
            "water_driver",
            d[0],
            partition_id,
            neighbor_partitions,
            par_data);
    // sim data: non adv data
    DefineData("sim_data", d[1], partition_id, neighbor_partitions, par_data);
    // adv data
    for (int i = 0; i < kAdvJobTypesNum; i++) {
        int num = adv_data_types[i].size();
        for (int j = 0; j < num; j++) {
            DefineData(
                    adv_data_types[i][j],
                    adv_data_ids[i][j],
                    partition_id,
                    neighbor_partitions,
                    par_data
                    );
        }
    }

    /* Spawn required jobs -- we have one init job for nonadvection data, and
     * multiple init jobs for advection data; one loop job.
     * Order:
     * - common init
     * - individual inits
     * - loop
     */
    int job_num = 2 + kAdvJobTypesNum;
    std::vector<job_id_t> j;
    GetNewJobID(&j, job_num);
    IDSet<data_id_t> read, write;
    IDSet<job_id_t> before, after;
    std::string par;
    // common init
    before.clear();
    after.clear();
    read.clear();
    write.clear();
    after.insert(j[1]);
    for (unsigned int i = 0; i < d.size(); i++) {
        read.insert(d[i]);
        write.insert(d[i]);
    }
    SpawnComputeJob("init", j[0], read, write, before, after, par);
    printf("Spawned init\n");
    before.clear();
    after.clear();
    read.clear();
    write.clear();
    before.insert(j[0]);
    for (unsigned int i = 0; i < d.size(); i++) {
        read.insert(d[i]);
        write.insert(d[i]);
    }
    SpawnComputeJob("loop", j[1], read, write, before, after, par);
    printf("Spawned loop\n");
    printf("Completed main\n");
};

Init::Init(Application *app) {
    set_application(app);
}

Job* Init::Clone() {
    printf("Cloning init job\n");
    return new Init(application());
};

void Init::Execute(std::string params, const DataArray& da) {
    printf("Executing init job\n");
    std::vector<int> types_list;
    types_list.push_back(driver_id);
    types_list.push_back(non_adv_id);
    for (int i = 0; i < 9; i++) {
        types_list.push_back(face_array_id);
    }
    GetJobData();
    int frame = 0;
    T_FACE_ARRAY *fv = new T_FACE_ARRAY();
    job_data.driver->face_velocities = job_data.central_vels[0]->data();
    job_data.driver->sim_data = job_data.sim_data;
    Add_Source(job_data.sim_data);
    job_data.sim_data->incompressible->
        Set_Custom_Boundary(*water_app->boundary());
    job_data.sim_data->initialize(
            job_data.driver,
            job_data.central_vels[0]->data(),
            frame);
    job_data.sim_data->incompressible->
        Set_Custom_Advection(*(water_app->advection_scalar()));
    job_data.sim_data->particle_levelset_evolution->Levelset_Advection(1).
        Set_Custom_Advection(*(water_app->advection_scalar()));
    delete(fv);
    printf("Successfully completed init job\n");
};

UptoAdvect::UptoAdvect(Application *app) {
  set_application(app);
};

Job* UptoAdvect::Clone() {
    printf("Cloning upto advect job\n");
    return new UptoAdvect(application());
};

void UptoAdvect::Execute(std::string params, const DataArray& da) {
    printf("@@ Running upto advect\n");
    std::vector<int> types_list;
    types_list.push_back(driver_id);
    types_list.push_back(non_adv_id);
    for (int i = 0; i < 9; i++) {
        types_list.push_back(face_array_id);
    }
    GetJobData();
    job_data.sim_data->BeforeAdvection(
            job_data.driver,
            job_data.central_vels[0]->data());
};

Advect::Advect(Application *app, ::application::JobRegion region) {
    set_application(app);
    set_region(region);
};

Job* Advect::Clone() {
    printf("Cloning advect job\n");
    return new Advect(application(), region());
};

void Advect::Execute(std::string params, const DataArray& da) {
    printf("@@ Running advect\n");
    std::vector<int> types_list;
    types_list.push_back(driver_id);
    for (int i = 0; i < 9; i++) {
        types_list.push_back(face_array_id);
    }
    GetJobData();
    ::parameters::AdvVelPar adv_vel_par_pb;
    adv_vel_par_pb.ParseFromString(params);
    T_FACE_ARRAY *face_vel_extended = 
        new T_FACE_ARRAY(*(job_data.central_vels[0]->grid()), kGhostSize);
    job_data.central_vels[0]->Extend_Array(
            job_data.central_vels[0]->data(),
            face_vel_extended,
            water_app->boundary(),
            kGhostSize,
            job_data.driver->dt + job_data.driver->time,
            true);
    Advect_Velocities(job_data.central_vels[0], face_vel_extended, water_app,
            adv_vel_par_pb.dt(), adv_vel_par_pb.time());
    delete(face_vel_extended);
};

AfterAdvect::AfterAdvect(Application *app) {
    set_application(app);
};

Job* AfterAdvect::Clone() {
    printf("Cloning after advect job\n");
    return new AfterAdvect(application());
};

void AfterAdvect::Execute(std::string params, const DataArray& da) {
    printf("@@ Running after advect\n");
    std::vector<int> types_list;
    types_list.push_back(driver_id);
    types_list.push_back(non_adv_id);
    for (int i = 0; i < 9; i++) {
        types_list.push_back(face_array_id);
    }
    GetJobData();
    job_data.sim_data->AfterAdvection(
            job_data.driver,
            job_data.central_vels[0]->data());
};

Loop::Loop(Application *app) {
    set_application(app);
};

Job* Loop::Clone() {
    printf("Cloning loop job\n");
    return new Loop(application());
};

void Loop::Execute(std::string params, const DataArray& da) {
    printf("Executing forloop job\n");
    std::vector<int> types_list;
    types_list.push_back(driver_id);
    types_list.push_back(non_adv_id);
    for (int i = 0; i < 9; i++) {
        types_list.push_back(face_array_id);
    }
    GetJobData();
    job_data.driver->IncreaseTime();
    if (!job_data.driver->CheckProceed()) {
        printf("... Simulation completed ...\n");
    }
    else {
        printf("Spawning new simulation jobs ...\n");

        std::vector<data_id_t> d;
        d.push_back(da[0]->id());
        d.push_back(da[1]->id());
        d.push_back(da[2]->id());
        d.push_back(da[3]->id());
        std::string par;
        IDSet<job_id_t> before, after;
        IDSet<data_id_t> read, write;
        std::vector<job_id_t> j;
        GetNewJobID(&j, 5);
        par = "";
        before.clear(); after.clear();
        read.clear(); write.clear();
        after.insert(j[1]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        SpawnComputeJob("uptoadvect", j[0], read, write, before, after, par);
        printf("Spawned upto advect\n");
        ::parameters::AdvVelPar adv_vel_par_pb;
        adv_vel_par_pb.set_dt(job_data.driver->dt);
        adv_vel_par_pb.set_time(job_data.driver->time);
        adv_vel_par_pb.SerializeToString(&par);
        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[0]);
        after.insert(j[2]);
        read.insert(d[0]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]);
        write.insert(d[2]); write.insert(d[3]);
        SpawnComputeJob("advect", j[1], read, write, before, after, par);
        printf("Spawned advect\n");
        par = "";
        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[1]);
        after.insert(j[3]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        SpawnComputeJob("afteradvect", j[2], read, write, before, after, par);
        printf("Spawned afteradvect\n");
        par = "";
        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[2]);
        after.insert(j[4]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        SpawnComputeJob("writeframe", j[3], read, write, before, after, par);
        printf("Spawned afteradvect\n");
        par = "";
        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[3]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        SpawnComputeJob("loop", j[4], read, write, before, after, par);
        printf("Spawned loop\n");
    }
};

WriteFrame::WriteFrame(Application *app) {
    set_application(app);
};

Job* WriteFrame::Clone() {
    printf("Cloning write frame job\n");
    return new WriteFrame(application());
};

void WriteFrame::Execute(std::string params, const DataArray& da) {
    printf( "@@ Executing write frame job\n");
    std::vector<int> types_list;
    types_list.push_back(driver_id);
    types_list.push_back(non_adv_id);
    for (int i = 0; i < 9; i++) {
        types_list.push_back(face_array_id);
    }
    GetJobData();
    if (job_data.driver->IsFrameDone()) {
        job_data.driver->Write_Output_Files(job_data.driver->current_frame);
    }
}
