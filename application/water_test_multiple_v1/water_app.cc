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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "advection.h"
#include "app_config.h"
#include "assert.h"
#include "data_face_arrays.h"
#include "proto_files/adv_vel_par.pb.h"
#include "shared/nimbus.h"
#include "stdio.h"
#include "stdlib.h"
#include "water_app.h"
#include "water_data_driver.h"
#include "water_driver.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

WaterApp::WaterApp() {
};

void WaterApp::Load() {

    printf("Worker beginning to load application\n");

    LOG::Initialize_Logging(false, false, 1<<30, true, 1);

    /* Declare and initialize data, jobs and policies. */

    RegisterJob("main", new Main(this));

    RegisterJob("init", new Init(this));
    RegisterJob("loop", new Loop(this));
    RegisterJob("uptoadvect", new UptoAdvect(this));
    RegisterJob("advect", new Advect(this));
    RegisterJob("afteradvect", new AfterAdvect(this));
    RegisterJob("writeframe", new WriteFrame(this));

    RegisterData("water_driver_1", new WaterDriver<TV>( STREAM_TYPE(T()) ) );
    RegisterData("face_velocities_1", new ::water_app_data
            ::FaceArray<TV>(kMainSize));
    RegisterData("face_velocities_ghost_1", new
            FaceArrayGhost<TV>(kGhostSize));
    RegisterData("sim_data_1", new NonAdvData<TV, T>(kMainSize));

    printf("Finished creating job and data definitions\n");
    printf("Finished loading application\n");

    /* Application data initialization. */

    set_advection_scalar(new ::PhysBAM::ADVECTION_SEMI_LAGRANGIAN_UNIFORM
            < ::PhysBAM::GRID<TV>,T>());

    set_boundary(new ::PhysBAM::BOUNDARY_UNIFORM< ::PhysBAM::GRID<TV>, T>());
    ::PhysBAM::VECTOR< ::PhysBAM::VECTOR<bool, 2>, TV::dimension>
        domain_boundary;
    for(int i=1;i<=TV::dimension;i++)
    {
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

void Main::Execute(std::string params, const DataArray& da)
{
    printf("Begin main\n");

    std::vector<job_id_t> j;
    std::vector<data_id_t> d;
    IDSet<data_id_t> read, write;
    IDSet<job_id_t> before, after;
    IDSet<partition_t> neighbor_partitions;
    partition_t partition_id = 0;
    std::string par;

    GetNewDataID(&d, 4);

    DefineData("water_driver_1", d[0], partition_id, neighbor_partitions, par);
    DefineData("face_velocities_1", d[1], partition_id, neighbor_partitions, par);
    DefineData("face_velocities_ghost_1", d[2], partition_id,
            neighbor_partitions, par);
    DefineData("sim_data_1", d[3], partition_id, neighbor_partitions, par);

    printf("Defined data\n");

    GetNewJobID(&j, 2);

    before.clear();
    after.clear();
    read.clear();
    write.clear();
    after.insert(j[1]);
    read.insert(d[0]);
    read.insert(d[1]);
    read.insert(d[2]);
    read.insert(d[3]);
    SpawnComputeJob("init", j[0], read, write, before, after, par);
    printf("Spawned init\n");

    before.clear();
    after.clear();
    read.clear();
    write.clear();
    before.insert(j[0]);
    read.insert(d[0]);
    read.insert(d[1]);
    read.insert(d[2]);
    read.insert(d[3]);
    write.insert(d[0]);
    write.insert(d[1]);
    write.insert(d[2]);
    write.insert(d[3]);
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

void Init::Execute(std::string params, const DataArray& da)
{
    printf("Executing init job\n");

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        printf("* Data with debug id %i\n", da[i]->get_debug_info());
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                printf("Error : unknown data!!\n");
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    int frame = 0;

    driver->face_velocities = face_velocities;
    driver->sim_data = sim_data;

    Add_Source(sim_data);

    WaterApp *water_app = (WaterApp *)application();
    sim_data->incompressible->Set_Custom_Boundary(*water_app->boundary());
    sim_data->initialize(driver, face_velocities, frame);

    printf("Successfully completed init job\n");
};

UptoAdvect::UptoAdvect(Application *app) {
  set_application(app);
};

Job* UptoAdvect::Clone() {
    printf("Cloning upto advect job\n");
    return new UptoAdvect(application());
};

void UptoAdvect::Execute(std::string params, const DataArray& da)
{

    printf("@@ Running upto advect\n");

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        printf("* Data with debug id %i\n", da[i]->get_debug_info());
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                printf("Error : unknown data!!\n");
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    WaterApp *water_app = (WaterApp *)application();
    sim_data->incompressible->
        Set_Custom_Advection(*(water_app->advection_scalar()));
    sim_data->particle_levelset_evolution->Levelset_Advection(1).
        Set_Custom_Advection(*(water_app->advection_scalar()));
    sim_data->incompressible->Set_Custom_Boundary(*water_app->boundary());

    sim_data->BeforeAdvection(driver, face_velocities);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    printf("Barrier %i\n", x);

};

Advect::Advect(Application *app) {
    set_application(app);
};

Job* Advect::Clone() {
    printf("Cloning advect job\n");
    return new Advect(application());
};

void Advect::Execute(std::string params, const DataArray& da)
{
    printf("@@ Running advect\n");

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        printf("* Data with debug id %i\n", da[i]->get_debug_info());
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                printf("Error : unknown data!!\n");
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    printf("### PARAMETERS \"%s\"\n", params.c_str());
//    ::parameters::AdvVelPar adv_vel_par_pb;
//    adv_vel_par_pb.ParseFromString(params);

    WaterApp *water_app = (WaterApp *)application();
    sim_data->incompressible->
        Set_Custom_Advection(*(water_app->advection_scalar()));
    sim_data->particle_levelset_evolution->Levelset_Advection(1).
        Set_Custom_Advection(*(water_app->advection_scalar()));
    sim_data->incompressible->Set_Custom_Boundary(*water_app->boundary());

    T_FACE_ARRAY *face_vel_extended = 
        new T_FACE_ARRAY(*(face_velocities->grid), kGhostSize);
    face_velocities->Initialize_Ghost_Regions(face_vel_extended,
            water_app->boundary(), kGhostSize,
            driver->dt + driver->time, true);

//    Advect_Velocities(face_velocities, face_vel_extended, water_app,
//            adv_vel_par_pb.dt(), adv_vel_par_pb.time());
    Advect_Velocities(face_velocities, face_vel_extended, water_app,
            driver->dt, driver->time, sim_data);

    delete(face_vel_extended);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    printf("Barrier %i\n", x);
};

AfterAdvect::AfterAdvect(Application *app) {
    set_application(app);
};

Job* AfterAdvect::Clone() {
    printf("Cloning after advect job\n");
    return new AfterAdvect(application());
};

void AfterAdvect::Execute(std::string params, const DataArray& da)
{
    printf("@@ Running after advect\n");

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        printf("* Data with debug id %i\n", da[i]->get_debug_info());
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                printf("Error : unknown data!!\n");
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    WaterApp *water_app = (WaterApp *)application();
    sim_data->incompressible->
        Set_Custom_Advection(*(water_app->advection_scalar()));
    sim_data->particle_levelset_evolution->Levelset_Advection(1).
        Set_Custom_Advection(*(water_app->advection_scalar()));
    sim_data->incompressible->Set_Custom_Boundary(*water_app->boundary());

    sim_data->AfterAdvection(driver, face_velocities);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    printf("Barrier %i\n", x);
};

Loop::Loop(Application *app) {
    set_application(app);
};

Job* Loop::Clone() {
    printf("Cloning loop job\n");
    return new Loop(application());
};

void Loop::Execute(std::string params, const DataArray& da)
{
    printf("Executing forloop job\n");

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        printf("* Data with debug id %i\n", da[i]->get_debug_info());
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                printf("Error : unknown data!!\n");
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    printf("Barrier %i\n", x);

    driver->IncreaseTime();

    if (!driver->CheckProceed())
    {
        printf("... Simulation completed ...\n");
    }
    else
    {
        printf("Spawning new simulation jobs ...\n");

        std::vector<data_id_t> d;
        d.push_back(16777217);
        d.push_back(16777218);
        d.push_back(16777219);
        d.push_back(16777220);
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
        adv_vel_par_pb.set_dt(driver->dt);
        adv_vel_par_pb.set_time(driver->time);
        adv_vel_par_pb.SerializeToString(&par);
        printf("*** PARAMETERS \"%s\"\n", par.c_str());
        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[0]);
        after.insert(j[2]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
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

void WriteFrame::Execute(std::string params, const DataArray& da)
{
    printf( "@@ Executing write frame job\n");

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        printf("* Data with debug id %i\n", da[i]->get_debug_info());
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                printf("Error : unknown data!!\n");
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    if (driver->IsFrameDone()) {
        driver->Write_Output_Files(driver->current_frame);
    }

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    printf("Barrier %i\n", x);
}
