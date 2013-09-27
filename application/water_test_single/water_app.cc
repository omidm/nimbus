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

#include <iostream>
#include "assert.h"
#include "shared/nimbus.h"
#include "./water_app.h"
#include "./water_data_types.h"
#include "./water_driver.h"
#include "stdlib.h"

static int ml = 200;
static int gl = 0;

using namespace PhysBAM;
using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

typedef float T;
typedef VECTOR<T, 2> TV;
typedef VECTOR<int, TV::dimension> TV_INT;

// TODO(chinmayee): alternative to assert?? pointer arithmetic error.

WaterApp::WaterApp()
{};

void WaterApp::Load() {

    std::cout << "Worker beginning to load application" << std::endl;

    LOG::Initialize_Logging(false, false, 1<<30, true, 1);

    /* Declare and initialize data and jobs. */

    RegisterJob("main", new Main(this));

    RegisterJob("init", new Init(this));
    RegisterJob("loop", new Loop(this));
    RegisterJob("uptoadvect", new UptoAdvect(this));
    RegisterJob("advect", new Advect(this));
    RegisterJob("afteradvect", new AfterAdvect(this));
    RegisterJob("writeframe", new WriteFrame(this));

    RegisterData("water_driver_1", new WaterDriver<TV>( STREAM_TYPE(T()) ) );
    RegisterData("face_velocities_1", new FaceArray<TV>(ml));
    RegisterData("face_velocities_ghost_1", new FaceArrayGhost<TV>(gl));
    RegisterData("sim_data_1", new NonAdvData<TV, T>(ml));

    std::cout << "Finished creating job and data definitions" << std::endl;
    std::cout << "Finished loading application" << std::endl;
}

Main::Main(Application *app) {
    set_application(app);
}

Job* Main::Clone() {
    std::cout << "Cloning main job\n";
    return new Main(application());
};

void Main::Execute(std::string params, const DataArray& da)
{
    std::cout << "Begin main\n" << std::endl;

    std::vector<int> j;
    std::vector<int> d;
    IDSet<data_id_t> read, write;
    IDSet<job_id_t> before, after;
    IDSet<partition_t> neighbor_partitions;
    partition_t partition_id = 0;
    std::string par;

    GetNewDataID(4, &d);

    DefineData("water_driver_1", d[0], partition_id, neighbor_partitions, par);
    DefineData("face_velocities_1", d[1], partition_id, neighbor_partitions, par);
    DefineData("face_velocities_ghost_1", d[2], partition_id, neighbor_partitions, par);
    DefineData("sim_data_1", d[3], partition_id, neighbor_partitions, par);

    std::cout << "Defined data\n";

    GetNewJobID(2, &j);

    before.clear();
    after.clear();
    read.clear();
    write.clear();
    after.insert(j[1]);
    read.insert(d[0]);
    read.insert(d[1]);
    read.insert(d[2]);
    read.insert(d[3]);
    std::cout << "Spawning init\n";
    SpawnComputeJob("init", j[0], read, write, before, after, par);
    std::cout << "Spawned init\n";

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
    std::cout << "Spawning loop\n";
    SpawnComputeJob("loop", j[1], read, write, before, after, par);
    std::cout << "Spawned loop\n";

    std::cout << "Completed main\n";
};

Init::Init(Application *app) {
    set_application(app);
}

Job* Init::Clone() {
    std::cout << "Cloning init job\n";
    return new Init(application());
};

void Init::Execute(std::string params, const DataArray& da)
{
    std::cout << "Executing init job\n";

    WaterDriver<TV> *driver = NULL;
    FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        std::cout<<"* Data with debug id "<<da[i]->get_debug_info()<<"\n";
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                std::cout << "Error : unknown data!!\n";
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
    sim_data->initialize(driver, face_velocities, frame);

    std::cout << "Successfully completed init job\n";
};

UptoAdvect::UptoAdvect(Application *app) {
    set_application(app);
};

Job* UptoAdvect::Clone() {
    std::cout << "Cloning upto advect job\n";
    return new UptoAdvect(application());
};

void UptoAdvect::Execute(std::string params, const DataArray& da)
{

    std::cout << "@@ Running upto advect\n";

    WaterDriver<TV> *driver = NULL;
    FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        std::cout<<"* Data with debug id "<<da[i]->get_debug_info()<<"\n";
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                std::cout << "Error : unknown data!!\n";
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    sim_data->BeforeAdvection(driver, face_velocities);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    std::cout << "Barrier "<<x<<"\n";

};

Advect::Advect(Application *app) {
    set_application(app);
};

Job* Advect::Clone() {
    std::cout << "Cloning advect job\n";
    return new Advect(application());
};

void Advect::Execute(std::string params, const DataArray& da)
{
    std::cout << "@@ Running advect\n";

    WaterDriver<TV> *driver = NULL;
    FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        std::cout<<"* Data with debug id "<<da[i]->get_debug_info()<<"\n";
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                std::cout << "Error : unknown data!!\n";
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    face_velocities->Advection(driver, sim_data);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    std::cout << "Barrier "<<x<<"\n";
};

AfterAdvect::AfterAdvect(Application *app) {
    set_application(app);
};

Job* AfterAdvect::Clone() {
    std::cout << "Cloning after advect job\n";
    return new AfterAdvect(application());
};

void AfterAdvect::Execute(std::string params, const DataArray& da)
{
    std::cout << "@@ Running after advect\n";

    WaterDriver<TV> *driver = NULL;
    FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        std::cout<<"* Data with debug id "<<da[i]->get_debug_info()<<"\n";
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                std::cout << "Error : unknown data!!\n";
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    sim_data->AfterAdvection(driver, face_velocities);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    std::cout << "Barrier "<<x<<"\n";
};

Loop::Loop(Application *app) {
    set_application(app);
};

Job* Loop::Clone() {
    std::cout << "Cloning loop job\n";
    return new Loop(application());
};

void Loop::Execute(std::string params, const DataArray& da)
{
    std::cout << "Executing forloop job\n";

    WaterDriver<TV> *driver = NULL;
    FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        std::cout<<"* Data with debug id "<<da[i]->get_debug_info()<<"\n";
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                std::cout << "Error : unknown data!!\n";
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    std::cout << "Barrier "<<x<<"\n";

    driver->IncreaseTime();

    if (!driver->CheckProceed())
    {
        std::cout << "... Simulation completed ...\n";
    }
    else
    {
        std::cout << "Spawning new simulation jobs ...\n";

        std::vector<int> d;
        d.push_back(1);
        d.push_back(2);
        d.push_back(3);
        d.push_back(4);
        std::string par;

        IDSet<job_id_t> before, after;
        IDSet<data_id_t> read, write;

        std::vector<int> j;
        GetNewJobID(5, &j);

        before.clear(); after.clear();
        read.clear(); write.clear();
        after.insert(j[1]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        std::cout << "Spawning upto advect\n";
        SpawnComputeJob("uptoadvect", j[0], read, write, before, after, par);
        std::cout << "Spawned upto advect\n";

        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[0]);
        after.insert(j[2]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        std::cout << "Spawning advect\n";
        SpawnComputeJob("advect", j[1], read, write, before, after, par);
        std::cout << "Spawned advect\n";

        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[1]);
        after.insert(j[3]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        std::cout << "Spawning afteradvect\n";
        SpawnComputeJob("afteradvect", j[2], read, write, before, after, par);
        std::cout << "Spawned afteradvect\n";

        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[2]);
        after.insert(j[4]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        std::cout << "Spawning writeframe\n";
        SpawnComputeJob("writeframe", j[3], read, write, before, after, par);
        std::cout << "Spawned afteradvect\n";

        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[3]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        std::cout << "Spawning loop\n";
        SpawnComputeJob("loop", j[4], read, write, before, after, par);
        std::cout << "Spawned loop\n";
    }
};

WriteFrame::WriteFrame(Application *app) {
    set_application(app);
};

Job* WriteFrame::Clone() {
    std::cout << "Cloning write frame job\n";
    return new WriteFrame(application());
};

void WriteFrame::Execute(std::string params, const DataArray& da)
{
    std::cout << "@@ Executing write frame job\n";

    WaterDriver<TV> *driver = NULL;
    FaceArray<TV> *face_velocities = NULL;
    NonAdvData<TV, T> *sim_data = NULL;

    for (int i = 0; i < 4; i++)
    {
        std::cout<<"* Data with debug id "<<da[i]->get_debug_info()<<"\n";
        switch (da[i]->get_debug_info())
        {
            case driver_id:
                driver = (WaterDriver<TV> *)da[i];
                break;
            case face_array_id:
                face_velocities = (FaceArray<TV> *)da[i];
                break;
            case face_array_ghost_id:
                break;
            case non_adv_id:
                sim_data = (NonAdvData<TV, T> *)da[i];
                break;
            default:
                std::cout << "Error : unknown data!!\n";
                exit(EXIT_FAILURE);
        }
    }

    assert(driver);
    assert(face_velocities);
    assert(sim_data);

    if (driver->IsFrameDone())
    {
        driver->Write_Output_Files(driver->current_frame);
    }

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    std::cout << "Barrier "<<x<<"\n";
}
