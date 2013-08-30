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
#include "assert.h"
#include "data_face_arrays.h"
#include <iostream>
#include "shared/nimbus.h"
#include "stdlib.h"
#include "types.h"
#include "water_app.h"
#include "water_data_driver.h"
#include "water_driver.h"

static int ml = 200;
static int gl = 0;

using namespace PhysBAM;
using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

// TODO(chinmayee): alternative to assert?? pointer arithmetic error.

WaterApp::WaterApp()
{};

void WaterApp::load() {

    std::cout << "Worker beginning to load application" << std::endl;

    LOG::Initialize_Logging(false, false, 1<<30, true, 1);

    /* Declare and initialize data and jobs. */

    registerJob("main", new Main(this, JOB_COMP));

    registerJob("init", new Init(this, JOB_COMP));
    registerJob("loop", new Loop(this, JOB_COMP));
    registerJob("uptoadvect", new UptoAdvect(this, JOB_COMP));
    registerJob("advect", new Advect(this, JOB_COMP));
    registerJob("afteradvect", new AfterAdvect(this, JOB_COMP));
    registerJob("writeframe", new WriteFrame(this, JOB_COMP));

    registerData("water_driver_1", new WaterDriver<TV>( STREAM_TYPE(T()) ) );
    registerData("face_velocities_1", new ::water_app_data::FaceArray<TV>(ml));
    registerData("face_velocities_ghost_1", new FaceArrayGhost<TV>(gl));
    registerData("sim_data_1", new NonAdvData<TV, T>(ml));

    std::cout << "Finished creating job and data definitions" << std::endl;
    std::cout << "Finished loading application" << std::endl;
}

Main::Main(Application *app, JobType type)
    : Job(app, type) {};

Job* Main::clone() {
    std::cout << "Cloning main job\n";
    return new Main(application_, type_);
};

void Main::execute(std::string params, const DataArray& da)
{
    std::cout << "Begin main\n" << std::endl;

    std::vector<int> j;
    std::vector<int> d;
    IDSet<data_id_t> read, write;
    IDSet<job_id_t> before, after;
    IDSet<partition_t> neighbor_partitions;
    partition_t partition_id = 0;
    std::string par;

    application_->getNewDataID(4, &d);

    application_->DefineData("water_driver_1", d[0], partition_id, neighbor_partitions, par);
    application_->DefineData("face_velocities_1", d[1], partition_id, neighbor_partitions, par);
    application_->DefineData("face_velocities_ghost_1", d[2], partition_id, neighbor_partitions, par);
    application_->DefineData("sim_data_1", d[3], partition_id, neighbor_partitions, par);

    std::cout << "Defined data\n";

    application_->getNewJobID(2, &j);

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
    application_->SpawnJob("init", j[0], read, write, before, after, JOB_COMP, par);
    std::cout << "Spawned init\n";

    before.clear();
    after.clear();
    read.clear();
    write.clear();
    read.insert(d[0]);
    read.insert(d[1]);
    read.insert(d[2]);
    read.insert(d[3]);
    write.insert(d[0]);
    write.insert(d[1]);
    write.insert(d[2]);
    write.insert(d[3]);
    std::cout << "Spawning loop\n";
    application_->SpawnJob("loop", j[1], read, write, before, after, JOB_COMP, par);
    std::cout << "Spawned loop\n";

    std::cout << "Completed main\n";
};

Init::Init(Application *app, JobType type)
    : Job(app, type) {};

Job* Init::clone() {
    std::cout << "Cloning init job\n";
    return new Init(application_, type_);
};

void Init::execute(std::string params, const DataArray& da)
{
    std::cout << "Executing init job\n";

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
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
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
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

UptoAdvect::UptoAdvect(Application *app, JobType type)
    : Job(app, type) {};

Job* UptoAdvect::clone() {
    std::cout << "Cloning upto advect job\n";
    return new UptoAdvect(application_, type_);
};

void UptoAdvect::execute(std::string params, const DataArray& da)
{

    std::cout << "@@ Running upto advect\n";

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
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
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
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

Advect::Advect(Application *app, JobType type)
    : Job(app, type) {};

Job* Advect::clone() {
    std::cout << "Cloning advect job\n";
    return new Advect(application_, type_);
};

void Advect::execute(std::string params, const DataArray& da)
{
    std::cout << "@@ Running advect\n";

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
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
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
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

    Advection(face_velocities, sim_data);

    int x = driver->get_debug_info() + face_velocities->get_debug_info() +
        sim_data->get_debug_info();
    std::cout << "Barrier "<<x<<"\n";
};

AfterAdvect::AfterAdvect(Application *app, JobType type)
    : Job(app, type) {};

Job* AfterAdvect::clone() {
    std::cout << "Cloning after advect job\n";
    return new AfterAdvect(application_, type_);
};

void AfterAdvect::execute(std::string params, const DataArray& da)
{
    std::cout << "@@ Running after advect\n";

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
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
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
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

Loop::Loop(Application *app, JobType type)
    : Job(app, type) {};

Job* Loop::clone() {
    std::cout << "Cloning loop job\n";
    return new Loop(application_, type_);
};

void Loop::execute(std::string params, const DataArray& da)
{
    std::cout << "Executing forloop job\n";

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
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
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
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
        application_->getNewJobID(5, &j);

        before.clear(); after.clear();
        read.clear(); write.clear();
        after.insert(j[1]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        std::cout << "Spawning upto advect\n";
        application_->SpawnJob("uptoadvect", j[0], read, write, before, after, JOB_COMP, par);
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
        application_->SpawnJob("advect", j[1], read, write, before, after, JOB_COMP, par);
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
        application_->SpawnJob("afteradvect", j[2], read, write, before, after, JOB_COMP, par);
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
        application_->SpawnJob("writeframe", j[3], read, write, before, after, JOB_COMP, par);
        std::cout << "Spawned afteradvect\n";

        before.clear(); after.clear();
        read.clear(); write.clear();
        before.insert(j[3]);
        read.insert(d[0]); read.insert(d[1]);
        read.insert(d[2]); read.insert(d[3]);
        write.insert(d[0]); write.insert(d[1]);
        write.insert(d[2]); write.insert(d[3]);
        std::cout << "Spawning loop\n";
        application_->SpawnJob("loop", j[4], read, write, before, after, JOB_COMP, par);
        std::cout << "Spawned loop\n";
    }
};

WriteFrame::WriteFrame(Application *app, JobType type)
    : Job(app, type) {};

Job* WriteFrame::clone() {
    std::cout << "Cloning write frame job\n";
    return new WriteFrame(application_, type_);
};

void WriteFrame::execute(std::string params, const DataArray& da)
{
    std::cout << "@@ Executing write frame job\n";

    WaterDriver<TV> *driver = NULL;
    ::water_app_data::FaceArray<TV> *face_velocities = NULL;
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
                face_velocities = (::water_app_data::FaceArray<TV> *)da[i];
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
