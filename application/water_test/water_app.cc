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
#include "./water_driver.h"

static int ml = 200;
static int gl = 0;

using namespace PhysBAM;
using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

typedef float T;
typedef VECTOR<T, 2> TV;
typedef VECTOR<int, TV::dimension> TV_INT;

WaterApp::WaterApp()
{};

void WaterApp::load() {

    WaterDriver<TV> *driver = new WaterDriver<TV> ( STREAM_TYPE(T()) );
    set_app_data(driver);

    /* Declare and initialize data and jobs. */
    std::cout << "Worker beginning to load application" << std::endl;

    registerJob("main", new Main(this, JOB_COMP));

    registerJob("init", new Init(this, JOB_COMP));
    registerJob("loop", new Loop(this, JOB_COMP));

    registerData("face_velocities_1", new FaceArray<TV>(ml));
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
    IDSet<job_id_t> before, after;
    IDSet<data_id_t> read, write;
    std::string par = "none";

    application_->getNewDataID(3, &d);

    application_->defineData("face_velocities_1", d[0]);
    application_->defineData("face_velocities_ghost_1", d[1]);
    application_->defineData("sim_data_1", d[2]);

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
    std::cout << "Spawning init\n";
    application_->spawnJob("init", j[0], before, after, read, write, par);
    std::cout << "Spawned init\n";

    before.clear();
    after.clear();
    read.clear();
    write.clear();
    std::cout << "Spawning loop\n";
    application_->spawnJob("loop", j[1], before, after, read, write, par);
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
    std::cout << "Executing init jon\n";
    void *data = application_->app_data();
    assert(data!=NULL);
    WaterDriver<TV> *driver = (WaterDriver<TV> *)data;
    // TODO(chinmayee): Implement driver initialize to set constant policies
    // and parameter values
    driver->initialize(false);
    // TODO(chinmayee): Initializa all the remaining data over here
    // initialize mac grid
};

Loop::Loop(Application *app, JobType type)
    : Job(app, type) {};

Job* Loop::clone() {
    std::cout << "Cloning loop job\n";
    return new Loop(application_, type_);
};

void Loop::execute(std::string params, const DataArray& da)
{};
