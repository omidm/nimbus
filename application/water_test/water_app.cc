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

#include "lib/nimbus.h"
#include "./water_app.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

WaterApp::WaterApp()
{};

void WaterApp::load() {

    /* Declare and initialize data and jobs. */
    std::cout << "Worker beginning to load application" << std::endl;

    registerJob("main", new Main(this, JOB_COMP));

    registerJob("init", new Init(this, JOB_COMP));
    registerJob("loop", new Loop(this, JOB_COMP));

    registerData("mac_grid_1", new Grid<TV>);
    registerData("mpi_grid_1", new MPIGrid<TV>);
    registerData("face_velocities_1", new FaceArray<TV>);
    registerData("face_velocities_ghost_1", new FaceArrayGhost<TV>);
    registerData("sim_data_1", new NonAdvData<TV, T>);

    std::cout << "Finished creating job and data definitions" << std::endl;
    std::cout << "Finished loading application" << std::endl;
}

Main::Main(Application *app, JobType type)
    : Job(app, type) {};

Job *Main::clone() {
    std::cout << "Cloning main job\n";
    return new Main(application_, type_);
};

void Main::execute(std::string params, const DataArray& da)
{
    std::cout << "Begin main\n" << std::endl;

    std::vector<int> j;
    std::vector<int> d;
    IDSet before, after, read, write;
    std::string par;

    application_->getNewDataID(5, &d);
    application_->defineData("mac_grid_1", d[0]);
    application_->defineData("mpi_grid_1", d[1]);
    application_->defineData("face_velocities_1", d[2]);
    application_->defineData("face_velocities_ghost_1", d[3]);
    application_->defineData("sim_data_1", d[4]);

    application_->getNewJobID(2, &d);

    before.clear();
    after.clear();
    read.clear();
    write.clear();
    par = "";
    after.insert(j[1]);
    read.insert(d[0]);
    read.insert(d[1]);
    read.insert(d[2]);
    read.insert(d[3]);
    read.insert(d[4]);
    application_->spawnJob("init", j[0], before, after, read, write, par);

    before.clear();
    after.clear();
    read.clear();
    write.clear();
    par = "";
    application_->spawnJob("loop", j[1], before, after, read, write, par);
};

Init::Init(Application *app, JobType type)
    : Job(app, type) {};

Job *Init::clone() {
    std::cout << "Cloning init job";
    return new Init(application_, type_);
};

void Init::execute(std::string params, const DataArray& da)
{
};

Loop::Loop(Application *app, JobType type)
    : Job(app, type) {};

Job *Loop::clone() {
    std::cout << "Cloning loop job";
    return new Loop(application_, type_);
};

void Loop::execute(std::string params, const DataArray& da)
{
};
