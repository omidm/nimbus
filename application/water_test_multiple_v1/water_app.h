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

#ifndef NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_WATER_APP_H_
#define NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_WATER_APP_H_

#include "app_config.h"
#include "data_utils.h"
#include "job_utils.h"
#include "physbam_include.h"
#include "shared/nimbus.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

/* Application class launched by Nimbus. Initialization of jobs, using
 * functions in water_driver, should be done here. Methods to initialize
 * simulation data and build the data map should also be called here.
 */
class WaterApp : public Application {

    private:
        ::PhysBAM::ADVECTION_SEMI_LAGRANGIAN_UNIFORM< ::PhysBAM::GRID<TV>, T>
            *advection_scalar_;
        ::PhysBAM::BOUNDARY_UNIFORM< ::PhysBAM::GRID<TV>, T> *boundary_;

    public:
        WaterApp();
        virtual void Load();

        ::PhysBAM::ADVECTION_SEMI_LAGRANGIAN_UNIFORM< ::PhysBAM::GRID<TV>, T>*
            advection_scalar() {
                return advection_scalar_;
            }

        void set_advection_scalar(
                ::PhysBAM::ADVECTION_SEMI_LAGRANGIAN_UNIFORM< ::PhysBAM::
                GRID<TV>, T> *advection_scalar) {
            advection_scalar_ = advection_scalar;
        }

        ::PhysBAM::BOUNDARY_UNIFORM< ::PhysBAM::GRID<TV>, T>* boundary() {
            return boundary_;
        }

        void set_boundary(
                ::PhysBAM::BOUNDARY_UNIFORM< ::PhysBAM::GRID<TV>, T> *boundary) {
            boundary_ = boundary;
        }
};

class Main : public Job {
    public:
        Main(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

class Init : public ::water_app_job::SimJob {
    public:
        Init(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

class Loop : public ::water_app_job::SimJob {
    public:
        Loop(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

class UptoAdvect : public ::water_app_job::SimJob {
    public:
        UptoAdvect(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

class Advect : public ::water_app_job::SimJob {
    public:
        Advect(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

class AfterAdvect : public ::water_app_job::SimJob {
    public:
        AfterAdvect(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

class WriteFrame : public ::water_app_job::SimJob {
    public:
        WriteFrame(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

#endif  // NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_WATER_APP_H_
