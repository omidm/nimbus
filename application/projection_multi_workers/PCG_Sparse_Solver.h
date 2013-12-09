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
 * Author: Zhihao Jia <zhihao@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_PROJECTION_TEST_PCG_SPARSE_SOLVER_H_
#define NIMBUS_APPLICATION_PROJECTION_TEST_PCG_SPARSE_SOLVER_H_
#include "shared/nimbus.h"
#include "shared/parser.h"
#include "shared/nimbus_types.h"
#include "worker/application.h"
#include "worker/job.h"
#include "worker/data.h"
#include "protocol_buffer/Sparse_Matrix_Float.pb.h"
#include "protocol_buffer/Vector_Float.pb.h"

#include "nimbus_pcg_sparse_mpi.h"
#define DESIRED_INTERATIONS 200

using nimbus::Job;
using nimbus::Data;
using nimbus::Application;
using namespace PhysBAM;
typedef float T;

class Init : public Job {
    public:
    	Init(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Condition : public Job {
    public:
    	Project_Forloop_Condition(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part1 : public Job {
    public:
    	Project_Forloop_Part1(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part2 : public Job {
    public:
    	Project_Forloop_Part2(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part3 : public Job {
    public:
    	Project_Forloop_Part3(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part4 : public Job {
    public:
    	Project_Forloop_Part4(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Global_Sum : public Job {
    public:
		Global_Sum(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Global_Max_Abs : public Job {
    public:
    	Global_Max_Abs(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Finish : public Job {
    public:
    	Finish(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};
