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

#include "Projection_App.h"

ProjectionApp::ProjectionApp() {
};

void ProjectionApp::Load() {
	printf("Worker beginning to load application\n");

	/* Declare and initialize data, jobs and policies. */

	RegisterJob("main", new Main(this));
	RegisterJob("Init", new Init(this));
	RegisterJob("Project_Forloop_Condition", new Project_Forloop_Condition(this));
	RegisterJob("Project_Forloop_Part1", new Project_Forloop_Part1(this));
	RegisterJob("Project_Forloop_Part2", new Project_Forloop_Part2(this));
	RegisterJob("Project_Forloop_Part3", new Project_Forloop_Part3(this));
	RegisterJob("Project_Forloop_Part4", new Project_Forloop_Part4(this));
	RegisterJob("Global_Sum", new Global_Sum(this));
	RegisterJob("Global_Max_Abs", new Global_Max_Abs(this));

	RegisterData("vector", new Vec(LEN));
	RegisterData("scalar", new Vec(1));
	RegisterData("matrix", new Vec((LEN*LEN)));

	printf("Finished creating job and data definitions\n");
}

Main::Main(Application *app) {
	set_application(app);
};

Job* Main::Clone() {
	printf("Cloning main job\n");
	return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
	printf("Begin main\n");
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	partition_id_t pid1 = 1, pid2 = 2;
	Parameter par;
	IDSet<param_id_t> param_idset;

	GetNewLogicalDataID(&d, NUM_OF_FORLOOP_INPUTS);

	GetNewJobID(&j, 3);
	
	// Proj_Initialize, pid = 1
	before.clear();
	after.clear();
	after.insert(j[1]);
	read.clear();
	write.clear();
	SpawnComputeJob("Proj_Initialize", j[0], read, write, before, after, par);

    for (int i = 1; i <= example.last_frame; i++) {
        before.clear();
        before.insert(j[i-1]);
        after.clear();
        after.insert(j[i+1]);
        read.clear();
        write.clear();
        SpwanComputeJob("Proj_AdvanceOneTime", j[i], read, write, before, after, par);
    }

	printf("Completed main\n");
};

Proj_Initialize::Proj_Initialize(Application *app) {
	set_application(app);
};

Job* Proj_Initialize::Clone() {
	printf("Cloning Proj_Initialize job\n");
	return new Proj_Initialize(application());
};

void Proj_Initialize::Execute(Parameter params, const DataArray& da) {
	printf("Begin Proj_Initialize\n");

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;
    time=example.Time_At_Frame(current_frame);

    // will map mpi_grid into Nimbus job dependency
    if(example.mpi_grid) example.mpi_grid->Initialize(example.domain_boundary);
    example.projection.elliptic_solver->mpi_grid=example.mpi_grid;
    if(example.mpi_grid) example.boundary=new BOUNDARY_MPI<GRID<TV>,T>(example.mpi_grid,example.boundary_scalar);
    else example.boundary=&example.boundary_scalar;

    // setup grids and velocities
    example.projection.Initialize_Grid(example.mac_grid);
    example.face_velocities.Resize(example.mac_grid);
    example.Initialize_Fields();

    // setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance(1e-11);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(200);
    example.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations=40;

    if(example.restart) example.Read_Output_Files(example.restart);

    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.Set_Boundary_Conditions(time); // get so CFL is correct

	printf("Completed Proj_Initialize\n");
};

Proj_AdvanceOnTime::Proj_AdvanceOnTime(Application *app) {
	set_application(app);
};

Job* Proj_AdvanceOnTime::Clone() {
	printf("Cloning Proj_AdvanceOnTime job\n");
	return new Proj_AdvanceOnTime(application());
};

void Proj_AdvanceOnTime::Execute(Parameter params, const DataArray& da) {
	printf("Begin Proj_AdvanceOnTime\n");


	printf("Completed Proj_AdvanceOnTime\n");
};

Proj_MainProjection::Proj_MainProjection(Application *app) {
	set_application(app);
};

Job* Proj_MainProjection::Clone() {
	printf("Cloning Proj_MainProjection job\n");
	return new Proj_MainProjection(application());
};

void Proj_MainProjection::Execute(Parameter params, const DataArray& da) {
	printf("Begin Proj_MainProjection\n");
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	partition_id_t pid1 = 1, pid2 = 2;
	Parameter par;
	IDSet<param_id_t> param_idset;

	GetNewLogicalDataID(&d, NUM_OF_FORLOOP_INPUTS);
	DefineData("scalar", d[0], pid1, neighbor_partitions, par); // residual
	DefineData("matrix", d[1], pid1, neighbor_partitions, par); // A_pid1
	DefineData("matrix", d[2], pid2, neighbor_partitions, par); // A_pid2
	DefineData("vector", d[3], pid1, neighbor_partitions, par); // b_interior_pid1
	DefineData("vector", d[4], pid2, neighbor_partitions, par); // b_interior_pid2
	DefineData("vector", d[5], pid1, neighbor_partitions, par); // x_interior_pid1
	DefineData("vector", d[6], pid2, neighbor_partitions, par); // x_interior_pid2
	DefineData("scalar", d[7], pid1, neighbor_partitions, par); // rho_old_pid1
	DefineData("scalar", d[8], pid2, neighbor_partitions, par); // rho_old_pid2
	assert(8 == NUM_OF_FORLOOP_INPUTS -1);

	GetNewJobID(&j, 3);

	// Init, pid = 1
	before.clear();
	after.clear();
	after.insert(j[2]);
	read.clear();
	write.clear();
	write.insert(d[1]);
	write.insert(d[3]);
	write.insert(d[5]);
	SpawnComputeJob("Init", j[0], read, write, before, after, par);

	// Init, pid = 2
	before.clear();
	after.clear();
	after.insert(j[2]);
	read.clear();
	write.clear();
	write.insert(d[2]);
	write.insert(d[4]);
	write.insert(d[6]);
	SpawnComputeJob("Init", j[1], read, write, before, after, par);

	// Porject_Forloop_Condition
	before.clear();
	before.insert(j[0]);
	before.insert(j[1]);
	after.clear();
	read.clear();
	read.insert(d[0]);
	write.clear();
	param_idset.clear();
	for (int i=0;i<NUM_OF_FORLOOP_INPUTS;i++) param_idset.insert(d[i]);
	param_idset.insert(1); // insert iteration
	par.set_idset(param_idset);
	SpawnComputeJob("Project_Forloop_Condition", j[2], read, write, before, after, par);

	printf("Completed Proj_MainProjection\n");
};
