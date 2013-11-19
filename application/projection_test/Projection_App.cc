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
}
;

void ProjectionApp::Load() {
	printf("Worker beginning to load application\n");

	/* Declare and initialize data, jobs and policies. */

	RegisterJob("main", new Main(this));
	RegisterJob("Init", new Init(this));
	RegisterJob("Project_Forloop_Condition",
			new Project_Forloop_Condition(this));
	RegisterJob("Project_Forloop_Part1", new Project_Forloop_Part1(this));
	RegisterJob("Project_Forloop_Part2", new Project_Forloop_Part2(this));
	RegisterJob("Project_Forloop_Part3", new Project_Forloop_Part3(this));
	RegisterJob("Project_Forloop_Part4", new Project_Forloop_Part4(this));
	RegisterJob("Global_Sum", new Global_Sum(this));
	RegisterJob("Global_Max_Abs", new Global_Max_Abs(this));

	RegisterData("vector", new Vec(LEN));
	RegisterData("scalar", new Vec(1));
	RegisterData("matrix", new Vec((LEN * LEN)));

	printf("Finished creating job and data definitions\n");
}

Main::Main(Application *app) {
	set_application(app);
}
;

Job* Main::Clone() {
	printf("Cloning main job\n");
	return new Main(application());
}
;

void Main::Execute(Parameter params, const DataArray& da) {
	printf("Begin main\n");
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet < logical_data_id_t > read, write;
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

	for (int i = 1; i <= 100; i++) {
		before.clear();
		before.insert(j[i - 1]);
		after.clear();
		after.insert(j[i + 1]);
		read.clear();
		write.clear();
		SpwanComputeJob("Proj_PrepareForProj", j[i], read, write, before,
				after, par);
	}

	printf("Completed main\n");
}
;

Proj_Initialize::Proj_Initialize(Application *app) {
	set_application(app);
}
;

Job* Proj_Initialize::Clone() {
	printf("Cloning Proj_Initialize job\n");
	return new Proj_Initialize(application());
}
;

void Proj_Initialize::Execute(Parameter params, const DataArray& da) {
	printf("Begin Proj_Initialize\n");
	printf("Initialize example\n");
	PROJECTION_EXAMPLE<TV>* example = new PROJECTION_EXAMPLE<TV> ();
	int scale = SCALE_VAL;
	RANGE<TV> range(TV(), TV::All_Ones_Vector() * 0.5);
	range.max_corner(2) = 1;
	TV_INT counts = TV_INT::All_Ones_Vector() * scale / 2;
	counts(2) = scale;
	example->Initialize_Grid(counts, range);
	example->restart = RESTART_VAL;
	example->last_frame = 100;

	printf("Initialize driver\n");
	// setup time
	if (example->restart)
		current_frame = example->restart;
	else
		current_frame = example->first_frame;
	time = example->Time_At_Frame(current_frame);

	// will map mpi_grid into Nimbus job dependency
	/*
	 if (example.mpi_grid)
	 example.mpi_grid->Initialize(example.domain_boundary);
	 example.projection.elliptic_solver->mpi_grid = example.mpi_grid;
	 if (example.mpi_grid)
	 example.boundary = new BOUNDARY_MPI<GRID<TV> , T> (example.mpi_grid,
	 example.boundary_scalar);
	 else
	 example.boundary = &example.boundary_scalar;
	 */

	// setup grids and velocities
	example->projection.Initialize_Grid(example.mac_grid);
	example->face_velocities.Resize(example.mac_grid);
	example->Initialize_Fields();

	// setup laplace
	example->projection.elliptic_solver->Set_Relative_Tolerance(1e-11);
	example->projection.elliptic_solver->pcg.Set_Maximum_Iterations(200);
	example->projection.elliptic_solver->pcg.evolution_solver_type
			= krylov_solver_cg;
	example->projection.elliptic_solver->pcg.cg_restart_iterations = 40;

	// setup domain boundaries
	VECTOR<VECTOR<bool, 2> , TV::dimension> constant_extrapolation;
	constant_extrapolation.Fill(VECTOR<bool, 2>::Constant_Vector(true));
	example->boundary->Set_Constant_Extrapolation(constant_extrapolation);
	example->Set_Boundary_Conditions(time); // get so CFL is correct

	printf("Completed Proj_Initialize\n");
}
;

Proj_PrepareForProj::Proj_PrepareForProj(Application *app) {
	set_application(app);
}
;

Job* Proj_PrepareForProj::Clone() {
	printf("Cloning Proj_PrepareForProj job\n");
	return new Proj_PrepareForProj(application());
}
;

// read set: current_frame, time
// write set: time
void Proj_PrepareForProj::Execute(Parameter params, const DataArray& da) {
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet < logical_data_id_t > read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	partition_id_t pid1 = 1, pid2 = 2;
	Parameter par;
	IDSet<param_id_t> param_idset;

	printf("Begin Proj_PrepareForProj\n");
	T target_time = example->Time_At_Frame(current_frame + 1);
	T dt = target_time - time;
	example->Set_Boundary_Conditions(time + dt);
	example->projection.p *= dt; // rescale pressure for guess
	example->projection.Compute_Divergence(typename INTERPOLATION_POLICY<GRID<
			TV> >::FACE_LOOKUP(example->face_velocities),
			example->projection.elliptic_solver); // find f - divergence of the velocity

	LAPLACE_UNIFORM<GRID<TV> >* laplace = example->projection.elliptic_solver;
	laplace->Find_Solution_Regions(); // flood fill
	typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR
			T_ARRAYS_SCALAR;
	typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
	typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
	typedef typename GRID<TV>::VECTOR_INT TV_INT;

	ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array(
			laplace->number_of_regions);
	T_ARRAYS_INT cell_index_to_matrix_index(laplace->grid.Domain_Indices(1));
	ARRAY<int, VECTOR<int, 1> > filled_region_cell_count(-1,
			laplace->number_of_regions);
	ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array(laplace->number_of_regions);
	ARRAY<VECTOR_ND<T> > b_array(laplace->number_of_regions);
	for (CELL_ITERATOR iterator(laplace->grid, 1); iterator.Valid(); iterator.Next())
		filled_region_cell_count(laplace->filled_region_colors(
				iterator.Cell_Index()))++;
	for (int color = 1; color <= laplace->number_of_regions; color++)
		if (laplace->filled_region_touches_dirichlet(color)
				|| laplace->solve_neumann_regions) {
			matrix_index_to_cell_index_array(color).Resize(
					filled_region_cell_count(color));
		}
	filled_region_cell_count.Fill(0); // reusing this array in order to make the indirection arrays
	DOMAIN_ITERATOR_THREADED_ALPHA<LAPLACE_UNIFORM<GRID<TV> > , TV>
			threaded_iterator(laplace->grid.Domain_Indices(1),
					laplace->thread_queue, 1, 1, 2, 1);

	ARRAY<int, TV_INT> domain_index(laplace->grid.Domain_Indices(1), false);

	for (int i = 1; i <= threaded_iterator.domains.m; i++) {
		RANGE<TV_INT> interior_domain(threaded_iterator.domains(i));
		interior_domain.max_corner -= TV_INT::All_Ones_Vector();
		interior_domain.min_corner += TV_INT::All_Ones_Vector();
		for (CELL_ITERATOR iterator(laplace->grid, interior_domain); iterator.Valid(); iterator.Next())
			domain_index(iterator.Cell_Index()) = i;
	}
	ARRAY<ARRAY<INTERVAL<int> > > interior_indices(laplace->number_of_regions);
	ARRAY<ARRAY<ARRAY<INTERVAL<int> > > > ghost_indices(
			laplace->number_of_regions);
	for (int color = 1; color <= laplace->number_of_regions; color++) {
		interior_indices(color).Resize(threaded_iterator.number_of_domains);
		ghost_indices(color).Resize(threaded_iterator.number_of_domains);
		for (int i = 1; i <= threaded_iterator.domains.m; i++)
			ghost_indices(color)(i).Resize(2 * TV::dimension);
	}
	//Begin laplace->laplace_mpi->Find_Matrix_Indices(filled_region_cell_count, cell_index_to_matrix_index, matrix_index_to_cell_index_array);

	RANGE<TV_INT> domain = laplace->grid.Domain_Indices(1);
	laplace->Find_A(domain, A_array, b_array, filled_region_cell_count,
			cell_index_to_matrix_index);

	GetNewJobID(&j, number_of_regions);
	for (int color = 1; color <= number_of_regions; color++)
		if (filled_region_cell_count(color) > 0
			&& (filled_region_touches_dirichlet(color)
				|| solve_neumann_regions)) {
			//Solve_Subregion(interior_indices(color), ghost_indices(color), matrix_index_to_cell_index_array(color), A_array(color), b_array(color), color, &domain_index);
			before.clear();
			after.clear();
			read.clear();
			write.clear();
			SpawnComputeJob("Prof_PrepareForOneRegion", j[color], read, write, before, after, par);
		}
	printf("Completed Proj_PrepareForProj\n");
}
;

Prof_PrepareForOneRegion::Prof_PrepareForOneRegion(Application *app) {
	set_application(app);
}
;

Job* Prof_PrepareForOneRegion::Clone() {
	printf("Cloning Prof_PrepareForOneRegion job\n");
	return new Prof_PrepareForOneRegion(application());
}
;

// read set:interior_indices, ghost_indices, matrix_index_to_cell_index, A, b, color, domain_index
// write set:
void Prof_PrepareForOneRegion::Execute(Parameter params, const DataArray& da) {
	printf("Begin Prof_PrepareForOneRegion\n");
	int number_of_unknowns = matrix_index_to_cell_index.m;
	A.Negate();
	b *= (T) - 1;
	VECTOR_ND<T> x(number_of_unknowns), q, s, r, k, z;
	for (int i = 1; i <= number_of_unknowns; i++)
		x(i) = u(matrix_index_to_cell_index(i));
	laplace->Find_Tolerance(b); // needs to happen after b is completely set up

	laplace_mpi->Solve(A, x, b, q, s, r, k, z, tolerance, color);

	for (int i = 1; i <= number_of_unknowns; i++) {
		TV_INT cell_index = matrix_index_to_cell_index(i);
		u(cell_index) = x(i);
	}
	printf("Completed Prof_PrepareForOneRegion\n");
}
;

Prof_AfterProj::Prof_AfterProj(Application *app) {
	set_application(app);
}
;

Job* Prof_AfterProj::Clone() {
	printf("Cloning Prof_AfterProj job\n");
	return new Prof_AfterProj(application());
}
;

// read set: example, dt, time
// write set:
void Prof_AfterProj::Execute(Parameter params, const DataArray& da) {
	printf("Begin Prof_AfterProj\n");

	example->Apply_Pressure(example->face_velocities, dt, time);

	example->projection.p *= (1 / dt); // unscale pressure
	time += dt;
	printf("Completed Prof_AfterProj\n");
}
;

Proj_MainProjection::Proj_MainProjection(Application *app) {
	set_application(app);
}
;

Job* Proj_MainProjection::Clone() {
	printf("Cloning Proj_MainProjection job\n");
	return new Proj_MainProjection(application());
}
;

void Proj_MainProjection::Execute(Parameter params, const DataArray& da) {
	printf("Begin Proj_MainProjection\n");
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet < logical_data_id_t > read, write;
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
	for (int i = 0; i < NUM_OF_FORLOOP_INPUTS; i++)
		param_idset.insert(d[i]);
	param_idset.insert(1); // insert iteration
	par.set_idset(param_idset);
	SpawnComputeJob("Project_Forloop_Condition", j[2], read, write, before,
			after, par);

	printf("Completed Proj_MainProjection\n");
}
;
