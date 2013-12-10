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

#include "shared/nimbus.h"
#include "projection_driver.h"
#include "projection_example.h"
#include "app.h"
#include "data_impl.h"
#include "job_impl.h"
#include "PCG_Sparse_Solver.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

#define LOCAL_N 1056
#define INTERIOR_N 1024
// HERE we set desired_iterations = pcg.maximum_iterations, since global_n >> pcg.maximum_iterations
#define DESIRED_ITERATIONS 200
#define GLOBAL_TOLERANCE 2.8e-10

Init::Init(Application *app) {
	set_application(app);
}
;

Job* Init::Clone() {
	printf("Cloning Init job\n");
	return new Init(application());
};

void Init::Execute(Parameter params, const DataArray& da) {	
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Init job starts on worker %d.\n", projection_app->_rankID);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	PhysBAM::ProjectionInternalData< PhysBAM::VECTOR<float,2> >* internal = app_driver->projection_internal_data;
	
	// load data
	Sparse_Matrix *A_out = reinterpret_cast<Sparse_Matrix*>(da[1]);
	Sparse_Matrix *AC_out = reinterpret_cast<Sparse_Matrix*>(da[2]);
	PCG_Vector *b_interior_out = reinterpret_cast<PCG_Vector*>(da[3]);
	PCG_Vector *z_interior_out = reinterpret_cast<PCG_Vector*>(da[4]);
	PCG_Vector *p_interior_out = reinterpret_cast<PCG_Vector*>(da[5]);
	PCG_Vector *p_ghost_out = reinterpret_cast<PCG_Vector*>(da[6]);
	PCG_Vector *temp_interior_out = reinterpret_cast<PCG_Vector*>(da[7]);
	PCG_Vector *x_interior_out = reinterpret_cast<PCG_Vector*>(da[8]);
	
	// execution
	app_driver->pcg_mpi->Initialize_Datatypes();
	int local_n=(*app_driver->projection_data->A_array)(1).n;
	internal->temp = new VECTOR_ND<T>(local_n,false);
	internal->p = new VECTOR_ND<T>(local_n,false);	
	VECTOR_ND<T>* x_interior_ptr = new VECTOR_ND<T>;
	x_interior_ptr->Set_Subvector_View(*app_driver->projection_data->x, app_driver->pcg_mpi->partition.interior_indices);
	VECTOR_ND<T>* b_interior_ptr = new VECTOR_ND<T>;
	b_interior_ptr->Set_Subvector_View((*app_driver->projection_data->b_array)(1), app_driver->pcg_mpi->partition.interior_indices);
	internal->p_interior = new VECTOR_ND<T>;
	internal->p_interior->Set_Subvector_View(*app_driver->projection_internal_data->p, app_driver->pcg_mpi->partition.interior_indices);	

	VECTOR_ND<float>& x = (*app_driver->projection_data->x); 
	VECTOR_ND<float>& temp = (*app_driver->projection_internal_data->temp);
	VECTOR_ND<float>& b_interior = (*b_interior_ptr);
	VECTOR_ND<float> temp_interior;
	temp_interior.Set_Subvector_View(temp, app_driver->pcg_mpi->partition.interior_indices);
	SPARSE_MATRIX_FLAT_NXN<float>& A = (*app_driver->projection_data->A_array)(1);
	app_driver->pcg_mpi->Fill_Ghost_Cells(x);
	A.Times(x,temp);b_interior-=temp_interior;
	delete A.C;A.C=A.Create_Submatrix(app_driver->pcg_mpi->partition.interior_indices);	   
	A.C->In_Place_Incomplete_Cholesky_Factorization(1, 0.97, 1e-8, 1e-8);
	
	// output data
	A_out->matrix_ = &A;
	AC_out->matrix_ = A.C;
	b_interior_out->vec_ = b_interior_ptr;
	z_interior_out->vec_ = new VECTOR_ND<T>(INTERIOR_N,false);
	p_interior_out->vec_ = internal->p_interior;
	p_ghost_out->vec_ = internal->p;
	temp_interior_out->vec_ = new VECTOR_ND<T>;	
	x_interior_out->vec_ = x_interior_ptr;
	
	dbg(DBG_PROJ, "||Init job finishes on worker %d.\n", projection_app->_rankID);
};

Project_Forloop_Condition::Project_Forloop_Condition(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Condition::Clone() {
	
	return new Project_Forloop_Condition(application());
};

void Project_Forloop_Condition::Execute(Parameter params, const DataArray& input_data) {
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Forloop_Condition job starts on worker %d.\n", projection_app->_rankID);	
	
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> params_data;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;		
	Parameter par;
	IDSet<param_id_t> param_idset;

	// parameters
	IDSet<param_id_t>::IDSetContainer::iterator it;
	IDSet<param_id_t> temp_set = params.idset();
	for (it = temp_set.begin(); it != temp_set.end(); it++) {
		params_data.push_back(*it);
	}
	
	// load data
	PartialNorm *residual_in = reinterpret_cast<PartialNorm*>(input_data[0]);
	
	// execution	
	int iteration = params_data.back();
	printf("Jia: forloop check passed, iter = %d, res = %lf", iteration, residual_in->norm_);	
	if(iteration == 1 || (iteration < DESIRED_ITERATIONS && residual_in->norm_ > GLOBAL_TOLERANCE)) {
		GetNewJobID(&j, 19);

		// Project_Forloop_Part1, pid = 1
		read.clear();
		read.insert(params_data[4]); // AC_pid1
		read.insert(params_data[6]); // b_interior_pid1
		write.clear();
		write.insert(params_data[8]); //z_interior_pid1
		write.insert(params_data[10]); // local_dot_sum_pid1
		before.clear();
		after.clear(); after.insert(j[2]);
		SpawnComputeJob("Project_Forloop_Part1", j[0], read, write, before, after, par);

		// Project_Forloop_Part1, pid = 2
		read.clear();
		read.insert(params_data[5]); // AC_pid2
		read.insert(params_data[7]); // b_interior_pid2
		write.clear();
		write.insert(params_data[9]); // z_interior_pid2
		write.insert(params_data[11]); // local_dot_sum_pid2
		before.clear();
		after.clear();
		after.insert(j[2]);
		after.insert(j[12]); // SpawnCopyJob
		SpawnComputeJob("Project_Forloop_Part1", j[1], read, write, before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[1]);
		after.clear();
		after.insert(j[2]);
		SpawnCopyJob(j[12], params_data[11], params_data[0], before, after, par);

		// Global_Sum
		read.clear();
		read.insert(params_data[10]); // local_dot_sum_pid1
		read.insert(params_data[0]); // local_dot_sum_pid2, CopyJob instance
		write.clear();
		write.insert(params_data[10]); // local_dot_sum
		before.clear();
		before.insert(j[0]); // Project_Forloop_Part1, pid = 1
		before.insert(j[1]); // Project_Forloop_Part1, pid = 2
		before.insert(j[12]); // SpawnCopyJob j[12]
		after.clear();
		after.insert(j[3]); // Project_Forloop_Part2, pid = 1
		after.insert(j[4]); // Project_Forloop_Part2, pid = 2
		after.insert(j[13]); // SpawnCopyJob j[13]
		SpawnComputeJob("Global_Sum", j[2], read, write, before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[2]);
		after.clear();
		after.insert(j[4]);
		SpawnCopyJob(j[13], params_data[10], params_data[11], before, after, par);
	
		// Project_Forloop_Part2, pid = 1
		read.clear();
		read.insert(params_data[8]); // z_interior_pid1
		read.insert(params_data[10]); // local_dot_sum_pid1
		read.insert(params_data[12]); // rho_pid1
		read.insert(params_data[14]); // p_interior_pid1
		write.clear();		
		write.insert(params_data[12]); // rho_pid1
		write.insert(params_data[14]); // p_interior_pid1
		before.clear();
		before.insert(j[2]); // Global_Sum
		after.clear();
		after.insert(j[5]); // Project_Forloop_Part3
		after.insert(j[14]); // SpawnCopyJob
		param_idset.clear(); param_idset.insert(iteration);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Part2", j[3], read, write, before, after, par);

		// Project_Forloop_Part2, pid = 2
		read.clear();
		read.insert(params_data[9]); // z_interior_pid2
		read.insert(params_data[11]); // local_dot_sum_pid2
		read.insert(params_data[13]); // rho_pid2
		read.insert(params_data[15]); // p_interior_pid2
		write.clear();		
		write.insert(params_data[13]); // rho_pid2
		write.insert(params_data[15]); // p_interior_pid2
		before.clear();
		before.insert(j[2]); // Global_Sum
		before.insert(j[13]); // SpawnCopyJob
		after.clear();
		after.insert(j[6]); // Project_Forloop_Part3
		after.insert(j[15]); // SpawnCopyJob
		param_idset.clear(); param_idset.insert(iteration);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Part2", j[4], read, write, before, after, par);
		
		// SpawnCopyJob
		before.clear();
		before.insert(j[3]);
		after.clear();
		after.insert(j[5]);
		after.insert(j[6]);
		SpawnCopyJob(j[14], params_data[14], params_data[17], before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[4]);
		after.clear();
		after.insert(j[5]);
		after.insert(j[6]);
		SpawnCopyJob(j[15], params_data[15], params_data[16], before, after, par);

		// Project_Forloop_Part3, pid = 1
		read.clear();
		read.insert(params_data[2]); // A_pid1
		read.insert(params_data[14]); // p_interior_pid1
		read.insert(params_data[16]); // p_ghost_pid1
		write.clear();
		write.insert(params_data[18]); // tmp_interior_pid1
		write.insert(params_data[10]); // local_dot_prod_p_temp_pid1
		before.clear();
		before.insert(j[3]);
		before.insert(j[4]);
		before.insert(j[14]); // SpawnCopyJob
		before.insert(j[15]); // SpawnCopyJob
		after.clear();
		after.insert(j[7]);
		SpawnComputeJob("Project_Forloop_Part3", j[5], read, write, before, after, par);

		// Project_Forloop_Part3, pid = 2
		read.clear();
		read.insert(params_data[3]); // A_pid2
		read.insert(params_data[15]); // p_interior_pid2
		read.insert(params_data[17]); // p_ghost_pid2
		write.clear();
		write.insert(params_data[19]); // tmp_interior_pid2
		write.insert(params_data[11]); // local_dot_prod_p_temp_pid2
		before.clear();
		before.insert(j[3]);
		before.insert(j[4]);
		before.insert(j[14]); // SpawnCopyJob
		before.insert(j[15]); // SpawnCopyJob
		after.clear();
		after.insert(j[7]);
		after.insert(j[16]);
		SpawnComputeJob("Project_Forloop_Part3", j[6], read, write, before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[6]);
		after.clear();
		after.insert(j[7]);
		SpawnCopyJob(j[16], params_data[11], params_data[0], before, after, par);

		// Global_Sum
		read.clear();
		read.insert(params_data[10]); // local_dot_sum_pid1
		read.insert(params_data[0]); // local_dot_sum_pid2, CopyJob instance
		write.clear();
		write.insert(params_data[10]); // global_sum
		before.clear();
		before.insert(j[5]); // Project_Forloop_Part3, pid = 1
		before.insert(j[6]); // Project_Forloop_Part3, pid = 2
		before.insert(j[16]); // SpawnCopyJob
		after.clear();
		after.insert(j[8]); // Project_Forloop_Part4, pid = 1
		after.insert(j[9]); // Project_Forloop_Part4, pid = 2
		after.insert(j[17]); // SpawnJobCopy
		SpawnComputeJob("Global_Sum", j[7], read, write, before, after, par);

		// SpawnCopyJob		
		before.clear();
		before.insert(j[7]);
		after.clear();
		after.insert(j[8]);
		after.insert(j[9]);
		SpawnCopyJob(j[17], params_data[10], params_data[11], before, after, par);

		// Project_Forloop_Part4, pid = 1
		read.clear();
		read.insert(params_data[20]); // x_interior_pid1
		read.insert(params_data[14]); // p_interior_pid1
		read.insert(params_data[6]); // b_interior_pid1
		read.insert(params_data[18]); // tmp_interior_pid1
		read.insert(params_data[12]); // rho_pid1
		read.insert(params_data[10]); // local_dot_sum_pid1
		write.clear();
		write.insert(params_data[22]); // local_max_abs_pid1
		write.insert(params_data[20]); // x_interior_pid1
		write.insert(params_data[6]); // b_interior_pid1		
		before.clear();
		before.insert(j[7]);
		before.insert(j[17]);
		after.clear();
		after.insert(j[10]);
		SpawnComputeJob("Project_Forloop_Part4", j[8], read, write, before, after, par);

		// Project_Forloop_Part4, pid = 2
		read.clear();
		read.insert(params_data[21]); // x_interior_pid2
		read.insert(params_data[15]); // p_interior_pid2		
		read.insert(params_data[7]); // b_interior_pid2
		read.insert(params_data[19]); // tmp_interior_pid2
		read.insert(params_data[13]); // rho_pid2
		read.insert(params_data[11]); // local_dot_sum_pid2				
		write.clear();
		write.insert(params_data[23]); // local_max_abs_pid2
		write.insert(params_data[21]); // x_interior_pid2
		write.insert(params_data[7]); // b_interior_pid2
		before.clear();
		before.insert(j[7]);
		before.insert(j[17]); // SpawnCopyJob
		after.clear();
		after.insert(j[10]);
		after.insert(j[18]); // SpawnCopyJobprojection_internal_data
		SpawnComputeJob("Project_Forloop_Part4", j[9], read, write, before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[9]);
		after.clear();
		after.insert(j[10]);
		SpawnCopyJob(j[18], params_data[23], params_data[0], before, after, par);

		// Global_Max_Abs
		read.clear();
		read.insert(params_data[22]); // local_max_abs_pid1
		read.insert(params_data[0]); // local_max_abs_pid2, Copy instance in pid1
		write.clear();
		write.insert(params_data[22]); // residual
		before.clear();
		before.insert(j[8]);
		before.insert(j[9]);
		before.insert(j[18]); // SpawnCopyJob
		after.clear();
		after.insert(j[11]);
		SpawnComputeJob("Global_Max", j[10], read, write, before, after, par);

		// Project_Forloop_Condition
		read.clear();
		read.insert(params_data[22]);
		write.clear();
		before.clear();
		before.insert(j[10]);
		after.clear();
		param_idset.clear();
		for (unsigned i=0;i<params_data.size()-1;i++) param_idset.insert(params_data[i]);
		param_idset.insert(iteration + 1);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Condition", j[11], read, write, before, after, par);
	}
	else {
		GetNewJobID(&j, 2);
		read.clear();
		read.insert(params_data[19]);
		write.clear();
		before.clear();
		after.clear();
		SpawnComputeJob("Finish", j[0], read, write, before, after, par);
		read.clear();
		read.insert(params_data[20]);
		write.clear();
		before.clear();
		after.clear();
		SpawnComputeJob("Finish", j[1], read, write, before, after, par);		
	}	
	dbg(DBG_PROJ, "||Forloop_Condition job finishes on worker %d.\n", projection_app->_rankID);
};

Project_Forloop_Part1::Project_Forloop_Part1(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part1::Clone() {
	std::cout << "Cloning Project_Forloop_Part1 job!\n";
	return new Project_Forloop_Part1(application());
};

void Project_Forloop_Part1::Execute(Parameter params, const DataArray& da) {	
	
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Forloop_Part1 job starts on worker %d.\n", projection_app->_rankID);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	
	// load data	
	Sparse_Matrix *AC_in = reinterpret_cast<Sparse_Matrix*>(da[0]);
	PCG_Vector *b_interior_in = reinterpret_cast<PCG_Vector*>(da[1]);
	PCG_Vector *z_interior_out = reinterpret_cast<PCG_Vector*>(da[2]);
	PartialNorm *local_dot_sum_out = reinterpret_cast<PartialNorm*>(da[3]);
	
	//execution
	VECTOR_ND<T>& z_interior = (*z_interior_out->vec_);
	VECTOR_ND<T>& b_interior = (*b_interior_in->vec_);
	VECTOR_ND<T> temp(LOCAL_N, false);
	VECTOR_ND<T> temp_interior;
	temp_interior.Set_Subvector_View(temp, app_driver->pcg_mpi->partition.interior_indices);
	AC_in->matrix_->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
	AC_in->matrix_->Solve_Backward_Substitution(temp_interior,z_interior,false,true); // diagonal is inverted to save on divides

	local_dot_sum_out->norm_ = VECTOR_ND<float>::Dot_Product_Double_Precision(z_interior, b_interior);
	dbg(DBG_PROJ, "||Forloop_Part1 job finishes on worker %d.\n", projection_app->_rankID);
};

Project_Forloop_Part2::Project_Forloop_Part2(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part2::Clone() {
	std::cout << "Cloning Project_Forloop_Part2 job!\n";
	return new Project_Forloop_Part2(application());
};

void Project_Forloop_Part2::Execute(Parameter params, const DataArray& da) {	
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Forloop_Part2 job starts on worker %d.\n", projection_app->_rankID);	
	
	// load data	
	PCG_Vector *z_interior_in = reinterpret_cast<PCG_Vector*>(da[0]);
	PartialNorm *global_sum_in = reinterpret_cast<PartialNorm*>(da[1]);
	PartialNorm *rho_inout = reinterpret_cast<PartialNorm*>(da[2]);
	PCG_Vector *p_interior_inout = reinterpret_cast<PCG_Vector*>(da[3]);
	
	// execution
	int iteration = *(params.idset().begin());
	VECTOR_ND<float>& z_interior = (*z_interior_in->vec_);
	//VECTOR_ND<float>& b_interior = (*internal->b_interior);
	VECTOR_ND<float>& p_interior = (*p_interior_inout->vec_);
	
	float beta = 0;	
	if (iteration == 1) {
		p_interior = z_interior;
	} else {
		beta = global_sum_in->norm_ / rho_inout->norm_;		
		for(int i=1;i<=INTERIOR_N;i++)
			p_interior(i) = z_interior(i) + beta * p_interior(i);
	}
	
	//write data
	rho_inout->norm_ = global_sum_in->norm_;
	dbg(DBG_PROJ, "||Forloop_Part2 job finishes on worker %d.\n", projection_app->_rankID);
};

Project_Forloop_Part3::Project_Forloop_Part3(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part3::Clone() {
	std::cout << "Cloning Project_Forloop_Part3 job!\n";
	return new Project_Forloop_Part3(application());
};

void Project_Forloop_Part3::Execute(Parameter params, const DataArray& da) {	
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Forloop_Part3 job starts on worker %d.\n", projection_app->_rankID);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	PhysBAM::ProjectionInternalData< PhysBAM::VECTOR<float,2> >* internal = app_driver->projection_internal_data;
	
	// load data
	Sparse_Matrix *A = reinterpret_cast<Sparse_Matrix*>(da[0]);
	PCG_Vector *p_interior_in = reinterpret_cast<PCG_Vector*>(da[1]);
	PCG_Vector *temp_interior_out = reinterpret_cast<PCG_Vector*>(da[3]);
	PartialNorm *local_dot_sum_out = reinterpret_cast<PartialNorm*>(da[4]);
	
	// execution
	VECTOR_ND<float>& p = (*internal->p);
	app_driver->pcg_mpi->Fill_Ghost_Cells(p);
	//int color = 1;
	//SPARSE_MATRIX_FLAT_NXN<float>& A = (*app_driver->projection_data->A_array)(color);
	VECTOR_ND<float>& temp = (*internal->temp);
	VECTOR_ND<float>& temp_interior = *(temp_interior_out->vec_);
	VECTOR_ND<float>& p_interior = *(p_interior_in->vec_);
	A->matrix_->Times(p, temp);
	temp_interior.Set_Subvector_View(temp, app_driver->pcg_mpi->partition.interior_indices);
	
	local_dot_sum_out->norm_ = VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior);
	
	dbg(DBG_PROJ, "||Forloop_Part3 job finishes on worker %d.\n", projection_app->_rankID);
};

Project_Forloop_Part4::Project_Forloop_Part4(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part4::Clone() {
	std::cout << "Cloning Project_Forloop_Part4 job!\n";
	return new Project_Forloop_Part4(application());
};

void Project_Forloop_Part4::Execute(Parameter params, const DataArray& da) {	
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Forloop_Part4 job starts on worker %d.\n", projection_app->_rankID);
	
	// load data
	PCG_Vector *x_interior_inout = reinterpret_cast<PCG_Vector*>(da[0]);
	PCG_Vector *p_interior_in = reinterpret_cast<PCG_Vector*>(da[1]);
	PCG_Vector *b_interior_inout = reinterpret_cast<PCG_Vector*>(da[2]);
	PCG_Vector *temp_interior_in = reinterpret_cast<PCG_Vector*>(da[3]);
	PartialNorm *rho_in = reinterpret_cast<PartialNorm*>(da[4]);
	PartialNorm *local_dot_sum_in = reinterpret_cast<PartialNorm*>(da[5]);
	PartialNorm *local_max_abs_out = reinterpret_cast<PartialNorm*>(da[6]);
		
	// execution
	VECTOR_ND<float>& p_interior = (*p_interior_in->vec_);
	VECTOR_ND<float>& x_interior = (*x_interior_inout->vec_);
	VECTOR_ND<float>& b_interior = (*b_interior_inout->vec_);
	VECTOR_ND<float>& temp_interior = (*temp_interior_in->vec_);
	
	float alpha=(float)(rho_in->norm_/local_dot_sum_in->norm_);
	for(int i=1;i<=INTERIOR_N;i++) {
		x_interior(i) += alpha * p_interior(i);
		b_interior(i) -= alpha * temp_interior(i);
	}
	local_max_abs_out->norm_ = b_interior.Max_Abs();	
	dbg(DBG_PROJ, "||Forloop_Part4 job finishes on worker %d.\n", projection_app->_rankID);
};

Global_Sum::Global_Sum(Application* app) {
	set_application(app);
};

Job * Global_Sum::Clone() {
	std::cout << "Cloning Global_Sum job!\n";
	return new Global_Sum(application());
};

void Global_Sum::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Global_Sum job\n";
	PartialNorm *d0 = reinterpret_cast<PartialNorm*>(da[0]); // part1
	PartialNorm *d1 = reinterpret_cast<PartialNorm*>(da[1]); // part2
	PartialNorm *d2 = reinterpret_cast<PartialNorm*>(da[2]); // sum
	
	d2->norm_ = d0->norm_ + d1->norm_;	
	std::cout << "Completed Global_Sum job\n";
};

Global_Max::Global_Max(Application* app) {
	set_application(app);
};

Job * Global_Max::Clone() {
	std::cout << "Cloning Global_Max job!\n";
	return new Global_Max(application());
};

void Global_Max::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Global_Max job\n";
	PartialNorm *d0 = reinterpret_cast<PartialNorm*>(da[0]); // part1
	PartialNorm *d1 = reinterpret_cast<PartialNorm*>(da[1]); // part2
	PartialNorm *d2 = reinterpret_cast<PartialNorm*>(da[2]); // max
	
	if(d0->norm_ > d1->norm_)
		d2->norm_ = d0->norm_;
	else
		d2->norm_ = d1->norm_;
	std::cout << "Completed Global_Max job\n";
};

Finish::Finish(Application* app) {
	set_application(app);
}

void Finish::Execute(Parameter params, const DataArray& da) {	
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Finish job starts on worker %d.\n", projection_app->_rankID);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	
	VECTOR_ND<T>& x = (*app_driver->projection_data->x); 
	app_driver->pcg_mpi->Fill_Ghost_Cells(x);
	
	//=============================== the below code MUST be executed after PCG_SPARSE_MPI::Parellel_Solve ===================
	app_driver->WindUpForOneRegion();
	app_driver->ApplyPressureAndFinish();
	projection_app->FinishMain();
	dbg(DBG_PROJ, "||Finish job finishes on worker %d.\n", projection_app->_rankID);
	dbg(DBG_PROJ, "||||Projection finishes successfully!!!\n");
}

Job* Finish::Clone() {
	return new Finish(application());
}
