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
#include "PCG_Sparse_Solver.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

Init::Init(Application *app) {
	set_application(app);
}
;

Job* Init::Clone() {
	printf("Cloning Init job\n");
	return new Init(application());
}
;


void Init::Execute(Parameter params, const DataArray& da) {
	printf("Begin Init\n");
	App* projection_app = dynamic_cast<App*>(application());
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	
	app_driver->pcg_mpi.Initialize_Datatypes();
	int local_n=(*projection_data->A_array)(1).n;
	app_driver->projection_internal_data->interior_n=partition.interior_indices.Size()+1;
	app_driver->projection_internal_data->x_interior = new VECTOR_ND<T>;
	app_driver->projection_internal_data->x_interior->Set_Subvector_View(*projection_data->x, partition.interior_indices);
	app_driver->projection_internal_data->b_interior = new VECTOR_ND<T>;
	app_driver->projection_internal_data->b_interior->Set_Subvector_View((*projection_data->b_array)(1), partition.interior_indices);
	app_driver->projection_internal_data->p_interior = new VECTOR_ND<T>;
	app_driver->projection_internal_data->p_interior->Set_Subvector_View(*projection_internal_data->p, partition.interior_indices);
	app_driver->projection_internal_data->temp_interior = new VECTOR_ND<T>;
	app_driver->projection_internal_data->temp_interior->Set_Subvector_View(*projection_internal_data->temp, partition.interior_indices);
	app_driver->projection_internal_data->rho=0;
	app_driver->projection_internal_data->rho_old=0;
	app_driver->projection_internal_data->alpha=0;
	app_driver->projection_internal_data->beta=0;

	app_driver->projection_internal_data->global_n = app_driver->pcg_mpi.Global_Sum(interior_n);
	app_driver->projection_internal_data->global_tolerance = app_driver->pcg_mpi.Global_Max(app_driver->projection_data->tolerance);
	int desired_iterations=global_n;
	if (pcg.maximum_iterations)
		app_driver->projection_internal_data->desired_iterations=min(desired_iterations, pcg.maximum_iterations);

	VECTOR_ND<float>& x = (*app_driver->projection_data->x); 
	VECTOR_ND<float>& temp = (*app_driver->projection_internal_data->temp);
	VECTOR_ND<float>& b_interior = (*app_driver->projection_internal_data->b_interior);
	VECTOR_ND<float>& temp_interior = (*app_driver->projection_internal_data->temp_interior);
	SPARSE_MATRIX_FLAT_NXN<float>& A = (*projection_data->A_array)(1);
	app_driver->pcg_mpi.Fill_Ghost_Cells(x);
	A.Times(x,temp);b_interior-=temp_interior;
	delete A.C;A.C=A.Create_Submatrix(partition.interior_indices);
	A.C->In_Place_Incomplete_Cholesky_Factorization(
			app_driver->pcg_mpi.pcg.modified_incomplete_cholesky, 
			app_driver->pcg_mpi.pcg.modified_incomplete_cholesky_coefficient, 
			app_driver->pcg_mpi.pcg.preconditioner_zero_tolerance,
			app_driver->pcg_mpi.pcg.preconditioner_zero_replacement);
	printf("Completed Init\n");
};

Project_Forloop_Condition::Project_Forloop_Condition(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Condition::Clone() {
	std::cout << "Cloning Project_Forloop_Condition job!\n";
	return new Project_Forloop_Condition(application());
};

void Project_Forloop_Condition::Execute(Parameter params,
		const DataArray& input_data) {
	std::cout << "Executing the Project_Forloop_Condition job\n";
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d, da;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	partition_id_t pid1 = 1, pid2 = 2;
	Parameter par;
	IDSet<param_id_t> param_idset;

	// parameters
	IDSet<param_id_t>::IDSetContainer::iterator it;
	IDSet<param_id_t> temp_set = params.idset();
	for (it = temp_set.begin(); it != temp_set.end(); it++) {
		da.push_back(*it);
	}

	// input_data
	Vec *d0 = reinterpret_cast<Vec*>(input_data[0]); // residual

	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	PhysBAM::ProjectionInternalData< PhysBAM::VECTOR<float,2> >* internal = app_driver->projection_internal_data;
	
	// execution
	VECTOR_ND<float>& b_interior = (*internal->b_interior);
	internal->residual = app_driver->pcg_mpi.Global_Max(b_interior.Max_Abs());
	if(internal->iteration == 1 || (internal->iteration < DESIRED_ITERATIONS && internal->residual > internal->global_tolerance)) {
		printf("Jia: forloop check passed, iter = %d, res = %f\n", iteration, residual);
		GetNewJobID(&j, 19);

		// Project_Forloop_Part1, pid = 1
		read.clear();
		read.insert(da[1]); // A_pid1
		read.insert(da[3]); // b_interior_pid1
		write.clear();
		write.insert(d[14]); //z_interior_pid1
		write.insert(d[2]); // local_dot_prod_zb_pid1
		before.clear();
		after.clear(); after.insert(j[2]);
		SpawnComputeJob("Project_Forloop_Part1", j[0], read, write, before, after, par);

		// Project_Forloop_Part1, pid = 2
		read.clear();
		read.insert(da[2]); // A_pid2
		read.insert(da[4]); // b_interior_pid2
		write.clear();
		write.insert(d[15]); //z_interior_pid2
		write.insert(d[3]); // local_dot_prod_zb_pid2
		before.clear();
		after.clear();
		after.insert(j[2]);
		after.insert(j[13]); // SpawnCopyJob
		SpawnComputeJob("Project_Forloop_Part1", j[1], read, write, before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[1]);
		after.clear();
		after.insert(j[2]);
		SpawnCopyJob(j[12], d[3], d[16], before, after, par);

		// Global_Sum
		read.clear();
		read.insert(d[2]); // local_dot_prod_zb_pid1
		read.insert(d[16]); // local_dot_prod_zb_pid2, CopyJob instance
		write.clear();
		write.insert(d[4]); // rho
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
		SpawnCopyJob(j[13], d[4], d[17], before, after, par);

		// Project_Forloop_Part2, pid = 1
		read.clear();
		read.insert(d[4]); // rho
		read.insert(da[7]); // rho_old_pid1
		read.insert(d[14]); // z_interior_pid1
		read.insert(d[5]); // p_interior_pid1
		write.clear();
		write.insert(d[5]); // p_interior_pid1
		write.insert(da[7]); // rho_old_pid1
		before.clear();
		before.insert(j[2]); // Global_Sum
		after.clear();
		after.insert(j[5]); // Project_Forloop_Part3
		after.insert(j[14]); // SpawnCopyJob
		param_idset.clear(); param_idset.insert(internal->iteration);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Part2", j[3], read, write, before, after, par);

		// Project_Forloop_Part2, pid = 2
		read.clear();
		read.insert(d[17]); // rho, CopyJob instance
		read.insert(da[8]); // rho_old_pid2
		read.insert(d[15]); // z_interior_pid2
		read.insert(d[6]); // p_interior_pid2
		write.clear();
		write.insert(d[6]); // p_interior_pid2
		write.insert(da[8]); // rho_old_pid2
		before.clear();
		before.insert(j[2]); // Global_Sum
		before.insert(j[13]); // SpawnCopyJob
		after.clear();
		after.insert(j[6]); // Project_Forloop_Part3
		after.insert(j[15]); // SpawnCopyJob
		param_idset.clear(); param_idset.insert(internal->iteration);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Part2", j[4], read, write, before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[3]);
		after.clear();
		after.insert(j[6]);
		SpawnCopyJob(j[14], d[5], d[19], before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[4]);
		after.clear();
		after.insert(j[5]);
		SpawnCopyJob(j[15], d[6], d[18], before, after, par);

		// Project_Forloop_Part3, pid = 1
		read.clear();
		read.insert(da[1]); // A_pid1
		read.insert(d[5]); // p_interior_pid1
		read.insert(d[18]); // p_interior_pid2, CopyJob instance
		write.clear();
		write.insert(d[0]); // tmp_interior_pid1
		write.insert(d[9]); // local_dot_prod_p_temp_pid1
		before.clear();
		before.insert(j[3]);
		before.insert(j[4]);
		before.insert(j[15]); // SpawnCopyJob
		after.clear();
		after.insert(j[7]);
		SpawnComputeJob("Project_Forloop_Part3", j[5], read, write, before, after, par);

		// Project_Forloop_Part3, pid = 2
		read.clear();
		read.insert(da[2]); // A_pid2
		read.insert(d[6]); // p_interior_pid2
		read.insert(d[19]); // p_interior_pid1, CopyJob instance
		write.clear();
		write.insert(d[1]); // tmp_interior_pid2
		write.insert(d[10]); // local_dot_prod_p_temp_pid2
		before.clear();
		before.insert(j[3]);
		before.insert(j[4]);
		before.insert(j[14]); // SpawnCopyJob
		after.clear();
		after.insert(j[7]);
		after.insert(j[16]);
		SpawnComputeJob("Project_Forloop_Part3", j[6], read, write, before, after, par);

		// SpawnCopyJob
		before.clear();
		before.insert(j[6]);
		after.clear();
		after.insert(j[7]);
		SpawnCopyJob(j[16], d[10], d[20], before, after, par);

		// Global_Sum
		read.clear();
		read.insert(d[9]); // local_dot_prod_p_temp_pid1
		read.insert(d[20]); // local_dot_prod_p_temp_pid2, CopyJob instance
		write.clear();
		write.insert(d[11]); // global_sum
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
		after.insert(j[9]);
		SpawnCopyJob(j[17], d[11], d[21], before, after, par);

		// Project_Forloop_Part4, pid = 1
		read.clear();
		read.insert(d[4]); // rho
		read.insert(d[11]); // global_sum
		read.insert(da[5]); // x_interior_pid1
		read.insert(d[5]); // p_interior_pid1
		read.insert(da[3]); // b_interior_pid1
		read.insert(d[0]); // tmp_interior_pid1
		write.clear();
		write.insert(da[5]); // x_interior_pid1
		write.insert(da[3]); // b_interior_pid1
		before.clear();
		before.insert(j[7]);
		after.clear();
		after.insert(j[10]);
		SpawnComputeJob("Project_Forloop_Part4", j[8], read, write, before, after, par);

		// Project_Forloop_Part4, pid = 2
		read.clear();
		read.insert(d[17]); // rho, CopyJob instance
		read.insert(d[21]); // global_sum, CopyJob instance
		read.insert(da[6]); // x_interior_pid2
		read.insert(d[6]); // p_interior_pid2
		read.insert(da[4]); // b_interior_pid2
		read.insert(d[1]); // tmp_interior_pid2
		write.clear();
		write.insert(da[6]); // x_interior_pid2
		write.insert(da[4]); // b_interior_pid2
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
		SpawnCopyJob(j[18], da[4], d[22], before, after, par);

		// Global_Max_Abs
		read.clear();
		read.insert(da[3]); // b_interior_pid1
		read.insert(d[22]); // b_interior_pid2, CopyJob instance
		write.clear();
		write.insert(da[0]); // residual
		before.clear();
		before.insert(j[8]);
		before.insert(j[9]);
		before.insert(j[18]); // SpawnCopyJob
		after.clear();
		after.insert(j[11]);
		SpawnComputeJob("Global_Max_Abs", j[10], read, write, before, after, par);

		// Project_Forloop_Condition
		read.clear();
		read.insert(da[0]);
		write.clear();
		before.clear();
		before.insert(j[10]);
		after.clear();
		param_idset.clear();
		for (int i=0;i<NUM_OF_FORLOOP_INPUTS;i++) param_idset.insert(da[i]);
		param_idset.insert(iteration+1);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Condition", j[11], read, write, before, after, par);
	}
	else {
		GetNewJobID(&j, 2);
		READ_0();
		WRITE_0();
		BEFORE_0();
		AFTER_0();
		SpawnComputeJob("finish", j[0], read, write, before, after, par);
		READ_0();
		WRITE_0();
		BEFORE_0();
		AFTER_0();
		SpawnComputeJob("finish", j[1], read, write, before, after, par);		
	}	
};

Project_Forloop_Part1::Project_Forloop_Part1(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part1::Clone() {
	std::cout << "Cloning Project_Forloop_Part1 job!\n";
	return new Project_Forloop_Part1(application());
};

void Project_Forloop_Part1::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Project_Forloop_Part1 job\n";
	
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;	
	
	//execution
	int color = 1;
	SPARSE_MATRIX_FLAT_NXN<T>& A = (*app_driver->projection_data->A_array)(color);
	VECTOR_ND<T>& z_interior = (*app_driver->projection_internal_data->z_interior);
	VECTOR_ND<T>& b_interior = (*app_driver->projection_internal_data->b_interior);
	VECTOR_ND<T>& temp_interior = (*app_driver->projection_internal_data->temp_interior);
	A.C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
	A.C->Solve_Backward_Substitution(temp_interior,z_interior,false,true); // diagonal is inverted to save on divides

	std::cout << "Completed Project_Forloop_Part1 job\n";
};

Project_Forloop_Part2::Project_Forloop_Part2(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part2::Clone() {
	std::cout << "Cloning Project_Forloop_Part2 job!\n";
	return new Project_Forloop_Part2(application());
};

void Project_Forloop_Part2::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Project_Forloop_Part2 job\n";
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	PhysBAM::ProjectionInternalData< PhysBAM::VECTOR<float,2> >* internal = app_driver->projection_internal_data;
	
	// execution
	int iteration = *(params.idset().begin());
	VECTOR_ND<float>& z_interior = (*internal->z_interior);
	VECTOR_ND<float>& b_interior = (*internal->b_interior);
	VECTOR_ND<float>& p_interior = (*internal->p_interior);
	
	internal->rho_old = internal->rho;
	internal->rho = app_driver->pcg_mpi.Global_Sum(VECTOR_ND<float>::Dot_Product_Double_Precision(z_interior, b_interior));
	internal->beta = 0;
	if (iteration == 1) {
		p_interior = z_interior;
	} else {
		internal->beta = internal->rho / internal->rho_old;		
		for(int i=1;i<=internal->interior_n;i++)
			p_interior(i) = z_interior(i) + internal->beta * p_interior(i);
	}
	std::cout << "Completed Project_Forloop_Part2 job\n";
};

Project_Forloop_Part3::Project_Forloop_Part3(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part3::Clone() {
	std::cout << "Cloning Project_Forloop_Part3 job!\n";
	return new Project_Forloop_Part3(application());
};

void Project_Forloop_Part3::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Project_Forloop_Part3 job\n";
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	PhysBAM::ProjectionInternalData< PhysBAM::VECTOR<float,2> >* internal = app_driver->projection_internal_data;
	
	// execution
	app_driver->pcg_mpi.Fill_Ghost_Cells(internal->p);
	int color = 1;
	SPARSE_MATRIX_FLAT_NXN<float>& A = (*projection_data->A_array)(color);
	VECTOR_ND<float>& temp = (*internal->temp);
	VECTOR_ND<float>& p = (*internal->p);
	A.Times(p, temp);
	
	std::cout << "Completed Project_Forloop_Part3 job\n";
};

Project_Forloop_Part4::Project_Forloop_Part4(Application* app) {
	set_application(app);
};

Job * Project_Forloop_Part4::Clone() {
	std::cout << "Cloning Project_Forloop_Part4 job!\n";
	return new Project_Forloop_Part4(application());
};

void Project_Forloop_Part4::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Project_Forloop_Part4 job\n";
	// load driver
	App* projection_app = dynamic_cast<App*>(application());
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	PhysBAM::ProjectionInternalData< PhysBAM::VECTOR<float,2> >* internal = app_driver->projection_internal_data;
	
	// execution
	VECTOR_ND<float>& p_interior = (*internal->p_interior);
	VECTOR_ND<float>& x_interior = (*internal->x_interior);
	VECTOR_ND<float>& b_interior = (*internal->b_interior);
	VECTOR_ND<float>& temp_interior = (*internal->temp_interior);
	
	float alpha=(float)(internal->rho/app_driver->pcg_mpi.Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior)));
	for(int i=1;i<=internal->interior_n;i++) {
		x_interior(i) += alpha * p_interior(i);
		b_interior(i) -= alpha * temp_interior(i);
	}
	std::cout << "Completed the Project_Forloop_Part4 job\n";
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
	std::cout << "Completed Global_Sum job\n";
};

Global_Max_Abs::Global_Max_Abs(Application* app) {
	set_application(app);
};

Job * Global_Max_Abs::Clone() {
	std::cout << "Cloning Global_Max_Abs job!\n";
	return new Global_Max_Abs(application());
};

void Global_Max_Abs::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Global_Max_Abs job\n";
	std::cout << "Completed Global_Max_Abs job\n";
};

Finish::Finish(Application* app) {
	set_application(app);
}

void Finish::Execute(Parameter params, const DataArray& da) {
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Finish job starts on worker %d.\n", projection_app->_rankID);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
	app_driver->pcg_mpi->ExchangePressure(app_driver->projection_internal_data, app_driver->projection_data);
	
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


