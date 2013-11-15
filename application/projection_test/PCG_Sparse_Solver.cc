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
#include "projection_app.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

using vector_msg::VectorMsg;

Vec::Vec(int size) {
  size_ = size;
};

Vec::~Vec() {
};

Data * Vec::Clone() {
  std::cout << "Cloning Vec data!\n";
  return new Vec(size_);
};

void Vec::Create() {
  arr_ = new T[size_];
};

void Vec::Destroy() {
  delete arr_;
};

void Vec::Copy(Data* from) {
  Vec *d = reinterpret_cast<Vec*>(from);
  for (int i = 0; i < size_; i++)
    arr_[i] = d->arr()[i];
}

bool Vec::Serialize(SerializedData* ser_data) {
  VectorMsg vec_msg;
  for (int i = 0; i < size_; i++)
    vec_msg.add_elem(arr_[i]);
  std::string str;
  vec_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());
  return true;
}

bool Vec::DeSerialize(const SerializedData& ser_data, Data** result) {
  VectorMsg vec_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  vec_msg.ParseFromString(str);
  Vec* vec = new Vec(size_);
  vec->Create();
  for (int i = 0; (i < size_) && (i < vec_msg.elem_size()); i++)
     vec->arr()[i] = vec_msg.elem(i);

  *result = vec;
  return true;
}

int Vec::size() {
  return size_;
}

T* Vec::arr() {
  return arr_;
}

Init::Init(Application *app) {
	set_application(app);
};

Job* Init::Clone() {
	printf("Cloning Init job\n");
	return new Init(application());
};

void Init::Execute(Parameter params, const DataArray& da) {
	printf("Begin Init\n");
	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // AC
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // b_interior
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // x_interior
	printf("Checkpoint #1\n");
	printf("d0->size() = %d\n", d0->size());
	printf("d1->size() = %d\n", d1->size());
	printf("d2->size() = %d\n", d2->size());
	for (int i=0; i<d0->size(); i++) d0->arr()[i] = i;
	printf("Checkpoint #2\n");
	for (int i=0; i<d1->size(); i++) d1->arr()[i] = 1;
	printf("Checkpoint #3\n");
	for (int i=0; i<d2->size(); i++) d2->arr()[i] = 0;
	printf("Completed Init\n");
};

Project_Forloop_Condition::Project_Forloop_Condition(Application* app) {
	set_application(app);
}
;

Job * Project_Forloop_Condition::Clone() {
	std::cout << "Cloning Project_Forloop_Condition job!\n";
	return new Project_Forloop_Condition(application());
}
;

void Project_Forloop_Condition::Execute(Parameter params, const DataArray& input_data) {
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

	// new data
	GetNewLogicalDataID(&d, 23);
	DefineData("vector", d[0], pid1, neighbor_partitions, par); // temp_interior_pid1
	DefineData("vector", d[1], pid2, neighbor_partitions, par); // temp_interior_pid2
	DefineData("vector", d[2], pid1, neighbor_partitions, par); // local_dot_prod_z_b_pid1
	DefineData("vector", d[3], pid2, neighbor_partitions, par); // local_dot_prod_z_b_pid2
	DefineData("scalar", d[4], pid1, neighbor_partitions, par); // rho
	DefineData("vector", d[5], pid1, neighbor_partitions, par); // p_interior_pid1
	DefineData("vector", d[6], pid2, neighbor_partitions, par); // p_interior_pid2
	DefineData("vector", d[7], pid1, neighbor_partitions, par); // p_boundary_pid1
	DefineData("vector", d[8], pid2, neighbor_partitions, par); // p_boundary_pid2
	DefineData("vector", d[9], pid1, neighbor_partitions, par); // local_dot_prod_p_temp_pid1
	DefineData("vector", d[10], pid2, neighbor_partitions, par); // local_dot_prod_p_temp_pid2
	DefineData("scalar", d[11], pid1, neighbor_partitions, par); // global_sum
	DefineData("scalar", d[12], pid1, neighbor_partitions, par); // rho_old_pid1
	DefineData("scalar", d[13], pid2, neighbor_partitions, par); // rho_old_pid2
	DefineData("vector", d[14], pid1, neighbor_partitions, par); // z_interior_pid1
	DefineData("vector", d[15], pid2, neighbor_partitions, par); // z_interior_pid2
	
	// copy related data
	DefineData("vector", d[16], pid1, neighbor_partitions, par); // local_dot_prod_z_b_pid2 => copy to pid1
	DefineData("vector", d[17], pid2, neighbor_partitions, par); // rho => copy to pid2
	DefineData("vector", d[18], pid1, neighbor_partitions, par); // p_interior_pid2 => copy to pid1
	DefineData("vector", d[19], pid2, neighbor_partitions, par); // p_interior_pid1 => copy to pid2
	DefineData("vector", d[20], pid1, neighbor_partitions, par); // local_dot_prod_p_temp_pid2 => copy to pid1
	DefineData("vector", d[21], pid2, neighbor_partitions, par); // global_sum => copy to pid2
	DefineData("vector", d[22], pid1, neighbor_partitions, par); // b_interior_pid2 => copy to pid1
	
	// executiong part
	int iteration = da[NUM_OF_FORLOOP_INPUTS];
	T residual = d0->arr()[0];
	printf("Jia: forloop before check, iter = %d, res = %f\n", iteration, residual);
	if(iteration == 1 || (iteration < DESIRED_ITERATIONS && residual> GLOBAL_TOLERANCE)) {
		
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
		param_idset.clear(); param_idset.insert(iteration);
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
		param_idset.clear(); param_idset.insert(iteration);
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
		after.insert(j[18]); // SpawnCopyJob
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

	}
	// TODO: Fill_Ghost_Cells(x);
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

	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // AC
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // b_interior
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // z_interior
	Vec *d3 = reinterpret_cast<Vec*>(da[3]); // local_dot_prod_zb

	printf("Begin SPARSE_MATRIX_FLAT_NXN\n");
	SPARSE_MATRIX_FLAT_NXN<T> AC;
	ARRAY<int> rows(d1->size());
	for (int i=1; i<=rows.m; i++) rows(i) = rows.m;
	AC.Set_Row_Lengths(rows);
	int pos = 0;
	for (int i=1; i<=rows.m; i++)
		for(int j=1; j<=rows.m; j++) {
			AC.Set_Element(i,j,d0->arr()[pos++]);
		}
	
	VECTOR_ND<T> b_interior(rows.m, false), temp_interior(rows.m, false), z_interior(rows.m, false);
	for (int i=1; i<=rows.m; i++) b_interior(i) = d1->arr()[i-1];
	
	printf("Finish SPARSE_MATRIX_FLAT_NXN init.\n");
	//AC.Solve_Forward_Substitution(b_interior, temp_interior, true); // diagonal should be treated as the identity
	printf("Checkpoint #4\n");
	//AC.Solve_Backward_Substitution(temp_interior, z_interior, false, true); // diagonal is inverted to save on divides
	printf("Finish SPARSE_MATRIX_FLAT_NXN solve.\n");
	
	for (int i=1; i<=rows.m; i++) z_interior(i) = b_interior(i); // TODO: switch back to incomplete_cholesky in the future
	for (int i=1; i<=rows.m; i++) d2->arr()[i-1] = z_interior(i);
	for (int i=1; i<=rows.m; i++) d3->arr()[i-1] = z_interior(i) * b_interior(i);
	printf("Jia: z(%f, %f, %f, %f)\n", z_interior(1), z_interior(2), z_interior(3), z_interior(4));
	printf("Jia: zb(%f, %f, %f, %f)\n", d3->arr()[0], d3->arr()[1], d3->arr()[2], d3->arr()[3]);
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
	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // rho
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // rho_old
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // z_interior
	Vec *d3 = reinterpret_cast<Vec*>(da[3]); // p_interior
	Vec *d4 = reinterpret_cast<Vec*>(da[4]); // p_interior
	Vec *d5 = reinterpret_cast<Vec*>(da[5]); // rho_old
	printf("Checkpoint #1\n");
	int interation = *(params.idset().begin());
	printf("interation = %d\n", interation);
	T rho = d0->arr()[0], rho_old = d1->arr()[0];
	printf("Jia: rho = %f, rho_old = %f\n", rho, rho_old);
	printf("d2->size() = %d\n", d2->size());
	printf("d4->size() = %d\n", d4->size());
	if (interation == 1 || rho_old == 0) {
		for (int i=0;i<d2->size();i++) d4->arr()[i] = d2->arr()[i];
	}
	else {
		T beta = rho/rho_old;
		for (int i=0;i<d2->size();i++) d4->arr()[i] = d3->arr()[i] * beta + d2->arr()[i];
	}
	d5->arr()[0] = rho;
	printf("Jia: rho_old = %f\n", d5->arr()[0]);
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
	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // A
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // p_interior_pid1
	//Vec *d2 = reinterpret_cast<Vec*>(da[2]); // p_interior_pid2
	Vec *d3 = reinterpret_cast<Vec*>(da[3]); // tmp_interior	
	Vec *d4 = reinterpret_cast<Vec*>(da[4]); // local_dot_prod_p_temp
	printf("Checkpoint #1\n");
	SPARSE_MATRIX_FLAT_NXN<T> AC;
	ARRAY<int> rows(d1->size());
	for (int i=1; i<=rows.m; i++) rows(i) = rows.m;
	AC.Set_Row_Lengths(rows);
	int pos = 0;
	for (int i=1; i<=rows.m; i++)
		for(int j=1; j<=rows.m; j++) {
			AC.Set_Element(i,j,d0->arr()[pos++]);
		}

	VECTOR_ND<T> p_interior(rows.m,false), tmp_interior(rows.m,false);
	for (int i=1;i<=rows.m;i++) p_interior(i) = d1->arr()[i-1];
	printf("Jia: p(%f, %f, %f, %f)\n", p_interior(1), p_interior(2), p_interior(3), p_interior(4));
	AC.Times(p_interior, tmp_interior);
	for (int i=1;i<=rows.m;i++) d3->arr()[i-1] = tmp_interior(i);	
	for (int i=1;i<=rows.m;i++) d4->arr()[i-1] = p_interior(i) * tmp_interior(i);
	printf("Jia: tmp(%f, %f, %f, %f)\n", tmp_interior(1), tmp_interior(2), tmp_interior(3), tmp_interior(4));
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
	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // rho
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // global_sum
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // x_interior
	Vec *d3 = reinterpret_cast<Vec*>(da[3]); // p_interior
	Vec *d4 = reinterpret_cast<Vec*>(da[4]); // b_interior
	Vec *d5 = reinterpret_cast<Vec*>(da[5]); // tmp_interior
	Vec *d6 = reinterpret_cast<Vec*>(da[6]); // x_interior
	Vec *d7 = reinterpret_cast<Vec*>(da[7]); // b_interior

	printf("Checkpoint #1\n");
	printf("d0->arr()[0] = %f, d1->arr()[0] = %f\n", d0->arr()[0], d1->arr()[0]);
	T alpha = d0->arr()[0] / d1->arr()[0];
	for (int i=0;i<d2->size();i++){
		d6->arr()[i] = d2->arr()[i] + alpha * d3->arr()[i];
		d7->arr()[i] = d4->arr()[i] - alpha * d5->arr()[i];
	}
	printf("Jia: x(%f, %f, %f, %f)\n", d6->arr()[0], d6->arr()[1], d6->arr()[2], d6->arr()[3]);
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
	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // part1
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // part2
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // sum
	
	assert(d0);assert(d1);assert(d2);
	printf("Jia: d0->arr(%f, %f, %f, %f)\n", d0->arr()[0], d0->arr()[1], d0->arr()[2], d0->arr()[3]);
	printf("Jia: d1->arr(%f, %f, %f, %f)\n", d1->arr()[0], d1->arr()[1], d1->arr()[2], d1->arr()[3]);
	T sum = 0;
	for (int i=0;i<d0->size();i++) sum += d0->arr()[i];
	for (int i=0;i<d1->size();i++) sum += d1->arr()[i];
	d2->arr()[0] = sum;
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
	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // part1
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // part2
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // max

	T maxmum = 0;
	printf("Jia: d0->arr(%f, %f, %f, %f)\n", d0->arr()[0], d0->arr()[1], d0->arr()[2], d0->arr()[3]);
	printf("Jia: d1->arr(%f, %f, %f, %f)\n", d1->arr()[0], d1->arr()[1], d1->arr()[2], d1->arr()[3]);
	for (int i=0;i<d0->size();i++)
		if (abs(d0->arr()[i]) > maxmum) maxmum = abs(d0->arr()[i]);
	for (int i=0;i<d1->size();i++)
		if (abs(d1->arr()[i]) > maxmum) maxmum = abs(d1->arr()[i]);
	d2->arr()[0] = maxmum;
	printf("Jia: max = %f\n", d2->arr()[0]);
	std::cout << "Completed Global_Max_Abs job\n";
};
