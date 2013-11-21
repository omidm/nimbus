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

#include "ProtoBuf_test.h"

using Sparse_Matrix::Sparse_Matrix_Float;

Sparse_Matrix::Sparse_Matrix() {
	//matrix_ = matrix;
};

Sparse_Matrix::~Sparse_Matrix() {
};

Data * Sparse_Matrix::Clone() {
	std::cout << "Cloning Sparse_Matrix data!\n";
	return new Sparse_Matrix(matrix_);
};

void Sparse_Matrix::Create() {
	matrix_ = new SPARSE_MATRIX_FLAT_NXN();
};

void Sparse_Matrix::Destroy() {
	delete matrix_;
};

void Sparse_Matrix::Copy(Data* from) {
	Sparse_Matrix_Float *d = reinterpret_cast<Sparse_Matrix_Float*>(from);
	matrix_->n = d->n;
	matrix_->offsets = new ARRAY<int>(d->offsets.m);
	for (int i = 1; i <= matrix_->offsets.m; i++) 
		matrix_->offsets(i) = d->offsets(i); 
	matrix_->A = new ARRAY<SPARSE_MATRIX_ENTRY<float>>(d->A.m);
	for (int i = 1; i <= matrix_->A.m; i++)
		matrix_->A(i) = d->A(i);
};

bool Sparse_Matrix::Serialize(SerializedData* ser_data) {
	Sparse_Matrix_Float matrix_msg;
	matrix_msg.add_n(matrix_->n);
	matrix_msg.set_allocated_offsets(matrix_->offsets);
	matrix_msg.set_allocated_a(matrix_->a);
	
	std::string str;
	matrix_msg.SerializeToString(&str);
	char* ptr = new char[str.length()];
	memcpy(ptr, str.c_str(), str.length());
	ser_data->set_data_ptr(ptr);
	ser_data->set_size(str.length());
	return true;
};

bool Sparse_Matrix::DeSerialize(const SerializedData& ser_data, Data** result) {
	Sparse_Matrix_Float matrix_msg;
	std::string str(ser_data.data_ptr_raw(), ser_data.size());
	matrix_msg.ParseFromString(str);
	Sparse_Matrix* sparseM = new Sparse_Matrix();
	sparseM->Create();
	sparseM->matrix_->n = matrix_msg.n();
	int _size;
	if (matrix_msg.offsets())
		_size = matrix_msg.offsets()->m();
	else
		_size = 0;
	sparseM->matrix_->offsets = new ARRAY<int>(_size);
	for (int i = 1; i <= _size; i++)
		sparseM->matrix_->offsets(i) = matrix_msg.offsets().elem(i-1);
	
	if (matrix_msg.a())
		_size = matrix_msg.a()->m();
	else
		_size = 0;
	sparseM->matrix_->a = new ARRAY<SPARSE_MATRIX_ENTRY<float>>(_size);
	for (int i = 1; i <= _size; i++)
		sparseM->matrix_->a(i) = matrix_msg.a().elem(i-1);

	*result = sparseM;
	return true;
};

TestApp::TestApp() {
}
;

void TestApp::Load() {
	printf("Worker beginning to load application\n");

	/* Declare and initialize data, jobs and policies. */

	RegisterJob("main", new Main(this));
	RegisterJob("init", new Init(this));
	RegisterJob("verify",	new Project_Forloop_Condition(this));	

	RegisterData("matrix", new Sparse_Matrix());

	printf("Finished creating job and data definitions\n");
}

Main::Main(Application *app) {
	set_application(app);
};

Job* Main::Clone() {
	printf("Cloning Main job\n");
	return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
	printf("Begin Main\n");
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	partition_id_t p_1 = 1;
	partition_id_t p_2 = 2;
	Parameter par;
	IDSet<param_id_t> param_idset;

	GetNewJobID(&j, 3);
	GetNewLogicalDataID(&d, 2);
	
	DefineData("matrix", d[0], p_1, neighbor_partitions, par);
	DefineData("matrix", d[1], p_2, neighbor_partitions, par);
	
	read.clear();
	write.clear();write.insert(d[0]);
	before.clear();
	after.clear();after.insert(j[1]);
	SpawnComputeJob("init", j[0], read, write, before, after, par);
	
	before.clear();before.insert(j[0]);
	after.clear();after.insert(j[2]);
	SpawnCopyJob(j[1], d[0], d[1], before, after, par);
	
	read.clear();read.insert(d[1]);
	write.clear();
	before.clear();before.clear(j[1]);
	after.clear();
	SpawnComputeJob("verify", j[2], read, write, before, after, par);
	printf("Completed Main\n");
};

Initialization::Initialization(Application *app) {
	set_application(app);
};

Job* Initialization::Clone() {
	printf("Cloning Initialization job\n");
	return new Initialization(application());
};

void Initialization::Execute(Parameter params, const DataArray& da) {
	printf("Begin Initialization\n");
	SparseMatrix *d = reinterpret_cast<SparseMatrix*>(da[0]);
	d->matrix_->n = ARRAY_SIZE;
	d->matrix_->offsets = new ARRAY<int>(ARRAY_SIZE);
	for (int i = 1; i <= ARRAY_SIZE; i++)
		d->matrix_->offsets(i) = i;
	d->matrix_->a = new ARRAY<SPARSE_MATRIX_ENTRY<float>>(ARRAY_SIZE);
	for (int i = 1; i <= ARRAY_SIZE; i++)
		d->matrix_->a(i) = new SPARSE_MATRIX_ENTRY(i, i+0.5);		
	printf("Completed Initialization\n");
};

Verification::Verification(Application *app) {
	set_application(app);
};

Job* Verification::Clone() {
	printf("Cloning Verification job\n");
	return new Verification(application());
};

void Verification::Execute(Parameter params, const DataArray& da) {
	printf("Begin Verification\n");
	SparseMatrix *d = reinterpret_cast<SparseMatrix*>(da[0]);
	printf("n = %d\n", d->matrix_->n);
	printf("offsets.m = %d\n", d->matrix_->offsets.m);
	for (int i = 1; i <= d->matrix_->offsets.m; i++)
		printf("offsets(%d) = %d\n", i, d->matrix_->offsets(i));
	
	printf("a.m = %d\n", d->matrix_->a.m);
	for (int i = 1; i <= d->matrix_->a.m)
		printf("a(%d) = (%d, %f)\n", i, d->matrix_->a(i).j, d->matrix_->a(i).a);
	
	printf("Completed Verification\n");
};
