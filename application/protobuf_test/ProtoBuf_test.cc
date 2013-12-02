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

using PhysBAM_Protocol::Sparse_Matrix_Float;
using PhysBAM_Protocol::Int_Array;
using PhysBAM_Protocol::Sparse_Matrix_Entry_Float_Array;
using PhysBAM_Protocol::Sparse_Matrix_Entry_Float;

Sparse_Matrix::Sparse_Matrix() {
	//matrix_ = matrix;
};

Sparse_Matrix::~Sparse_Matrix() {
};

Data * Sparse_Matrix::Clone() {
	std::cout << "Cloning Sparse_Matrix data!\n";
	return new Sparse_Matrix();
};

void Sparse_Matrix::Create() {
	matrix_ = new SPARSE_MATRIX_FLAT_NXN<float>();
};

void Sparse_Matrix::Destroy() {
	delete matrix_;
};

void Sparse_Matrix::Copy(Data* from) {
	Sparse_Matrix *d = reinterpret_cast<Sparse_Matrix*>(from);	
	matrix_->n = d->matrix_->n;
	matrix_->offsets = ARRAY<int>(d->matrix_->offsets.m);
	for (int i = 1; i <= matrix_->offsets.m; i++) {
		matrix_->offsets(i) = d->matrix_->offsets(i);
		printf("Checkpoint #3, offsets(%d) = %d\n", i, d->matrix_->offsets(i));
	}
	matrix_->A = ARRAY<SPARSE_MATRIX_ENTRY<float> >(d->matrix_->A.m);
	for (int i = 1; i <= matrix_->A.m; i++) {
		matrix_->A(i).j = d->matrix_->A(i).j;
		matrix_->A(i).a = d->matrix_->A(i).a;
	}
	/*
	printf("Checkpoint #1\n");
	matrix_->n = d->n();
	printf("Checkpoint #2, n = %d\n", d->n());
	matrix_->offsets = ARRAY<int>(d->offsets().m());
	printf("Checkpoint #3\n");
	for (int i = 1; i <= matrix_->offsets.m; i++) 
		matrix_->offsets(i) = d->offsets().elem(i); 
	printf("Checkpoint #4\n");
	matrix_->A = ARRAY<SPARSE_MATRIX_ENTRY<float> >(d->a().m());
	printf("Checkpoint #5\n");
	for (int i = 1; i <= matrix_->A.m; i++) {
		matrix_->A(i).j = d->a().elem(i).j();
		matrix_->A(i).a = d->a().elem(i).a();
	}
	printf("Checkpoint #6\n");
	*/
};

bool Sparse_Matrix::Serialize(SerializedData* ser_data) {
	Sparse_Matrix_Float msg_SparseMatrix;
	msg_SparseMatrix.set_n(matrix_->n);
	Int_Array* msg_IntArray = msg_SparseMatrix.mutable_offsets();
	msg_IntArray->set_m(matrix_->offsets.m);
	for (int i = 1;i <= matrix_->offsets.m; i++) {
		msg_IntArray->add_elem(matrix_->offsets(i));
	}
	Sparse_Matrix_Entry_Float_Array* msg_EntryArray = msg_SparseMatrix.mutable_a();
	msg_EntryArray->set_m(matrix_->A.m);
	for (int i = 1; i <= matrix_->A.m; i++) {
		Sparse_Matrix_Entry_Float* entry = msg_EntryArray->add_elem();		
		entry->set_j(matrix_->A(i).j);
		entry->set_a(matrix_->A(i).a);
	}
	std::string str;
	msg_SparseMatrix.SerializeToString(&str);
	char* ptr = new char[str.length()];
	for (unsigned i = 0; i < str.length(); i++) printf("str[%u] = %u\n", i, (unsigned) str[i]);
	memcpy(ptr, str.c_str(), str.length());
	ser_data->set_data_ptr(ptr);
	ser_data->set_size(str.length());
	Sparse_Matrix_Float msg;
	msg.ParseFromString(str);
	return true;
};

bool Sparse_Matrix::DeSerialize(const SerializedData& ser_data, Data** result) {
	Sparse_Matrix_Float msg_SparseMatrix;
	std::string str(ser_data.data_ptr_raw(), ser_data.size());
	msg_SparseMatrix.ParseFromString(str);
	Sparse_Matrix* sparseM = new Sparse_Matrix();
	sparseM->Create();
	sparseM->matrix_->n = msg_SparseMatrix.n();
	int _size;
	if (msg_SparseMatrix.has_offsets())
		_size = msg_SparseMatrix.offsets().m();
	else
		_size = 0;
	sparseM->matrix_->offsets = ARRAY<int>(_size);
	for (int i = 1; i <= _size; i++)
		sparseM->matrix_->offsets(i) = msg_SparseMatrix.offsets().elem(i-1);
	
	if (msg_SparseMatrix.has_a())
		_size = msg_SparseMatrix.a().m();
	else
		_size = 0;
	sparseM->matrix_->A = ARRAY<SPARSE_MATRIX_ENTRY<float> >(_size);
	for (int i = 1; i <= _size; i++) {
		sparseM->matrix_->A(i).j = msg_SparseMatrix.a().elem(i-1).j();
		sparseM->matrix_->A(i).a = msg_SparseMatrix.a().elem(i-1).a();
	}

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
	RegisterJob("init", new Initialization(this));
	RegisterJob("verify", new Verification(this));	

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
	before.clear();before.insert(j[1]);
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
	Sparse_Matrix *d = reinterpret_cast<Sparse_Matrix*>(da[0]);
	d->matrix_->n = ARRAY_SIZE;
	d->matrix_->offsets = ARRAY<int>(ARRAY_SIZE);
	for (int i = 1; i <= ARRAY_SIZE; i++)
		d->matrix_->offsets(i) = i;
	d->matrix_->A = ARRAY<SPARSE_MATRIX_ENTRY<float> >(ARRAY_SIZE);
	for (int i = 1; i <= ARRAY_SIZE; i++)
		d->matrix_->A(i) = SPARSE_MATRIX_ENTRY<float>(i, i+0.5);		
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
	Sparse_Matrix *d = reinterpret_cast<Sparse_Matrix*>(da[0]);
	printf("n = %d\n", d->matrix_->n);
	printf("offsets.m = %d\n", d->matrix_->offsets.m);
	for (int i = 1; i <= d->matrix_->offsets.m; i++)
		printf("offsets(%d) = %d\n", i, d->matrix_->offsets(i));
	
	printf("a.m = %d\n", d->matrix_->A.m);
	for (int i = 1; i <= d->matrix_->A.m; i++)
		printf("a(%d) = (%d, %f)\n", i, d->matrix_->A(i).j, d->matrix_->A(i).a);
	
	printf("Completed Verification\n");
};
