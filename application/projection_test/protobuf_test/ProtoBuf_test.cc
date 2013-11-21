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

Sparse_Matrix::Sparse_Matrix(SPARSE_MATRIX_FLAT_NXN matrix) {
	matrix_ = matrix;
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
	for (int i = 1; i <= matrix->offsets.m; i++) 
		matrix->offsets(i) = d->offsets(i); 
	matrix_->A = new ARRAY<SPARSE_MATRIX_ENTRY<float>>(d->A.m);
	for (int i = 1; i <= matrix_->A.m; i++)
		matrix->A(i) = d->A(i);
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
	Vec* vec = new Vec(size_);
	vec->Create();
	for (int i = 0; (i < size_) && (i < vec_msg.elem_size()); i++)
		vec->arr()[i] = vec_msg.elem(i);

	*result = vec;
	return true;
};
