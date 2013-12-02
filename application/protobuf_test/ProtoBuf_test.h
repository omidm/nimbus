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
#ifndef PROTOBUF_TEST
#define PROTOBUF_TEST

#include "shared/nimbus.h"
#include "shared/parser.h"
#include "shared/nimbus_types.h"
#include "worker/application.h"
#include "worker/job.h"
#include "worker/data.h"
#include "protocol_buffer/Sparse_Matrix_Float.pb.h"
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>

using nimbus::Job;
using nimbus::Data;
using nimbus::Application;
using namespace PhysBAM;
#define ARRAY_SIZE 10

class TestApp : public Application {
public:
	TestApp();
	virtual void Load();
};

class Main : public Job {
public:
	Main(Application *app);
	virtual void Execute(Parameter params, const DataArray& da);
	virtual Job* Clone();
};

class Initialization : public Job {
public:
	Initialization(Application *app);
	virtual void Execute(Parameter params, const DataArray& da);
	virtual Job* Clone();
};

class Verification : public Job {
public:
	Verification(Application *app);
	virtual void Execute(Parameter params, const DataArray& da);
	virtual Job* Clone();
};

class Sparse_Matrix : public Data {
public:
	explicit Sparse_Matrix();
	virtual ~Sparse_Matrix();

	virtual void Create();
	virtual void Destroy();
	virtual Data * Clone();
	virtual void Copy(Data* from);
	virtual bool Serialize(SerializedData* ser_data);
	virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

public:
	SPARSE_MATRIX_FLAT_NXN<float>* matrix_;
};

#endif
