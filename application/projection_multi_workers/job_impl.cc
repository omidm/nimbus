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

#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include "projection_driver.h"
#include "projection_example.h"

#include "app.h"
#include "data_impl.h"
#include "job_impl.h"
#include "PCG_Sparse_Solver.h"

#define NUM_OF_FORLOOP_PARAS 10

Main::Main(Application* app) {
	set_application(app);
}

Job* Main::Clone() {
	printf("Cloning Main job\n");
	return new Main(application());
}

void Main::Execute(Parameter params, const DataArray& da) {	
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Main job starts on worker %d.\n", projection_app->_rankID);
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	IDSet<param_id_t> param_idset;
	partition_id_t pid1 = 1, pid2 = 2;
	Parameter par;	
	GetNewLogicalDataID(&d, 30);
	GetNewJobID(&j, 3);
	DefineData("partial_norm", d[0], pid1, neighbor_partitions, par);
	DefineData("partial_norm", d[1], pid2, neighbor_partitions, par);
	
	// intermidate data
	DefineData("matrix", d[2], pid1, neighbor_partitions, par); // A_pid1
	DefineData("matrix", d[3], pid2, neighbor_partitions, par); // A_pid2
	DefineData("matrix", d[4], pid1, neighbor_partitions, par); // AC_pid1
	DefineData("matrix", d[5], pid2, neighbor_partitions, par); // AC_pid2
	DefineData("vector", d[6], pid1, neighbor_partitions, par); // b_interior_pid1
	DefineData("vector", d[7], pid2, neighbor_partitions, par); // b_interior_pid2
	DefineData("vector", d[8], pid1, neighbor_partitions, par); // z_interior_pid1
	DefineData("vector", d[9], pid2, neighbor_partitions, par); // z_interior_pid2	
	
	// execution
	read.clear();
	read.insert(d[2]);
	write.clear();
	write.insert(d[2]); // A_pid1
	write.insert(d[4]); // AC_pid1
	write.insert(d[6]); // b_interior_pid1
	write.insert(d[8]); // z_interior_pid1
	before.clear();
	after.clear();
	after.insert(j[2]);
	SpawnComputeJob("Init", j[0], read, write, before, after, par);
	read.clear();
	read.insert(d[3]);
	write.clear();
	write.insert(d[3]); // A_pid2
	write.insert(d[5]); // AC_pid2
	write.insert(d[7]); // b_interior_pid2
	write.insert(d[9]); // z_interior_pid2
	before.clear();
	after.clear();
	after.insert(j[2]);
	SpawnComputeJob("Init", j[1], read, write, before, after, par);
	read.clear();
	read.insert(d[2]);
	read.insert(d[3]);
	write.clear();
	before.clear();
	before.insert(j[0]);
	before.insert(j[1]);
	after.clear();
	param_idset.clear();
	for (int i = 0; i < NUM_OF_FORLOOP_PARAS; i++)
		param_idset.insert(d[i]);
	param_idset.insert(1); // insert iteration
	par.set_idset(param_idset);
	SpawnComputeJob("Project_Forloop_Condition", j[2], read, write, before,	after, par);
	dbg(DBG_PROJ, "||Main job finishes on worker %d.\n", projection_app->_rankID);
}
