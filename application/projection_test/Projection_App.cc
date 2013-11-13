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
	std::vector<data_id_t> d;
	IDSet<data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	partition_id_t pid1 = 1, pid2 = 2;
	Parameter par;
	IDSet<param_id_t> param_idset;

	GetNewDataID(&d, NUM_OF_FORLOOP_INPUTS);
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

	printf("Completed main\n");
};
