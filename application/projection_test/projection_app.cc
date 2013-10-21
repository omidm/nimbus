#include "shared/nimbus.h"
#include "projection_app.h"
#define float T

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

ProjectionApp::ProjectionApp() {
}
;

void ProjectionApp::Load() {
	printf("Worker beginning to load application\n");

	/* Declare and initialize data, jobs and policies. */

	RegisterJob("main", new Main(this));

	printf("Finished creating job and data definitions\n");
	printf("Finished loading application\n");
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

void Main::Execute(std::string params, const DataArray& da) {
	printf("Begin main\n");
	std::vector<job_id_t> j;
	std::vector<data_id_t> d;
	IDSet<data_id_t> read, write;
	IDSet<job_id_t> before, after;
	std::string par;

	GetNewJobID(&j, 2);
	before.clear();
	after.clear();
	read.clear();
	write.clear();
	after.insert(j[1]);
	SpawnComputeJob("init", j[0], read, write, before, after, par);
	printf("Spawned init\n");

	before.clear();
	after.clear();
	read.clear();
	write.clear();
	before.insert(j[0]);
	SpawnComputeJob("init", j[0], read, write, before, after, par);

	printf("Completed main\n");
}
;

Project_Forloop_Condition::Project_Forloop_Condition(Application* app) {
	set_application(app);
}
;

Job * Project_Forloop_Condition::Clone() {
	std::cout << "Cloning Project_Forloop_Part1 job!\n";
	return new Project_Forloop_Part1(application());
}
;

void Project_Forloop_Condition::Execute(std::string params, const DataArray& da) {
	std::cout << "Executing the Project_Forloop_Part1 job\n";
	std::vector<job_id_t> j;
	std::vector<data_id_t> d;
	IDSet<data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_t> neighbor_partitions;
	partition_t part1 = 1, part2 = 2, part3 = 3, part4 = 4;
	std::string par;

	char_separator<char> separator("-");
	tokenizer<char_separator<char> > tokens(params, separator);
	tokenizer<char_separator<char> >::iterator iter = tokens.begin();
	int iteration, desired_iterations;
	float global_tolerance;
	nimbus::ParseID(*iter, iteration);
	iter++;
	nimbus::ParseID(*iter, desired_iterations);
	iter++;
	nimbus::ParseID(*iter, global_tolerance);
	
	// input data
	T* residual = (T*) da[0]
	SPARSE_MATRIX_FLAT_NXN<T>* AC = (SPARSE_MATRIX_FLAT_NXN<T>*) da[1];
	VECTOR_ND<T>* b_interior = (VECTOR_ND<T>*) da[2];
	VECTOR_ND<T>* z_interior = (VECTOR_ND<T>*) da[3];
	
	// new data
	GetNewDataID(&d, 4);
	DefineData("temp_interior_part1", d[0], part1, neighbor_partitions, par);
	DefineData("temp_interior_part2", d[1], part2, neighbor_partitions, par);
	DefineData("local_dot_prod_zb_part1", d[2], part1, neighbor_partitions, par);
	DefineData("local_dot_prod_zb_part2", d[3], part2, neighbor_partitions, par);
	
	if(iteration == 1 || iteration < desired_iterations && residual> global_tolerance) {
		// Project_Forloop_Part1, pid = 1
		read.clear();
		read.insert(da[1]); // AC
		read.insert(da[2]); // b_interior
		read.insert(da[3]);  // tmp_interior
		read.insert(da[4]);  // z_interior
		write.clear();
		write.insert(d[1]);
		write.insert(d[2]);
		before.clear(); before.insert(j[0]); before.insert(j[1]);
		after.clear(); after.insert(j[4]);
		SpawnComputeJob("Project_Forloop_Part1", j[2], read, write, before, after, par);
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

void Project_Forloop_Part1::Execute(std::string params, const DataArray& da) {
	std::cout << "Executing the Project_Forloop_Part1 job\n";
	std::vector<job_id_t> j;
	std::vector<data_id_t> d;
	IDSet<data_id_t> read, write;
	IDSet<job_id_t> before, after;
	std::string par;

	AC->Solve_Forward_Substitution(b_interior, temp_interior, true); // diagonal should be treated as the identity
	AC->Solve_Backward_Substitution(temp_interior, z_interior, false, true); // diagonal is inverted to save on divides
}
;

template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::Initialize_Datatypes() {
	MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);
	MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);
	boundary_datatypes.Resize(partition.number_of_sides);
	ghost_datatypes.Resize(partition.number_of_sides);
	for (int s=1; s<=partition.number_of_sides; s++)
		if (partition.neighbor_ranks(s)!=MPI::PROC_NULL) {
			if (partition.boundary_indices(s).m) {
				const ARRAY<int>& displacements=partition.boundary_indices(s);
				ARRAY<int> block_lengths(displacements.m, false);
				ARRAYS_COMPUTATIONS::Fill(block_lengths, 1);
				boundary_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(displacements.m,&block_lengths(1),&displacements(1)); // TODO: collapse consecutive elements into blocks
				boundary_datatypes(s).Commit();
			}
			int ghost_indices_length=partition.ghost_indices(s).Size()+1;
			if (ghost_indices_length) {
				ghost_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(1,&ghost_indices_length,&partition.ghost_indices(s).min_corner);
				ghost_datatypes(s).Commit();
			}
		}
}

