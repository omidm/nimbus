#include "shared/nimbus.h"
#include "projection_app.h"
#define float T

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
  arr_ = new int[size_];
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

int* Vec::arr() {
  return arr_;
}

ProjectionApp::ProjectionApp() {
};

void ProjectionApp::Load() {
	printf("Worker beginning to load application\n");

	/* Declare and initialize data, jobs and policies. */

	RegisterJob("main", new Main(this));
	RegisterJob("Project_Forloop_Condition", new Project_Forloop_Condition(this));
	RegisterJob("Project_Forloop_Part1", new Project_Forloop_Part1(this));
	RegisterJob("Project_Forloop_Part2", new Project_Forloop_Part2(this));
	RegisterJob("Project_Forloop_Part3", new Project_Forloop_Part3(this));
	RegisterJob("Project_Forloop_Part4", new Project_Forloop_Part4(this));
	RegisterJob("Global_Sum", new Global_Sum(this));
	RegisterJob("Global_Max", new Global_Max(this));

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
	IDSet<partition_t> neighbor_partitions;
	partition_t pid0 = 0, pid1 = 1, pid2 = 2;
	Parameter par;
	IDSet<param_id_t> param_idset;

	GetNewDataID(&d, 7);
	DefineData("scalar", d[0], pid0, neighbor_partitions, par); // residual
	DefineData("matrix", d[1], pid1, neighbor_partitions, par); // A_pid1
	DefineData("matrix", d[2], pid2, neighbor_partitions, par); // A_pid2
	DefineData("vector", d[3], pid1, neighbor_partitions, par); // b_interior_pid1
	DefineData("vector", d[4], pid2, neighbor_partitions, par); // b_interior_pid2
	DefineData("vector", d[5], pid1, neighbor_partitions, par); // x_interior_pid1
	DefineData("vector", d[6], pid2, neighbor_partitions, par); // x_interior_pid2

	GetNewJobID(&j, 1);
	before.clear();
	after.clear();
	read.clear();
	write.clear();
	param_idset.clear();
	for (int i=0;i<7;i++) param_idset.insert(d[i]);
	param_idset.insert(0);
	par.set_idset(param_idset);
	SpawnComputeJob("Project_Forloop_Condition", j[0], read, write, before, after, par);

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

void Project_Forloop_Condition::Execute(Parameter params, const DataArray& input_data) {
	std::cout << "Executing the Project_Forloop_Part1 job\n";
	std::vector<job_id_t> j;
	std::vector<data_id_t> d, da;
	IDSet<data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_t> neighbor_partitions;
	partition_t pid0 = 0, pid1 = 1, pid2 = 2, pid3 = 3, pid4 = 4;
	Parameter par;
	IDSet<param_id_t> param_idset;
	
	// input data
	IDSet<param_id_t>::IDSetContainer::iterator it;
	IDSet<param_id_t> temp_set = params.idset();
	for (it = temp_set.begin(); it != temp_set.end(); it++) {
		da.push_back(*it);
	}
	int iteration = da[7];

	// new data
	GetNewDataID(&d, 16);
	DefineData("vector", d[0], pid1, neighbor_partitions, par); // temp_interior_pid1
	DefineData("vector", d[1], pid2, neighbor_partitions, par); // temp_interior_pid2
	DefineData("vector", d[2], pid1, neighbor_partitions, par); // local_dot_prod_z_b_pid1
	DefineData("vector", d[3], pid2, neighbor_partitions, par); // local_dot_prod_z_b_pid2
	DefineData("scalar", d[4], pid0, neighbor_partitions, par); // rho
	DefineData("vector", d[5], pid1, neighbor_partitions, par); // p_interior_pid1
	DefineData("vector", d[6], pid2, neighbor_partitions, par); // p_interior_pid2
	DefineData("vector", d[7], pid1, neighbor_partitions, par); // p_boundary_pid1
	DefineData("vector", d[8], pid2, neighbor_partitions, par); // p_boundary_pid2
	DefineData("vector", d[9], pid1, neighbor_partitions, par); // local_dot_prod_p_temp_pid1
	DefineData("vector", d[10], pid2, neighbor_partitions, par); // local_dot_prod_p_temp_pid2
	DefineData("scalar", d[11], pid0, neighbor_partitions, par); // global_sum
	DefineData("scalar", d[12], pid1, neighbor_partitions, par); // rho_old_pid1
	DefineData("scalar", d[13], pid2, neighbor_partitions, par); // rho_old_pid2
	DefineData("vector", d[14], pid1, neighbor_partitions, par); // z_interior_pid1
	DefineData("vector", d[15], pid2, neighbor_partitions, par); // z_interior_pid2
	
	if(iteration == 1 || iteration < DESIRED_ITERATIONS && residual> GLOBAL_TOLERANCE) {

		GetNewJobID(&j, 12);

		// Project_Forloop_Part1, pid = 1
		read.clear();
		read.insert(da[1]); // A_pid1
		read.insert(da[3]); // b_interior_pid1
		write.clear();
		write.insert(d[15]); //z_interior_pid1
		write.insert(d[2]); // local_dot_prod_zb_pid1
		before.clear();
		after.clear(); after.insert(j[2]);
		SpawnComputeJob("Project_Forloop_Part1", j[0], read, write, before, after, par);

		// Project_Forloop_Part1, pid = 2
		read.clear();
		read.insert(da[2]); // A_pid2
		read.insert(da[4]); // b_interior_pid2
		write.clear();
		write.insert(d[16]); //z_interior_pid2
		write.insert(d[3]); // local_dot_prod_zb_pid2
		before.clear();
		after.clear(); after.insert(j[2]);
		SpawnComputeJob("Project_Forloop_Part1", j[1], read, write, before, after, par);

		// Global_Sum
		read.clear();
		read.insert(d[2]); // local_dot_prod_zb_pid1
		read.insert(d[3]); // local_dot_prod_zb_pid2
		write.clear();
		write.insert(d[4]); // rho
		before.clear();
		before.insert(j[0]); // Project_Forloop_Part1, pid = 1
		before.insert(j[1]); // Project_Forloop_Part1, pid = 2
		after.clear();
		after.insert(j[3]); // Project_Forloop_Part2, pid = 1
		after.insert(j[4]); // Project_Forloop_Part2, pid = 2
		SpawnComputeJob("Global_Sum", j[2], read, write, before, after, par);

		// Project_Forloop_Part2, pid = 1
		read.clear();
		read.insert(d[4]); // rho
		read.insert(d[12]); // rho_old_pid1
		read.insert(d[15]); // z_interior_pid1
		read.insert(d[5]); // p_interior_pid1
		write.clear();
		write.insert(d[5]); // p_interior_pid1
		before.clear();
		before.insert(j[2]); // Global_Sum
		after.clear();
		after.insert(j[5]); // Project_Forloop_Part3
		param_idset.clear(); param_idset.insert(iteration);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Part2", j[3], read, write, before, after, par);

		// Project_Forloop_Part2, pid = 2
		read.clear();
		read.insert(d[4]); // rho
		read.insert(d[13]); // rho_old_pid2
		read.insert(d[16]); // z_interior_pid2
		read.insert(d[6]); // p_interior_pid2
		write.clear();
		write.insert(d[6]); // p_interior_pid2
		before.clear();
		before.insert(j[2]); // Global_Sum
		after.clear();
		after.insert(j[6]); // Project_Forloop_Part3
		param_idset.clear(); param_idset.insert(iteration);
		par.set_idset(param_idset);
		SpawnComputeJob("Project_Forloop_Part2", j[4], read, write, before, after, par);

		// Project_Forloop_Part3, pid = 1
		read.clear();
		read.insert(da[1]); // A_pid1
		read.insert(d[5]); // p_interior_pid1
		read.insert(d[6]); // p_interior_pid2
		write.clear();
		write.insert(d[0]); // tmp_interior_pid1
		write.insert(d[9]); // local_dot_prod_p_temp_pid1
		before.clear();
		before.insert(j[3]);
		before.insert(j[4]);
		after.clear();
		after.insert(j[7]);
		SpawnComputeJob("Project_Forloop_Part3", j[5], read, write, before, after, par);

		// Project_Forloop_Part3, pid = 2
		read.clear();
		read.insert(da[2]); // A_pid2
		read.insert(d[6]); // p_interior_pid2
		read.insert(d[5]); // p_interior_pid1
		write.clear();
		write.insert(d[1]); // tmp_interior_pid2
		write.insert(d[10]); // local_dot_prod_p_temp_pid2
		before.clear();
		before.insert(j[3]);
		before.insert(j[4]);
		after.clear();
		after.insert(j[7]);
		SpawnComputeJob("Project_Forloop_Part3", j[6], read, write, before, after, par);

		// Global_Sum
		read.clear();
		read.insert(d[9]); // local_dot_prod_p_temp_pid1
		read.insert(d[10]); // local_dot_prod_p_temp_pid2
		write.clear();
		write.insert(d[11]); // global_sum
		before.clear();
		before.insert(j[5]); // Project_Forloop_Part3, pid = 1
		before.insert(j[6]); // Project_Forloop_Part3, pid = 2
		after.clear();
		after.insert(j[8]); // Project_Forloop_Part4, pid = 1
		after.insert(j[9]); // Project_Forloop_Part4, pid = 2
		SpawnComputeJob("Global_Sum", j[7], read, write, before, after, par);

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
		after.clear(j[10]);
		SpawnComputeJob("Project_Forloop_Part4", j[8], read, write, before, after, par);

		// Project_Forloop_Part4, pid = 2
		read.clear();
		read.insert(d[4]); // rho
		read.insert(d[11]); // global_sum
		read.insert(da[6]); // x_interior_pid2
		read.insert(d[6]); // p_interior_pid2
		read.insert(da[4]); // b_interior_pid2
		read.insert(d[1]); // tmp_interior_pid2
		write.clear();
		write.insert(da[6]); // x_interior_pid2
		write.insert(da[4]); // b_interior_pid2
		before.clear();
		before.insert(j[7]);
		after.clear();
		after.clear(j[10]);
		SpawnComputeJob("Project_Forloop_Part4", j[9], read, write, before, after, par);

		// Global_Max
		read.clear();
		read.insert(da[3]); // b_interior_pid1
		read.insert(da[4]); // b_interior_pid2
		write.clear();
		write.insert(da[0]); // residual
		before.clear();
		before.insert(j[8]);
		before.insert(j[9]);
		after.clear();
		after.insert(j[11]);
		SpawnComputeJob("Global_Max", j[10], read, write, before, after, par);

		// Project_Forloop_Condition
		read.clear();
		write.clear();
		before.clear();
		before.insert(j[10]);
		after.clear();
		param_idset.clear();
		for (int i=0;i<7;i++) param_idset.insert(da[i]);
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

	SPARSE_MATRIX_FLAT_NXN<T> AC;
	ARRAY<int> rows;rows.m = d1->size();
	for (int i=1; i<=rows.m; i++) rows(i) = rows.m;
	AC.Set_Row_Lengths(rows);
	int pos = 0;
	for (int i=1; i<=rows.m; i++)
		for(int j=1; j<=rows.m; j++) {
			AC.Set_Element(i,j,d0->arr()[pos++]);
		}
	VECTOR_ND<T> b_interior(rows.m, false), temp_interior(rows.m, false), z_interior(rows.m, false);
	for (int i=1; i<=rows.m; i++) b_interior(i) = d1->arr()[i-1];

	AC.Solve_Forward_Substitution(b_interior, temp_interior, true); // diagonal should be treated as the identity
	AC.Solve_Backward_Substitution(temp_interior, z_interior, false, true); // diagonal is inverted to save on divides
	for (int i=1; i<=rows.m; i++) d2->arr()[i-1] = z_interior(i);
	for (int i=1; i<=rows.m; i++) d3->arr()[i-1] = z_interior(i) * b_interior(i);
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
	int interation = *(params.idset().begin());
	T rho = d0->arr()[0], rho_old = d1->arr()[0];
	if (interation == 1 || rho_old == 0) {
		for (int i=0;i<d2->size();i++) d4->arr()[i] = d2->arr()[i];
	}
	else {
		T beta = rho/rho_old;
		for (int i=0;i<d2->size();i++) d4->arr()[i] = d3->arr()[i] * beta + d2->arr()[i];
	}
	d1->arr()[0] = rho;
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
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // p_interior_pid2
	Vec *d3 = reinterpret_cast<Vec*>(da[3]); // tmp_interior
	Vec *d4 = reinterpret_cast<Vec*>(da[4]); // local_dot_prod_p_temp
	SPARSE_MATRIX_FLAT_NXN<T> AC;
	ARRAY<int> rows;rows.m = d1->size();
	for (int i=1; i<=rows.m; i++) rows(i) = rows.m;
	AC.Set_Row_Lengths(rows);
	int pos = 0;
	for (int i=1; i<=rows.m; i++)
		for(int j=1; j<=rows.m; j++) {
			AC.Set_Element(i,j,d0->arr()[pos++]);
		}

	VECTOR_ND<T> p_interior(rows.m,false), tmp_interior(rows.m,false);
	for (int i=1;i<=rows.m;i++) p_interior(i) = d1->arr()[i-1];
	AC.Times(p_interior, tmp_interior);
	for (int i=1;i<=rows.m;i++) d3->arr()[i-1] = tmp_interior(i);
	for (int i=1;i<=rows.m;i++) d4->arr()[i-1] = p_interior(i) * tmp_interior(i);
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

	T alpha = d0->arr()[0] / d1->arr()[0];
	for (int i=0;i<d2->size();i++){
		d6->arr()[i] = d2->arr()[i] + alpha * d3->arr()[i];
		d7->arr()[i] = d4->arr()[i] - alpha * d5->arr()[i];
	}
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

	T sum = 0;
	for (int i=0;i<d0->size();i++){
		sum += d0->arr()[i] + d1->arr()[i];
	}
	d2->arr()[0] = sum;
};

Global_Max::Global_Max(Application* app) {
	set_application(app);
};

Job * Global_Max::Clone() {
	std::cout << "Cloning Global_Max job!\n";
	return new Global_Max(application());
};

void Global_Sum::Execute(Parameter params, const DataArray& da) {
	std::cout << "Executing the Global_Max job\n";
	Vec *d0 = reinterpret_cast<Vec*>(da[0]); // part1
	Vec *d1 = reinterpret_cast<Vec*>(da[1]); // part2
	Vec *d2 = reinterpret_cast<Vec*>(da[2]); // max

	T maxmum = 0;
	for (int i=0;i<d0->size();i++){
		if (abs(d0->arr()[i]) > maxmum) maxmum = abs(d0->arr()[i]);
		if (abs(d1->arr()[i]) > maxmum) maxmum = abs(d1->arr()[i]);
	}
	d2->arr()[0] = maxmum;
};


/*
 * template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::Initialize_Datatypes() {
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
*/
