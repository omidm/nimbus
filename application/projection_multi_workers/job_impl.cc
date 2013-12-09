/*
 * The job specification of PhysBAM projection processing.
 * The job specification for projection inner loop
 * is in "job_impl_internal.cc".
 * MPI is heavily used in this version.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include "projection_driver.h"
#include "projection_example.h"

#include "app.h"
#include "data_impl.h"
#include "job_impl.h"
#include "PCG_Sparse_Solver.h"

// TODO The idset implementation should offer the helper function.
#define READ_0() read.clear()
#define READ_1(x) read.clear(); read.insert(x)
#define READ_2(x,y) READ_1(x); read.insert(y)
#define READ_3(x,y,z) READ_2(x,y); read.insert(z)
#define READ_4(x,y,z,e) READ_3(x,y,z); read.insert(e)
#define READ_5(x,y,z,e,f) READ_4(x,y,z,e); read.insert(f)
#define WRITE_0() write.clear()
#define WRITE_1(x) write.clear(); write.insert(x)
#define WRITE_2(x,y) WRITE_1(x); write.insert(y)
#define WRITE_3(x,y,z) WRITE_2(x,y); write.insert(z)
#define BEFORE_0() before.clear();
#define BEFORE_1(x) before.clear(); before.insert(x)
#define BEFORE_2(x,y) BEFORE_1(x); before.insert(y)
#define BEFORE_3(x,y,z) BEFORE_2(x,y); before.insert(z)
#define AFTER_0() after.clear();
#define AFTER_1(x) after.clear(); after.insert(x)
#define AFTER_2(x,y) AFTER_1(x); after.insert(y)
#define AFTER_3(x,y,z) AFTER_2(x,y); after.insert(z)

Main::Main(Application* app) {
	set_application(app);
}

Job* Main::Clone() {
	printf("Cloning Main job\n");
	return new Main(application());
}

void Main::Execute(Parameter params, const DataArray& da) {	
	
}
/*

Main::Main(Application* app) {
	set_application(app);
}

void Main::Execute(Parameter params, const DataArray& da) {
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Main job starts on worker %d.\n", projection_app->_rankID);
	std::vector<job_id_t> j;
	std::vector<logical_data_id_t> d;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	partition_id_t partition_id1 = 1;
	partition_id_t partition_id2 = 2;
	Parameter par;
	GetNewLogicalDataID(&d, 5);
	GetNewJobID(&j, 5);
	DefineData("region", d[0], partition_id1, neighbor_partitions, par);
	DefineData("region", d[1], partition_id2, neighbor_partitions, par);
	DefineData("partial_norm", d[2], partition_id1, neighbor_partitions, par);
	DefineData("partial_norm", d[3], partition_id2, neighbor_partitions, par);
	DefineData("partial_norm", d[4], partition_id1, neighbor_partitions, par);
	// Initializes on each partition,
	// and calculates the partial norm on each part.
	READ_2(d[0], d[2]);
	WRITE_1(d[2]);
	BEFORE_0();
	AFTER_1(j[3]);
	SpawnComputeJob("initialization", j[0], read, write, before, after, par);
	READ_2(d[1], d[3]);
	WRITE_1(d[3]);
	BEFORE_0();
	AFTER_1(j[2]);
	SpawnComputeJob("initialization", j[1], read, write, before, after, par);
	BEFORE_1(j[1]);
	AFTER_1(j[3]);
	SpawnCopyJob(j[2], d[3], d[4], before, after, par);
	READ_5(d[0], d[1], d[2], d[3], d[4]);
	WRITE_0();
	BEFORE_2(j[0], j[2]);
	AFTER_0();
	IDSet<logical_data_id_t> temp;
	temp.insert(d[0]);
	temp.insert(d[1]);
	temp.insert(d[2]);
	temp.insert(d[3]);
	temp.insert(d[4]);
	par.set_idset(temp);
	// TODO(quhang) Includes all the data object in the read set, not necessary.
	SpawnComputeJob("spawn_one_iteration_if_needed", j[3], read, write, before,	after, par);
	dbg(DBG_PROJ, "||Main job finishes on worker %d.\n", projection_app->_rankID);
}

Job* Main::Clone() {
	return new Main(application());
}

Initialization::Initialization(Application* app) {
	set_application(app);
}

void Initialization::Execute(Parameter params, const DataArray& da) {
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Init job starts on worker %d.\n", projection_app->_rankID);
	IDSet<logical_data_id_t> temp = params.idset();
	PartialNorm *partial_norm = reinterpret_cast<PartialNorm*>(da[1]);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver =
			projection_app->app_driver;
	std::vector<job_id_t> j;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	IDSet<partition_id_t> neighbor_partitions;
	Parameter par;
	printf("Jia: Check with new design.\n");
	app_driver->pcg_mpi->Initialize(app_driver->projection_internal_data, app_driver->projection_data);
	app_driver->pcg_mpi->CommunicateConfig(app_driver->projection_internal_data, app_driver->projection_data);
	app_driver->pcg_mpi->ExchangePressure(app_driver->projection_internal_data,	app_driver->projection_data);
	app_driver->pcg_mpi->InitializeResidual(app_driver->projection_internal_data, app_driver->projection_data);
	app_driver->pcg_mpi->SpawnFirstIteration(app_driver->projection_internal_data, app_driver->projection_data);
	app_driver->projection_internal_data->iteration=0;
	partial_norm->norm_ = app_driver->projection_internal_data->partial_norm;
	dbg(DBG_PROJ, "||Init job finishes on worker %d.\n", projection_app->_rankID);
}

Job* Initialization::Clone() {
	return new Initialization(application());
}

SpawnOneIterationIfNeeded::SpawnOneIterationIfNeeded(Application* app) {
	set_application(app);
}

double max(double a, double b) {
	return a>b ? a : b;
}

void SpawnOneIterationIfNeeded::Execute(Parameter params, const DataArray& da) {
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Spawning job starts on worker %d.\n", projection_app->_rankID);
	IDSet<logical_data_id_t> temp_par = params.idset();
	std::vector<logical_data_id_t> temp;
	for (IDSet<logical_data_id_t>::IDSetIter i = temp_par.begin(); i
			!= temp_par.end(); i++) {
		temp.push_back(*i);
	}
	//ProfileData *profile_data_1 = reinterpret_cast<ProfileData*>(da[0]);
	//ProfileData *profile_data_2 = reinterpret_cast<ProfileData*>(da[1]);
	PartialNorm *partial_norm_1 = reinterpret_cast<PartialNorm*>(da[2]);
	//PartialNorm *partial_norm_2 = reinterpret_cast<PartialNorm*>(da[3]);
	PartialNorm *partial_norm_buffer = reinterpret_cast<PartialNorm*>(da[4]);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver =
			projection_app->app_driver;
	std::vector<job_id_t> j;
	IDSet<logical_data_id_t> read, write;
	IDSet<job_id_t> before, after;
	Parameter par;
	bool move_on = true;
	{
		dbg(DBG_PROJ, "||norm1=%f norm2=%f tolerance=%f.\n",
				partial_norm_1->norm_, partial_norm_buffer->norm_,
				app_driver->projection_internal_data->global_tolerance);
		int desired_iterations=app_driver->projection_internal_data->global_n;
		app_driver->projection_internal_data->residual = max(
				partial_norm_1->norm_, partial_norm_buffer->norm_);
		if (app_driver->projection_internal_data->residual
				<= app_driver->projection_internal_data->global_tolerance) {
			move_on = false;
		}
		if (app_driver->projection_internal_data->iteration==desired_iterations) {
			move_on = false;
		}
	}
	if (move_on) {
		GetNewJobID(&j, 4);
		READ_2(temp[0], temp[2]);
		WRITE_1(temp[2]);
		BEFORE_0();
		AFTER_1(j[3]);
		SpawnComputeJob("one_iteration", j[0], read, write, before, after, par);
		READ_2(temp[1], temp[3]);
		WRITE_1(temp[3]);
		BEFORE_0();
		AFTER_1(j[2]);
		SpawnComputeJob("one_iteration", j[1], read, write, before, after, par);
		BEFORE_1(j[1]);
		AFTER_1(j[3]);
		SpawnCopyJob(j[2], temp[3], temp[4], before, after, par);
		READ_5(temp[0], temp[1], temp[2], temp[3], temp[4]);
		WRITE_0();
		BEFORE_2(j[0], j[2]);
		AFTER_0();
		// TODO Add a method to clean the parameter.
		par.set_idset(temp_par);
		SpawnComputeJob("spawn_one_iteration_if_needed", j[3], read, write, before, after, par);
	} else {
		GetNewJobID(&j, 2);
		READ_1(temp[0]);
		WRITE_0();
		BEFORE_0();
		AFTER_0();
		SpawnComputeJob("finish", j[0], read, write, before, after, par);
		READ_1(temp[1]);
		WRITE_0();
		BEFORE_0();
		AFTER_0();
		SpawnComputeJob("finish", j[1], read, write, before, after, par);
	}
	dbg(DBG_PROJ, "||Spawning job finishes on worker %d.\n", projection_app->_rankID);
}

Job* SpawnOneIterationIfNeeded::Clone() {
	return new SpawnOneIterationIfNeeded(application());
}

Finish::Finish(Application* app) {
	set_application(app);
}
void Finish::Execute(Parameter params, const DataArray& da) {
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Finish job starts on worker %d.\n", projection_app->_rankID);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver =
			projection_app->app_driver;
	app_driver->pcg_mpi->ExchangePressure(app_driver->projection_internal_data,
			app_driver->projection_data);
	
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


*/