#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include "projection_driver.h"
#include "projection_example.h"
#include "app.h"
#include "data_impl.h"
#include "job_impl.h"

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

void Main::Execute(Parameter params, const DataArray& da) {
  App* projection_app = dynamic_cast<App*>(application());
  printf("[QUH] Main Job starts.%d\n", projection_app->app_mpi_world->rank);
  std::vector<job_id_t> j;
  std::vector<logical_data_id_t> d;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;  // Don't use.
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
  READ_2(d[0],d[2]); WRITE_1(d[2]);
  BEFORE_0(); AFTER_1(j[3]);
  SpawnComputeJob("initialization", j[0], read, write, before, after, par);
  READ_2(d[1],d[3]); WRITE_1(d[3]);
  BEFORE_0(); AFTER_1(j[2]);
  SpawnComputeJob("initialization", j[1], read, write, before, after, par);
  // Spawn copy job.
  BEFORE_1(j[1]); AFTER_1(j[3]);
  SpawnCopyJob(j[2], d[3], d[4], before, after, par);
  // Spawn one itertion if needed.
  READ_5(d[0],d[1],d[2],d[4],d[3]); WRITE_0();
  BEFORE_2(j[0],j[2]); AFTER_0();
  SpawnComputeJob("spawn_one_iteration_if_needed", j[3], read, write, before, after, par);
  printf("[QUH] Main Job finishes.\n");
}

Job* Main::Clone() {
  return new Main(application());
}

Initialization::Initialization(Application* app) {
  set_application(app);
}

// Each initialization job will spawn corresponding jobs on that partition.
// [TODO] Change it.
void Initialization::Execute(Parameter params, const DataArray& da) {
  PartialNorm *partial_norm = reinterpret_cast<PartialNorm*>(da[1]);
  App* projection_app = dynamic_cast<App*>(application());
  printf("[QUH] Init Job starts.%d\n",projection_app->app_mpi_world->rank);
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  std::vector<job_id_t> j;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;  // Don't use.
  Parameter par;
  // GetNewJobID(&j, 1);
  // logical_data_id_t temp_d = *read_set().begin();
  app_driver->pcg_mpi->ExchangePressure(app_driver->projection_internal_data,
      app_driver->projection_data);
  app_driver->pcg_mpi->InitializeResidual(app_driver->projection_internal_data,
      app_driver->projection_data);
  app_driver->pcg_mpi->SpawnFirstIteration(app_driver->projection_internal_data,
      app_driver->projection_data);
  app_driver->projection_internal_data->iteration=0;
  partial_norm->norm_ = app_driver->projection_internal_data->partial_norm; 
  printf("[QUH] Init Job finishes.\n");

  /*
  if (app_driver->projection_internal_data->move_on) {
    app_driver->projection_internal_data->iteration=0;
    read.clear(); read.insert(temp_d); write.clear(); 
    before.clear(); after.clear();
    SpawnComputeJob("spawn_one_iteration_if_needed", j[0], read, write, before, after, par);
  } else {
    read.clear(); read.insert(temp_d); write.clear(); 
    before.clear(); after.clear();
    SpawnComputeJob("finish", j[0], read, write, before, after, par);
  }*/
}

Job* Initialization::Clone() {
  return new Initialization(application());
}

SpawnOneIterationIfNeeded::SpawnOneIterationIfNeeded(Application* app) {
  set_application(app);
}

double max(double a, double b) { return a>b?a:b;}
void SpawnOneIterationIfNeeded::Execute(Parameter params, const DataArray& da) {
  ProfileData *profile_data_1 = reinterpret_cast<ProfileData*>(da[0]);
  ProfileData *profile_data_2 = reinterpret_cast<ProfileData*>(da[1]);
  PartialNorm *partial_norm_1 = reinterpret_cast<PartialNorm*>(da[2]);
  PartialNorm *partial_norm_2 = reinterpret_cast<PartialNorm*>(da[3]);
  PartialNorm *partial_norm_buffer = reinterpret_cast<PartialNorm*>(da[4]);
  App* projection_app = dynamic_cast<App*>(application());
  printf("[QUH] SpawOneIteration starts.%d\n",projection_app->app_mpi_world->rank);
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  std::vector<job_id_t> j;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;
  printf("[QUH] Before calculate norm.\n");
  bool move_on = true;
  {
    int desired_iterations=app_driver->projection_internal_data->global_n;
    app_driver->projection_internal_data->residual = max(partial_norm_1->norm_, partial_norm_2->norm_);
    if(app_driver->projection_internal_data->residual<=
      app_driver->projection_internal_data->global_tolerance) {
      app_driver->projection_internal_data->move_on = false;
    }
    if(app_driver->projection_internal_data->iteration==desired_iterations) {
      app_driver->projection_internal_data->move_on = false;
    }
  }
  printf("[QUH] After calculate norm.\n");
  if (move_on) {
    GetNewJobID(&j, 4);
    app_driver->projection_internal_data->iteration++;
    // Spawn iteration on each worker.
    READ_1(profile_data_1->logical_id()); WRITE_1(partial_norm_1->logical_id());
    BEFORE_0(); AFTER_1(j[3]);
    SpawnComputeJob("one_iteration", j[0], read, write, before, after, par);
    READ_1(profile_data_2->logical_id()); WRITE_1(partial_norm_buffer->logical_id());
    BEFORE_0(); AFTER_1(j[2]);
    SpawnComputeJob("one_iteration", j[1], read, write, before, after, par);
    // Spawn copy.
    BEFORE_0(); AFTER_1(j[3]);
    SpawnCopyJob(j[2], partial_norm_buffer->logical_id(), partial_norm_2->logical_id(), before, after, par);

    READ_5(
	profile_data_1->logical_id(),
	profile_data_2->logical_id(),
        partial_norm_1->logical_id(),
        partial_norm_2->logical_id(),
        partial_norm_buffer->logical_id());
    WRITE_0();
    BEFORE_2(j[0],j[2]); AFTER_0();
    SpawnComputeJob("spawn_one_iteration_if_needed", j[3], read, write, before, after, par);
  } else {
    GetNewJobID(&j, 2);
    // Do finishing work on each worker.
    read.clear(); read.insert(profile_data_1->logical_id()); write.clear(); 
    before.clear(); after.clear();
    SpawnComputeJob("finish", j[0], read, write, before, after, par);
    read.clear(); read.insert(profile_data_2->logical_id()); write.clear(); 
    before.clear(); after.clear();
    SpawnComputeJob("finish", j[1], read, write, before, after, par);
  }
  printf("[QUH] SpawOneIteration finishes.\n");
}

Job* SpawnOneIterationIfNeeded::Clone() {
  return new SpawnOneIterationIfNeeded(application());
}

OneIteration::OneIteration(Application *app) {
  set_application(app);
}

void OneIteration::Execute(Parameter params, const DataArray& da) {
  PartialNorm *partial_norm = reinterpret_cast<PartialNorm*>(da[1]);
  App* projection_app = dynamic_cast<App*>(application());
  printf("[QUH] OneIteration starts.%d\n",projection_app->app_mpi_world->rank);
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  app_driver->pcg_mpi->DoPrecondition(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->CalculateBeta(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->UpdateSearchVector(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->ExchangeSearchVector(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->UpdateTempVector(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->CalculateAlpha(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->UpdateOtherVectors(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->CalculateResidual(app_driver->projection_internal_data, app_driver->projection_data);
  // Data is calculated and stored in projection_internal_data.
  partial_norm->norm_ = app_driver->projection_internal_data->partial_norm;
  //app_driver->pcg_mpi->DecideToSpawnNextIteration(app_driver->projection_internal_data,
  //    app_driver->projection_data);
  printf("[QUH] OneIteration finishes.\n");
}

Job* OneIteration::Clone() {
  return new OneIteration(application());
}

Finish::Finish(Application* app) {
  set_application(app);
}
void Finish::Execute(Parameter params, const DataArray& da) {
  printf("[QUH] Finish Job starts.\n");
  App* projection_app = dynamic_cast<App*>(application());
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  app_driver->pcg_mpi->ExchangePressure(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->WindUpForOneRegion();
  app_driver->ApplyPressureAndFinish();
  projection_app->FinishMain();
  printf("All is over!!! Nimbus can do projection!!!\n");
  printf("[QUH] Finish Job finishes.\n");
}

Job* Finish::Clone() {
  return new Finish(application());
}
