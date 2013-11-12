#include "projection_driver.h"
#include "projection_example.h"
#include "app.h"
#include "job_impl.h"


Main::Main(Application* app) {
  set_application(app);
}

void Main::Execute(Parameter params, const DataArray& da) {
  std::vector<job_id_t> j;
  std::vector<logical_data_id_t> d;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;  // Don't use.
  partition_id_t partition_id1 = 1;
  partition_id_t partition_id2 = 2;
  Parameter par;
  GetNewLogicalDataID(&d, 2);
  GetNewJobID(&j, 2);
  DefineData("region", d[0], partition_id1, neighbor_partitions, par); 
  DefineData("region", d[1], partition_id2, neighbor_partitions, par);
  // Initializes on each partition.
  read.clear(); read.insert(d[0]); write.clear();
  before.clear(); after.clear();
  SpawnComputeJob("initialization", j[0], read, write, before, after, par);
  read.clear(); read.insert(d[1]); write.clear();
  before.clear(); after.clear(); 
  SpawnComputeJob("initialization", j[1], read, write, before, after, par);
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
  App* projection_app = dynamic_cast<App*>(application());
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  std::vector<job_id_t> j;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;  // Don't use.
  Parameter par;
  GetNewJobID(&j, 1);
  job_id_t temp_d = *read_set().begin();
  app_driver->pcg_mpi->ExchangePressure(app_driver->projection_internal_data,
      app_driver->projection_data);
  app_driver->pcg_mpi->InitializeResidual(app_driver->projection_internal_data,
      app_driver->projection_data);
  app_driver->pcg_mpi->SpawnFirstIteration(app_driver->projection_internal_data,
      app_driver->projection_data);
  if (app_driver->projection_internal_data->move_on) {
    app_driver->projection_internal_data->iteration=0;
    read.clear(); read.insert(temp_d); write.clear(); 
    before.clear(); after.clear();
    SpawnComputeJob("spawn_one_iteration_if_needed", j[0], read, write, before, after, par);
  } else {
    read.clear(); read.insert(temp_d); write.clear(); 
    before.clear(); after.clear();
    SpawnComputeJob("finish", j[0], read, write, before, after, par);
  }
}

Job* Initialization::Clone() {
  return new Initialization(application());
}

SpawnOneIterationIfNeeded::SpawnOneIterationIfNeeded(Application* app) {
  set_application(app);
}

void SpawnOneIterationIfNeeded::Execute(Parameter params, const DataArray& da) {
  App* projection_app = dynamic_cast<App*>(application());
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  std::vector<job_id_t> j;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;
  job_id_t temp_d = *read_set().begin();
  // [TODO] This is clearly not the right way to spawn loops.
  if (app_driver->projection_internal_data->move_on) {
    GetNewJobID(&j, 2);
    app_driver->projection_internal_data->iteration++;
    read.clear(); read.insert(temp_d); write.clear(); 
    before.clear(); after.clear(); after.insert(j[1]);
    SpawnComputeJob("one_iteration", j[0], read, write, before, after, par);
    read.clear(); read.insert(temp_d); write.clear(); 
    before.clear(); before.insert(j[0]); after.clear();
    SpawnComputeJob("spawn_one_iteration_if_needed", j[1], read, write, before, after, par);
  } else {
    GetNewJobID(&j, 1);
    read.clear(); read.insert(temp_d); write.clear(); 
    before.clear(); after.clear();
    SpawnComputeJob("finish", j[0], read, write, before, after, par);
  }
}

Job* SpawnOneIterationIfNeeded::Clone() {
  return new SpawnOneIterationIfNeeded(application());
}

OneIteration::OneIteration(Application *app) {
  set_application(app);
}

void OneIteration::Execute(Parameter params, const DataArray& da) {
  App* projection_app = dynamic_cast<App*>(application());
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  app_driver->pcg_mpi->DoPrecondition(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->CalculateBeta(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->UpdateSearchVector(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->ExchangeSearchVector(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->UpdateTempVector(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->CalculateAlpha(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->UpdateOtherVectors(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->CalculateResidual(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->pcg_mpi->DecideToSpawnNextIteration(app_driver->projection_internal_data, app_driver->projection_data);
}

Job* OneIteration::Clone() {
  return new OneIteration(application());
}

Finish::Finish(Application* app) {
  set_application(app);
}
void Finish::Execute(Parameter params, const DataArray& da) {
  App* projection_app = dynamic_cast<App*>(application());
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver = projection_app->app_driver;
  app_driver->pcg_mpi->ExchangePressure(app_driver->projection_internal_data, app_driver->projection_data);
  app_driver->WindUpForOneRegion();
  app_driver->ApplyPressureAndFinish();
  projection_app->FinishMain();
  printf("All is over!!! Nimbus can do projection!!!\n");
}

Job* Finish::Clone() {
  return new Finish(application());
}
