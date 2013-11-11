#include "projection_driver.h"
#include "projection_example.h"
#include "init_main.h"
#include "app.h"

typedef float T;

/*
  // The workflow for projection.
  PROJECTION_DRIVER<TV> &driver;
  driver.PrepareForProjection();
  driver.PrepareForOneRegion();

  driver.pcg_mpi->ExchangePressure(driver.projection_internal_data, driver.projection_data);
  driver.pcg_mpi->InitializeResidual(driver.projection_internal_data, driver.projection_data);
  driver.pcg_mpi->SpawnFirstIteration(driver.projection_internal_data, driver.projection_data);
  if (driver.projection_internal_data->move_on) {
    driver.projection_internal_data->iteration = 0;
    do {
      driver.projection_internal_data->iteration++;
      driver.pcg_mpi->DoPrecondition(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->CalculateBeta(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->UpdateSearchVector(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->ExchangeSearchVector(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->UpdateTempVector(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->CalculateAlpha(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->UpdateOtherVectors(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->CalculateResidual(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->DecideToSpawnNextIteration(driver.projection_internal_data, driver.projection_data);
    } while (driver.projection_internal_data->move_on);
    driver.pcg_mpi->ExchangePressure(driver.projection_internal_data, driver.projection_data);
  }
  driver.WindUpForOneRegion();
  driver.ApplyPressureAndFinish();
*/

// "driver" is a global pointer to PROJECTION_DRIVER.
// [TODO] Figure out a better way.
using PhysBAM::driver;

void App::Load() {
  std::cout << "Start application loading!" << std::endl;
  RegisterJob("main", new Main(this));
  RegisterJob("initialization", new Initialization(this));
  RegisterJob("spawn_one_iteration_if_needed", new SpawnOneIterationIfNeeded(this));
  RegisterJob("one_iteration", new OneIteration(this));
  RegisterJob("finish", new Finish(this));
  // [TODO] Only dumb data for now.
  RegisterData("region", new ProfileData);

  // Simulates parameter passing mechanism.
  // [TODO] Revise or delete.
  char a[] = "multiple_projection"; 
  char* temp[2];
  temp[0] = a;
  temp[1] = NULL;
  InitMain(1, temp);

  driver->PrepareForProjection();
  driver->PrepareForOneRegion();
  std::cout << "Finish application loading!" << std::endl;
}

Main::Main(Application* app) {
  set_application(app);
}

void Main::Execute(Parameter params, const DataArray& da) {
  std::vector<job_id_t> j;
  std::vector<data_id_t> d;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;  // Don't use.
  partition_id_t partition_id1 = 1;
  partition_id_t partition_id2 = 2;
  Parameter par;
  GetNewDataID(&d, 2);
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
  std::vector<job_id_t> j;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;  // Don't use.
  Parameter par;
  GetNewJobID(&j, 1);
  job_id_t temp_d = *read_set().begin();
  driver->pcg_mpi->ExchangePressure(driver->projection_internal_data, driver->projection_data);
  driver->pcg_mpi->InitializeResidual(driver->projection_internal_data, driver->projection_data);
  driver->pcg_mpi->SpawnFirstIteration(driver->projection_internal_data, driver->projection_data);
  if (driver->projection_internal_data->move_on) {
    driver->projection_internal_data->iteration=0;
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
  std::vector<job_id_t> j;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;
  job_id_t temp_d = *read_set().begin();
  // [TODO] This is clearly not the right way to spawn loops.
  if (driver->projection_internal_data->move_on) {
    GetNewJobID(&j, 2);
    driver->projection_internal_data->iteration++;
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
    driver->pcg_mpi->DoPrecondition(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->CalculateBeta(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->UpdateSearchVector(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->ExchangeSearchVector(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->UpdateTempVector(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->CalculateAlpha(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->UpdateOtherVectors(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->CalculateResidual(driver->projection_internal_data, driver->projection_data);
    driver->pcg_mpi->DecideToSpawnNextIteration(driver->projection_internal_data, driver->projection_data);
}

Job* OneIteration::Clone() {
  return new OneIteration(application());
}

Finish::Finish(Application* app) {
  set_application(app);
}
void Finish::Execute(Parameter params, const DataArray& da) {
  driver->pcg_mpi->ExchangePressure(driver->projection_internal_data, driver->projection_data);
  driver->WindUpForOneRegion();
  driver->ApplyPressureAndFinish();
  FinishMain();
  printf("All is over!!! Nimbus can do projection!!!\n");
}

Job* Finish::Clone() {
  return new Finish(application());
}
