/*
 * Copyright 2013 Stanford University.
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
 * Job classes for job spawner application.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include "./job.h"
#include "./data.h"
#include "./utils.h"
#include "shared/helpers.h"

#define LOOP_COUNTER static_cast<Stencil1DApp*>(application())->counter_
#define LOOP_CONDITION 0
#define PART_NUM static_cast<Stencil1DApp*>(application())->part_num_
#define CHUNK_PER_PART static_cast<Stencil1DApp*>(application())->chunk_per_part_
#define CHUNK_SIZE static_cast<Stencil1DApp*>(application())->chunk_size_
#define BANDWIDTH static_cast<Stencil1DApp*>(application())->bandwidth_

#define STENCIL_SIZE (2*BANDWIDTH)+1
#define PART_SIZE (CHUNK_PER_PART)*CHUNK_SIZE
#define CHUNK_NUM PART_NUM*CHUNK_PER_PART

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  std::cout << "Cloning main job!\n";
  return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the main job\n";
  assert(CHUNK_SIZE > (2 * BANDWIDTH));
  assert(CHUNK_NUM >= PART_NUM);
  assert(CHUNK_NUM % PART_NUM == 0);

  std::vector<job_id_t> job_ids;
  std::vector<logical_data_id_t> d;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;
  Parameter par;



  /*
   * Defining partition and data.
   */
  GetNewLogicalDataID(&d, CHUNK_NUM * 3);

  for (size_t i = 0; i < CHUNK_NUM; ++i) {
    GeometricRegion r_l(i * CHUNK_SIZE, 0, 0,
                        BANDWIDTH, 1, 1);
    ID<partition_id_t> p_l(i * 3);
    DefinePartition(p_l, r_l);
    DefineData(DATA_NAME, d[i * 3], p_l.elem(), neighbor_partitions);

    GeometricRegion r_m(i * CHUNK_SIZE + BANDWIDTH, 0, 0,
                        CHUNK_SIZE - 2 * BANDWIDTH, 1, 1);
    ID<partition_id_t> p_m(i * 3 + 1);
    DefinePartition(p_m, r_m);
    DefineData(DATA_NAME, d[i * 3 + 1], p_m.elem(), neighbor_partitions);

    GeometricRegion r_r(i * CHUNK_SIZE + CHUNK_SIZE - BANDWIDTH, 0, 0,
                        BANDWIDTH, 1, 1);
    ID<partition_id_t> p_r(i * 3 + 2);
    DefinePartition(p_r, r_r);
    DefineData(DATA_NAME, d[i * 3 + 2], p_r.elem(), neighbor_partitions);
  }

  /*
   * Spawning jobs
   */
  GetNewJobID(&job_ids, CHUNK_NUM + 1);

  // StartTemplate("main_job");

  // Spawning inti jobs
  for (size_t i = 0; i < CHUNK_NUM; ++i) {
    GeometricRegion region(i * CHUNK_SIZE, 0, 0, CHUNK_SIZE, 1, 1);
    read.clear();
    LoadLdoIdsInSet(&read, region, DATA_NAME, NULL);
    write.clear();
    LoadLdoIdsInSet(&write, region, DATA_NAME, NULL);
    before.clear();
    StageJobAndLoadBeforeSet(&before, INIT_JOB_NAME, job_ids[i], read, write);
    SerializeParameter(&par, 0);
    SpawnComputeJob(INIT_JOB_NAME, job_ids[i], read, write, before, after, par, true, region);
  }

  MarkEndOfStage();

  // Spawning loop job
  read.clear();
  write.clear();
  before.clear();
  StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, job_ids[CHUNK_NUM], read, write, true);
  SerializeParameter(&par, LOOP_COUNTER);
  SpawnComputeJob(LOOP_JOB_NAME, job_ids[CHUNK_NUM], read, write, before, after, par);

  // EndTemplate("main_job");
};

ForLoop::ForLoop(Application* app) {
  set_application(app);
};

Job * ForLoop::Clone() {
  std::cout << "Cloning forLoop job!\n";
  return new ForLoop(application());
};

void ForLoop::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the forLoop job: " << id().elem() << std::endl;

  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;

  size_t loop_counter;
  LoadParameter(&params, &loop_counter);

  if (loop_counter > LOOP_CONDITION) {
    StartTemplate("for_loop");

    // Spawn the batch of jobs in each stencil
    std::vector<job_id_t> stencil_job_ids;
    GetNewJobID(&stencil_job_ids, PART_NUM);
    for (size_t i = 0; i < PART_NUM; ++i) {
      read.clear();
      GeometricRegion r_r(i * PART_SIZE - BANDWIDTH, 0, 0, PART_SIZE + 2 * BANDWIDTH, 1, 1);
      LoadLdoIdsInSet(&read, r_r, DATA_NAME, NULL);
      write.clear();
      GeometricRegion r_w(i * PART_SIZE, 0, 0, PART_SIZE, 1, 1);
      LoadLdoIdsInSet(&write, r_w, DATA_NAME, NULL);
      before.clear();
      StageJobAndLoadBeforeSet(&before, STENCIL_JOB_NAME, stencil_job_ids[i], read, write);
      after.clear();
      SerializeParameter(&par, i);
      SpawnComputeJob(STENCIL_JOB_NAME, stencil_job_ids[i], read, write, before, after, par, true, r_w); // NOLINT
    }

    MarkEndOfStage();

    // Spawning the print jobs at the end of each loop
    std::vector<job_id_t> print_job_ids;
    GetNewJobID(&print_job_ids, PART_NUM);
    for (size_t i = 0; i < PART_NUM; ++i) {
      read.clear();
      GeometricRegion region(i * PART_SIZE, 0, 0, PART_SIZE, 1, 1);
      LoadLdoIdsInSet(&read, region, DATA_NAME, NULL);
      write.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, PRINT_JOB_NAME, print_job_ids[i], read, write);
      after.clear();
      SpawnComputeJob(PRINT_JOB_NAME, print_job_ids[i], read, write, before, after, par, true, region); // NOLINT
    }

    MarkEndOfStage();

    // Spawning the next for loop job
    std::vector<job_id_t> forloop_job_id;
    GetNewJobID(&forloop_job_id, 1);
    read.clear();
    write.clear();
    before.clear();
    StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, forloop_job_id[0], read, write, true);
    after.clear();
    SerializeParameter(&par, loop_counter - 1);
    SpawnComputeJob(LOOP_JOB_NAME, forloop_job_id[0], read, write, before, after, par);

    EndTemplate("for_loop");
  } else {
    StartTemplate("for_loop_end");

    std::vector<job_id_t> print_job_id;
    GetNewJobID(&print_job_id, 1);

    // Spawning final print job
    read.clear();
    GeometricRegion region(0, 0, 0, CHUNK_SIZE * CHUNK_NUM, 1, 1);
    LoadLdoIdsInSet(&read, region, DATA_NAME, NULL);
    write.clear();
    before.clear();
    after.clear();
    SpawnComputeJob(PRINT_JOB_NAME, print_job_id[0], read, write, before, after, par, true);

    EndTemplate("for_loop_end");

    TerminateApplication();
  }
};

Init::Init() {
};

Job * Init::Clone() {
  std::cout << "Cloning init job!\n";
  return new Init();
};

void Init::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the init job: " << id().elem() << std::endl;
  std::vector<int> read_data;
  std::vector<int> write_data;
  LoadDataFromNimbus(this, da, &read_data);

  size_t base_val;
  LoadParameter(&params, &base_val);
  for (size_t i = 0; i < read_data.size() ; ++i) {
    write_data.push_back(base_val + i);
  }

  SaveDataToNimbus(this, da, &write_data);
};


Print::Print() {
};

Job * Print::Clone() {
  std::cout << "Cloning print job!\n";
  return new Print();
};

void Print::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the print job: " << id().elem() << std::endl;
  std::vector<int> read_data;
  std::vector<int> write_data;
  LoadDataFromNimbus(this, da, &read_data);

  std::cout << "OUTPUT: ";
  for (size_t i = 0; i < read_data.size(); ++i) {
      std::cout << read_data[i] << ", ";
  }
  std::cout << std::endl;

  SaveDataToNimbus(this, da, &write_data);
};

Stencil::Stencil(Application* app) {
  set_application(app);
};

Job * Stencil::Clone() {
  std::cout << "Cloning init job!\n";
  return new Stencil(application());
};

void Stencil::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the stencil job: " << id().elem() << std::endl;
  std::vector<int> read_data;
  std::vector<int> write_data;
  LoadDataFromNimbus(this, da, &read_data);

  for (size_t i = BANDWIDTH; i < (read_data.size() - BANDWIDTH); ++i) {
    int temp = 2 * BANDWIDTH * read_data[i];
    for (size_t j = 1; j <= BANDWIDTH; ++j) {
      temp -= read_data[i - j];
      temp -= read_data[i + j];
    }
    write_data.push_back(temp);
  }

  // compensate for the left-most and right-most partitions.
  size_t part_num;
  LoadParameter(&params, &part_num);
  if (part_num == 0) {
    write_data.insert(write_data.begin(), read_data.begin(), read_data.begin() + BANDWIDTH);
  }
  if (part_num == (PART_NUM - 1)) {
    write_data.insert(write_data.end(), read_data.end() - BANDWIDTH, read_data.end());
  }

  SaveDataToNimbus(this, da, &write_data);
};





