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

#define LOOP_COUNTER 3
#define LOOP_CONDITION 0
#define STAGE_NUM 10
#define JOB_LENGTH_SEC 0
#define PART_NUM 100
#define CHUNK_NUM 100
#define CHUNK_SIZE 50
#define BANDWIDTH 10
#define STERILE_FLAG true

#define STENCIL_SIZE (2*BANDWIDTH)+1
#define PART_SIZE (CHUNK_NUM/PART_NUM)*CHUNK_SIZE

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
  IDSet<param_id_t> param_idset;



  /*
   * Defining partition and data.
   */
  GetNewLogicalDataID(&d, CHUNK_NUM * 3);

  for (int i = 0; i < CHUNK_NUM; ++i) {
    GeometricRegion r_l(i * CHUNK_SIZE, 0, 0,
                        BANDWIDTH, 1, 1);
    ID<partition_id_t> p_l(i * 3);
    DefinePartition(p_l, r_l, par);
    DefineData(DATA_NAME, d[i * 3], p_l.elem(), neighbor_partitions, par);

    GeometricRegion r_m(i * CHUNK_SIZE + BANDWIDTH, 0, 0,
                        CHUNK_SIZE - 2 * BANDWIDTH, 1, 1);
    ID<partition_id_t> p_m(i * 3 + 1);
    DefinePartition(p_m, r_m, par);
    DefineData(DATA_NAME, d[i * 3 + 1], p_m.elem(), neighbor_partitions, par);

    GeometricRegion r_r(i * CHUNK_SIZE + CHUNK_SIZE - BANDWIDTH, 0, 0,
                        BANDWIDTH, 1, 1);
    ID<partition_id_t> p_r(i * 3 + 2);
    DefinePartition(p_r, r_r, par);
    DefineData(DATA_NAME, d[i * 3 + 2], p_r.elem(), neighbor_partitions, par);
  }

  /*
   * Spawning jobs
   */
  GetNewJobID(&job_ids, CHUNK_NUM * 3 + 1);

  // Spawning inti jobs
  for (int i = 0; i < CHUNK_NUM; ++i) {
    read.clear(); read.insert(d[i * 3]);
    write.clear(); write.insert(d[i * 3]);
    before.clear();
    param_idset.clear(); param_idset.insert(0);
    par.set_idset(param_idset);
    SpawnComputeJob(INIT_JOB_NAME, job_ids[i * 3],
        read, write, before, after, par, STERILE_FLAG);

    read.clear(); read.insert(d[i * 3 + 1]);
    write.clear(); write.insert(d[i * 3 + 1]);
    before.clear();
    param_idset.clear(); param_idset.insert(BANDWIDTH);
    par.set_idset(param_idset);
    SpawnComputeJob(INIT_JOB_NAME, job_ids[i * 3 + 1],
        read, write, before, after, par, STERILE_FLAG);

    read.clear(); read.insert(d[i * 3 + 2]);
    write.clear(); write.insert(d[i * 3 + 2]);
    before.clear();
    param_idset.clear(); param_idset.insert(CHUNK_SIZE - BANDWIDTH);
    par.set_idset(param_idset);
    SpawnComputeJob(INIT_JOB_NAME, job_ids[i * 3 + 2],
        read, write, before, after, par, STERILE_FLAG);
  }

  // Spawning loop job
  read.clear();
  write.clear();
  before.clear();
  for (int j = 0; j < CHUNK_NUM * 3; ++j) {
    before.insert(job_ids[j]);
  }
  after.clear();
  param_idset.clear();
  param_idset.insert(LOOP_COUNTER);
  par.set_idset(param_idset);
  SpawnComputeJob(LOOP_JOB_NAME, job_ids[CHUNK_NUM * 3], read, write, before, after, par);
};

ForLoop::ForLoop(Application* app) {
  set_application(app);
};

Job * ForLoop::Clone() {
  std::cout << "Cloning forLoop job!\n";
  return new ForLoop(application());
};

void ForLoop::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the forLoop job\n";

  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;
  IDSet<param_id_t> param_idset;

  param_id_t loop_counter = *(params.idset().begin());

  if (loop_counter > LOOP_CONDITION) {
    // Spawn the batch of jobs in each stage
    std::vector<job_id_t> stage_job_ids;
    GetNewJobID(&stage_job_ids, STAGE_NUM * PART_NUM);
    std::vector<job_id_t> connector_job_ids;
    GetNewJobID(&connector_job_ids, STAGE_NUM - 1);
    for (int s = 0; s < STAGE_NUM; ++s) {
      for (int i = 0; i < PART_NUM; ++i) {
        read.clear();
        GeometricRegion r_r(i * PART_SIZE - BANDWIDTH, 0, 0,
            PART_SIZE + 2 * BANDWIDTH, 1, 1);
        LoadLogicalIdsInSet(this, &read, r_r, DATA_NAME, NULL);
        write.clear();
        GeometricRegion r_w(i * PART_SIZE, 0, 0,
            PART_SIZE, 1, 1);
        LoadLogicalIdsInSet(this, &write, r_w, DATA_NAME, NULL);
        before.clear();
        if (s > 0) {
          before.insert(connector_job_ids[s - 1]);
        }
        after.clear();
        SpawnComputeJob(STAGE_JOB_NAME, stage_job_ids[s * PART_NUM + i],
            read, write, before, after, par, STERILE_FLAG);
      }
      if (s < (STAGE_NUM - 1)) {
        read.clear();
        write.clear();
        before.clear();
        for (int j = 0; j < PART_NUM; ++j) {
          before.insert(stage_job_ids[s * PART_NUM + j]);
        }
        after.clear();
        SpawnComputeJob(CONNECTOR_JOB_NAME, connector_job_ids[s],
            read, write, before, after, par, STERILE_FLAG);
      }
    }

    // Spawning the print jobs at the end of each loop
    std::vector<job_id_t> print_job_ids;
    GetNewJobID(&print_job_ids, PART_NUM);
    for (int i = 0; i < PART_NUM; ++i) {
      read.clear();
      GeometricRegion r(i * PART_SIZE, 0, 0, PART_SIZE, 1, 1);
      LoadLogicalIdsInSet(this, &read, r, DATA_NAME, NULL);
      write.clear();
      before.clear();
      before.insert(stage_job_ids[(STAGE_NUM - 1) * PART_NUM + i]);
      after.clear();
      SpawnComputeJob(PRINT_JOB_NAME, print_job_ids[i],
          read, write, before, after, par, STERILE_FLAG);
    }

    // Spawning the next for loop job
    std::vector<job_id_t> forloop_job_id;
    GetNewJobID(&forloop_job_id, 1);
    read.clear();
    write.clear();
    before.clear();
    for (int j = 0; j < PART_NUM; ++j) {
      before.insert(print_job_ids[j]);
    }
    after.clear();
    param_idset.clear();
    param_idset.insert(loop_counter - 1);
    par.set_idset(param_idset);
    SpawnComputeJob(LOOP_JOB_NAME, forloop_job_id[0], read, write, before, after, par);
  } else {
    std::vector<job_id_t> print_job_id;
    GetNewJobID(&print_job_id, 1);

    // Spawning final print job
    read.clear();
    GeometricRegion r(0, 0, 0, CHUNK_SIZE * CHUNK_NUM, 1, 1);
    LoadLogicalIdsInSet(this, &read, r, DATA_NAME, NULL);
    write.clear();
    before.clear();
    after.clear();
    SpawnComputeJob(PRINT_JOB_NAME, print_job_id[0],
        read, write, before, after, par, STERILE_FLAG);


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
/*
  std::cout << "Executing the init job\n";
  uint32_t base_val;
  base_val = *(params.idset().begin());
  Vec *d = reinterpret_cast<Vec*>(da[0]);
  for (int i = 0; i < d->size() ; i++)
    d->arr()[i] = base_val + i;
*/
  std::cout << "Executing the init job\n";
  std::vector<int> read_data;
  std::vector<int> write_data;
  LoadDataFromNimbus(this, da, &read_data);

  param_id_t base_val;
  base_val = *(params.idset().begin());
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
/*
  std::cout << "Executing the print job\n";
  std::cout << "OUTPUT: ";
  for (size_t i = 0; i < da.size(); ++i) {
    Vec *d = reinterpret_cast<Vec*>(da[i]);
    for (int j = 0; j < d->size(); ++j)
      std::cout << d->arr()[j] << ", ";
  }
  std::cout << std::endl;
*/
  std::cout << "Executing the print job\n";
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

Stage::Stage(Application* app) {
  set_application(app);
};

Job * Stage::Clone() {
  std::cout << "Cloning stage job!\n";
  return new Stage(application());
};

void Stage::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the stage job\n";

  usleep(1000000 * JOB_LENGTH_SEC);
};

Connector::Connector(Application* app) {
  set_application(app);
};

Job * Connector::Clone() {
  std::cout << "Cloning connector job!\n";
  return new Connector(application());
};

void Connector::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the connector job\n";
  // This job is empty and is meant to just form the job graph.
};





