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

#define APP_LOOP_COUNTER 15
#define APP_LOOP_CONDITION 0
#define APP_CHUNK_NUM 10
#define APP_CHUNK_SIZE 5
#define APP_BANDWIDTH 1

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  std::cout << "Cloning main job!\n";
  return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the main job\n";
  assert(APP_CHUNK_SIZE > (2 * APP_BANDWIDTH));

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
  GetNewLogicalDataID(&d, APP_CHUNK_NUM * 3);

  for (int i = 0; i < APP_CHUNK_NUM; ++i) {
    GeometricRegion r_l(i * APP_CHUNK_SIZE, 0, 0,
                        APP_BANDWIDTH, 1, 1);
    ID<partition_id_t> p_l(i * 3);
    DefinePartition(p_l, r_l, par);
    DefineData(DATA_NAME, d[i * 3], p_l.elem(), neighbor_partitions, par);



    GeometricRegion r_m(i * APP_CHUNK_SIZE + APP_BANDWIDTH, 0, 0,
                        APP_CHUNK_SIZE - 2 * APP_BANDWIDTH, 1, 1);
    ID<partition_id_t> p_m(i * 3 + 1);
    DefinePartition(p_m, r_m, par);
    DefineData(DATA_NAME, d[i * 3 + 1], p_m.elem(), neighbor_partitions, par);

    GeometricRegion r_r(i * APP_CHUNK_SIZE + APP_CHUNK_SIZE - APP_BANDWIDTH, 0, 0,
                        APP_BANDWIDTH, 1, 1);
    ID<partition_id_t> p_r(i * 3 + 2);
    DefinePartition(p_r, r_r, par);
    DefineData(DATA_NAME, d[i * 3 + 2], p_r.elem(), neighbor_partitions, par);
  }

  /*
   * Spawning jobs
   */
  GetNewJobID(&job_ids, APP_CHUNK_NUM * 3 + 1);

  for (int i = 0; i < APP_CHUNK_NUM; ++i) {
    read.clear(); read.insert(d[i * 3]);
    write.clear(); write.insert(d[i * 3]);
    before.clear();
    param_idset.clear(); param_idset.insert(0);
    par.set_idset(param_idset);
    SpawnComputeJob(INIT_JOB_NAME, job_ids[i * 3], read, write, before, after, par);

    read.clear(); read.insert(d[i * 3 + 1]);
    write.clear(); write.insert(d[i * 3 + 1]);
    before.clear();
    param_idset.clear(); param_idset.insert(APP_BANDWIDTH);
    par.set_idset(param_idset);
    SpawnComputeJob(INIT_JOB_NAME, job_ids[i * 3 + 1], read, write, before, after, par);

    read.clear(); read.insert(d[i * 3 + 2]);
    write.clear(); write.insert(d[i * 3 + 2]);
    before.clear();
    param_idset.clear(); param_idset.insert(APP_CHUNK_SIZE - APP_BANDWIDTH);
    par.set_idset(param_idset);
    SpawnComputeJob(INIT_JOB_NAME, job_ids[i * 3 + 2], read, write, before, after, par);
  }

  read.clear();
  write.clear();
  before.clear();
  for (int j = 0; j < APP_CHUNK_NUM * 3; ++j) {
    before.insert(job_ids[j]);
  }
  after.clear();
  param_idset.clear();
  param_idset.insert(APP_LOOP_COUNTER);
  par.set_idset(param_idset);
  SpawnComputeJob(LOOP_JOB_NAME, job_ids[APP_CHUNK_NUM * 3], read, write, before, after, par);
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

  std::vector<job_id_t> job_ids;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;
  IDSet<param_id_t> param_idset;

  param_id_t loop_counter = *(params.idset().begin());

  if (loop_counter > APP_LOOP_CONDITION) {
    GetNewJobID(&job_ids, 1);

    read.clear();
    write.clear();
    before.clear();
    after.clear();
    param_idset.clear();
    param_idset.insert(loop_counter - 1);
    par.set_idset(param_idset);
    SpawnComputeJob(LOOP_JOB_NAME, job_ids[0], read, write, before, after, par);
  } else {
    GetNewJobID(&job_ids, APP_CHUNK_NUM);

    for (int i = 0; i < APP_CHUNK_NUM; ++i) {
      read.clear();
      GeometricRegion r(i * APP_CHUNK_SIZE, 0, 0, APP_CHUNK_SIZE, 1, 1);
      LoadLogicalIdsInSet(this, &read, r, DATA_NAME, NULL);
      write.clear();
      before.clear();
      after.clear();
      SpawnComputeJob(PRINT_JOB_NAME, job_ids[i], read, write, before, after, par);
    }

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
  uint32_t base_val;
  base_val = *(params.idset().begin());
  std::cout << "Executing the init job\n";
  Vec *d = reinterpret_cast<Vec*>(da[0]);
  for (int i = 0; i < d->size() ; i++)
    d->arr()[i] = base_val + i;
};


Print::Print() {
};

Job * Print::Clone() {
  std::cout << "Cloning print job!\n";
  return new Print();
};

void Print::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the print job\n";
  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);
  Vec *d3 = reinterpret_cast<Vec*>(da[2]);
  std::cout << "OUTPUT: ";
  std::cout << " d1: ";
  for (int i = 0; i < d1->size(); i++)
    std::cout << d1->arr()[i] << ", ";
  std::cout << " d2: ";
  for (int i = 0; i < d2->size(); i++)
    std::cout << d2->arr()[i] << ", ";
  std::cout << " d3: ";
  for (int i = 0; i < d3->size(); i++)
    std::cout << d3->arr()[i] << ", ";
  std::cout << std::endl;
};




