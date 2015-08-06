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

#include <math.h>
#include "./app.h"
#include "./job.h"
#include "./data.h"
#include "./utils.h"
#include "shared/helpers.h"

#define ITERATION_NUM static_cast<LogisticRegression*>(application())->iteration_num()
#define PARTITION_NUM static_cast<LogisticRegression*>(application())->partition_num()
#define PARTITION_SIZE static_cast<LogisticRegression*>(application())->sample_num_per_partition()

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  std::cout << "Cloning main job!\n";
  return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the main job\n";
  assert(PARTITION_NUM > 0);

  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;
  Parameter par;

  /*
   * Defining partition and data.
   */
  std::vector<logical_data_id_t> weight_data_ids;
  std::vector<logical_data_id_t> sample_data_ids;
  GetNewLogicalDataID(&sample_data_ids, PARTITION_NUM);
  GetNewLogicalDataID(&weight_data_ids, PARTITION_NUM);

  for (size_t i = 0; i < PARTITION_NUM; ++i) {
    GeometricRegion r(i * PARTITION_SIZE, 0, 0, PARTITION_SIZE, 1, 1);

    ID<partition_id_t> p(i);
    DefinePartition(p, r);

    DefineData(WEIGHT_DATA_NAME, weight_data_ids[i], p.elem(), neighbor_partitions);
    DefineData(SAMPLE_BATCH_DATA_NAME, sample_data_ids[i], p.elem(), neighbor_partitions);
  }

  /*
   * Spawning the loop job
   */
  std::vector<job_id_t> loop_job_id;
  GetNewJobID(&loop_job_id, 1);

  // StartTemplate("main_job");

  // Spawning loop job
  read.clear();
  write.clear();
  before.clear();
  StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, loop_job_id[0], read, write, true);
  SerializeParameter(&par, ITERATION_NUM);
  SpawnComputeJob(LOOP_JOB_NAME, loop_job_id[0], read, write, before, after, par);

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

  if (loop_counter > 0) {
    StartTemplate("for_loop");

    // Spawn the batch of jobs for gradient stage
    std::vector<job_id_t> gradient_job_ids;
    GetNewJobID(&gradient_job_ids, PARTITION_NUM);
    for (size_t i = 0; i < PARTITION_NUM; ++i) {
      read.clear();
      GeometricRegion r(i * PARTITION_SIZE, 0, 0, PARTITION_SIZE, 1, 1);
      LoadLdoIdsInSet(&read, r, SAMPLE_BATCH_DATA_NAME, NULL);
      write.clear();
      LoadLdoIdsInSet(&write, r, WEIGHT_DATA_NAME, NULL);
      before.clear();
      StageJobAndLoadBeforeSet(&before, GRADIENT_JOB_NAME, gradient_job_ids[i], read, write);
      after.clear();
      SpawnComputeJob(GRADIENT_JOB_NAME, gradient_job_ids[i], read, write, before, after, par, true, r); // NOLINT
    }

    MarkEndOfStage();

    // Spawning the reduction job
    std::vector<job_id_t> reduction_job_id;
    GetNewJobID(&reduction_job_id, 1);
    {
      read.clear();
      GeometricRegion r(0, 0, 0, PARTITION_NUM * PARTITION_SIZE, 1, 1);
      LoadLdoIdsInSet(&read, r, WEIGHT_DATA_NAME, NULL);
      write.clear();
      LoadLdoIdsInSet(&write, r, WEIGHT_DATA_NAME, NULL);
      before.clear();
      StageJobAndLoadBeforeSet(&before, REDUCE_JOB_NAME, reduction_job_id[0], read, write);
      after.clear();
      SpawnComputeJob(REDUCE_JOB_NAME, reduction_job_id[0], read, write, before, after, par, true, r); // NOLINT
    }

    MarkEndOfStage();

    // Spawning the next for loop job
    std::vector<job_id_t> forloop_job_id;
    GetNewJobID(&forloop_job_id, 1);
    {
      read.clear();
      write.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, forloop_job_id[0], read, write, true);
      after.clear();
      SerializeParameter(&par, loop_counter - 1);
      SpawnComputeJob(LOOP_JOB_NAME, forloop_job_id[0], read, write, before, after, par);
    }
    EndTemplate("for_loop");
  } else {
    // StartTemplate("for_loop_end");

    // EndTemplate("for_loop_end");

    TerminateApplication();
  }
};

Gradient::Gradient(Application* app) {
  set_application(app);
};

Job * Gradient::Clone() {
  std::cout << "Cloning gradient job!\n";
  return new Gradient(application());
};

void Gradient::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the gradient job: " << id().elem() << std::endl;
  assert(da.size() == 2);
  SampleBatch *sb = reinterpret_cast<SampleBatch*>(da[0]); // NOLINT
  assert(sb);
  Weight *w = reinterpret_cast<Weight*>(da[1]); // NOLINT
  assert(w);

  std::vector<float> gradient(w->dimension(), 0);
  std::vector<Sample>::iterator iter = sb->samples()->begin();
  for (; iter != sb->samples()->end(); ++iter) {
    float l = iter->label();
    std::vector<float>* x = iter->vector();
    float scale = (1 / (1 + (exp(l * VectorDotProduct(x, w->vector())))) - 1) * l;
    VectorAddWithScale(&gradient, x, scale);
  }

  w->set_vector(gradient);
};

Reduce::Reduce(Application *app) {
  set_application(app);
};

Job * Reduce::Clone() {
  std::cout << "Cloning reduce job!\n";
  return new Reduce(application());
};

void Reduce::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the reduce job: " << id().elem() << std::endl;
  assert(da.size() == 2 * PARTITION_NUM);
  DataArray::const_iterator iter = da.begin();
  for (; iter != da.end(); ++iter) {
    Weight *w = reinterpret_cast<Weight*>(*iter); // NOLINT
    assert(w);
    w->set_vector(std::vector<float>(w->dimension()));
  }
};






