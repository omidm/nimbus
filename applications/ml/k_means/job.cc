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

#include <inttypes.h>
#include <math.h>
#include "./app.h"
#include "./job.h"
#include "./data.h"
#include "./utils.h"
#include "src/shared/helpers.h"

#define CLUSTER_NUM static_cast<KMeans*>(application())->cluster_num()
#define ITERATION_NUM static_cast<KMeans*>(application())->iteration_num()
#define PARTITION_NUM static_cast<KMeans*>(application())->partition_num()
#define SPIN_WAIT_US static_cast<KMeans*>(application())->spin_wait_us()
#define PARTITION_SIZE 1  // to fix the problem with the limit of int_dimension_t
// #define PARTITION_SIZE static_cast<KMeans*>(application())->sample_num_per_partition() // NOLINT
#define REDUCTION_PARTITION_NUM static_cast<KMeans*>(application())->reduction_partition_num() // NOLINT
#define AUTOMATIC_REDUCTION_ACTIVE static_cast<KMeans*>(application())->automatic_reduction_active() // NOLINT
#define REDUCTION_COMBINER_ACTIVE static_cast<KMeans*>(application())->reduction_combiner_active() // NOLINT

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  dbg(DBG_APP, "Cloning main job!\n");
  return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the main job!\n");
  assert(PARTITION_NUM > 0);

  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;
  Parameter par;

  /*
   * Defining partition and data.
   */
  std::vector<logical_data_id_t> sample_data_ids;
  GetNewLogicalDataID(&sample_data_ids, PARTITION_NUM);

  for (size_t i = 0; i < PARTITION_NUM; ++i) {
    GeometricRegion r(i * PARTITION_SIZE, 0, 0, PARTITION_SIZE, 1, 1);

    ID<partition_id_t> p(i);
    DefinePartition(p, r);

    DefineData(SAMPLE_BATCH_DATA_NAME, sample_data_ids[i], p.elem(), neighbor_partitions);
  }

  // differentiate the automatic and manual reduction cases!
  if (AUTOMATIC_REDUCTION_ACTIVE) {
    std::vector<logical_data_id_t> means_data_ids;
    GetNewLogicalDataID(&means_data_ids, 1);

    {
      GeometricRegion r(0, 0, 0, PARTITION_SIZE * PARTITION_NUM, 1, 1);

      ID<partition_id_t> p(PARTITION_NUM + 0);
      DefinePartition(p, r);

      DefineData(MEANS_DATA_NAME, means_data_ids[0], p.elem(), neighbor_partitions);
    }
  } else {
    std::vector<logical_data_id_t> means_data_ids;
    GetNewLogicalDataID(&means_data_ids, PARTITION_NUM);

    for (size_t i = 0; i < PARTITION_NUM; ++i) {
      GeometricRegion r(i * PARTITION_SIZE, 0, 0, PARTITION_SIZE, 1, 1);

      ID<partition_id_t> p(PARTITION_NUM + i);
      DefinePartition(p, r);

      DefineData(MEANS_DATA_NAME, means_data_ids[i], p.elem(), neighbor_partitions);
    }

    assert((PARTITION_NUM % REDUCTION_PARTITION_NUM) == 0);
    size_t REDUCTION_PARTITION_SIZE = PARTITION_SIZE * PARTITION_NUM / REDUCTION_PARTITION_NUM;
    std::vector<logical_data_id_t> scratch_means_data_ids;
    GetNewLogicalDataID(&scratch_means_data_ids, REDUCTION_PARTITION_NUM);

    for (size_t i = 0; i < REDUCTION_PARTITION_NUM; ++i) {
      GeometricRegion r(i * REDUCTION_PARTITION_SIZE, 0, 0, REDUCTION_PARTITION_SIZE, 1, 1);

      ID<partition_id_t> p(2 * PARTITION_NUM + i);
      DefinePartition(p, r);

      DefineData(SCRATCH_MEANS_DATA_NAME, scratch_means_data_ids[i], p.elem(), neighbor_partitions); // NOLINT
    }
  }


  /*
   * Spawning the init jobs and loop job
   */

  // Spawn the batch of jobs for init stage
  std::vector<job_id_t> init_job_ids;
  GetNewJobID(&init_job_ids, 2 * PARTITION_NUM);
  for (size_t i = 0; i < PARTITION_NUM; ++i) {
    GeometricRegion r(i * PARTITION_SIZE, 0, 0, PARTITION_SIZE, 1, 1);

    read.clear();
    write.clear();
    LoadLdoIdsInSet(&write, r, SAMPLE_BATCH_DATA_NAME, NULL);
    before.clear();
    StageJobAndLoadBeforeSet(&before, INIT_SAMPLES_JOB_NAME, init_job_ids[2*i+1], read, write);
    SerializeParameter(&par, i);
    SpawnComputeJob(INIT_SAMPLES_JOB_NAME, init_job_ids[2*i+1], read, write, before, after, par, true, r); // NOLINT

    if (AUTOMATIC_REDUCTION_ACTIVE && (i > 0)) {
      // if in automatic reduction mode, there is only one means data!
      continue;
    }

    read.clear();
    write.clear();
    LoadLdoIdsInSet(&write, r, MEANS_DATA_NAME, NULL);
    before.clear();
    StageJobAndLoadBeforeSet(&before, INIT_MEANS_JOB_NAME, init_job_ids[2*i], read, write);
    SerializeParameter(&par, i);
    SpawnComputeJob(INIT_MEANS_JOB_NAME, init_job_ids[2*i], read, write, before, after, par, true, r); // NOLINT
  }

  MarkEndOfStage();

  // Spawning loop job
  std::vector<job_id_t> loop_job_id;
  GetNewJobID(&loop_job_id, 1);
  {
    read.clear();
    write.clear();
    before.clear();
    StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, loop_job_id[0], read, write, true);
    SerializeParameter(&par, ITERATION_NUM);
    SpawnComputeJob(LOOP_JOB_NAME, loop_job_id[0], read, write, before, after, par);
  }
};


InitMeans::InitMeans(Application* app) {
  set_application(app);
};

Job * InitMeans::Clone() {
  dbg(DBG_APP, "Cloning init mean job!\n");
  return new InitMeans(application());
};

void InitMeans::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the init means job: %lu\n", id().elem());
  size_t base_val;
  LoadParameter(&params, &base_val);

  assert(da.size() == 1);
  Means *m = static_cast<Means*>(da[0]);
  assert(m->name() == MEANS_DATA_NAME);
  assert(m->cluster_num() == CLUSTER_NUM);

  size_t dimension = m->dimension();
  std::vector<Mean>::iterator iter = m->means()->begin();
  size_t cluster = 0;
  assert(CLUSTER_NUM <= dimension);
  for (; iter != m->means()->end(); ++iter, ++cluster) {
    iter->set_scratch_weight(0);
    for (size_t i = 0; i < dimension; i++) {
      iter->scratch()->operator[](i) = 0;
      // omidm
      iter->vector()->operator[](i) = (i == cluster) ? 1 : 0;
    }
  }
};


InitSamples::InitSamples(Application* app) {
  set_application(app);
};

Job * InitSamples::Clone() {
  dbg(DBG_APP, "Cloning init job!\n");
  return new InitSamples(application());
};

void InitSamples::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the init samples job: %lu\n", id().elem());
  size_t base_val;
  LoadParameter(&params, &base_val);

  assert(da.size() == 1);
  SampleBatch *sb = static_cast<SampleBatch*>(da[0]);
  assert(sb->name() == SAMPLE_BATCH_DATA_NAME);

  size_t dimension = sb->dimension();
  size_t generated_value = 0;
  std::vector<Sample>::iterator iter = sb->samples()->begin();
  for (; iter != sb->samples()->end(); ++iter) {
    for (size_t i = 0; i < dimension; i++) {
      // omidm
      iter->vector()->operator[](i) = generated_value % CLUSTER_NUM;
      generated_value++;
    }
  }
};


ForLoop::ForLoop(Application* app) {
  set_application(app);
};

Job * ForLoop::Clone() {
  dbg(DBG_APP, "Cloning forLoop job!\n");
  return new ForLoop(application());
};

void ForLoop::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the forLoop job: %lu\n", id().elem());

  IDSet<logical_data_id_t> read, write, scratch, reduce;
  IDSet<job_id_t> before, after;
  Parameter par;

  size_t loop_counter;
  LoadParameter(&params, &loop_counter);

  if (loop_counter > 0) {
    StartTemplate("__MARK_STAT_for_loop");

    // differentiate the automatic and manual reduction cases!
    if (AUTOMATIC_REDUCTION_ACTIVE) {
      // Spawn the batch of jobs for cluster stage
      std::vector<job_id_t> cluster_job_ids;
      GetNewJobID(&cluster_job_ids, PARTITION_NUM);
      for (size_t i = 0; i < PARTITION_NUM; ++i) {
        GeometricRegion r(i * PARTITION_SIZE, 0, 0, PARTITION_SIZE, 1, 1);
        read.clear();
        LoadLdoIdsInSet(&read, r, SAMPLE_BATCH_DATA_NAME, MEANS_DATA_NAME, NULL);
        reduce.clear();
        write.clear();
        scratch.clear();
        LoadLdoIdsInSet(&scratch, r, MEANS_DATA_NAME, NULL);
        before.clear();
        StageJobAndLoadBeforeSet(&before, CLUSTER_JOB_NAME, cluster_job_ids[i], read, write, scratch, reduce); // NOLINT
        SpawnComputeJob(CLUSTER_JOB_NAME, cluster_job_ids[i], read, write, scratch, reduce, before, after, par, true, r); // NOLINT
      }

      MarkEndOfStage();

      std::vector<job_id_t> reduction_job_id;
      GetNewJobID(&reduction_job_id, 1);
      {
        GeometricRegion r(0, 0, 0, PARTITION_NUM * PARTITION_SIZE, 1, 1);
        read.clear();
        reduce.clear();
        LoadLdoIdsInSet(&reduce, r, MEANS_DATA_NAME, NULL);
        write.clear();
        LoadLdoIdsInSet(&write, r, MEANS_DATA_NAME, NULL);
        scratch.clear();
        before.clear();
        SerializeParameter(&par, loop_counter - 1);
        StageJobAndLoadBeforeSet(&before, CLUSTER_JOB_NAME, reduction_job_id[0], read, write, scratch, reduce); // NOLINT
        if (REDUCTION_COMBINER_ACTIVE) {
          SpawnComputeJob(REDUCE_JOB_NAME, reduction_job_id[0], read, write, scratch, reduce, before, after, par, COMBINE_JOB_NAME, true, r); // NOLINT
        } else {
          SpawnComputeJob(REDUCE_JOB_NAME, reduction_job_id[0], read, write, scratch, reduce, before, after, par, true, r); // NOLINT
        }
      }

      MarkEndOfStage();
    } else {
      // Spawn the batch of jobs for cluster stage
      std::vector<job_id_t> cluster_job_ids;
      GetNewJobID(&cluster_job_ids, PARTITION_NUM);
      for (size_t i = 0; i < PARTITION_NUM; ++i) {
        GeometricRegion r(i * PARTITION_SIZE, 0, 0, PARTITION_SIZE, 1, 1);
        read.clear();
        LoadLdoIdsInSet(&read, r, SAMPLE_BATCH_DATA_NAME, MEANS_DATA_NAME, NULL);
        write.clear();
        LoadLdoIdsInSet(&write, r, MEANS_DATA_NAME, NULL);
        before.clear();
        StageJobAndLoadBeforeSet(&before, CLUSTER_JOB_NAME, cluster_job_ids[i], read, write);
        SpawnComputeJob(CLUSTER_JOB_NAME, cluster_job_ids[i], read, write, before, after, par, true, r); // NOLINT
      }

      MarkEndOfStage();

      // Spawning the reduction stages
      assert((PARTITION_NUM % REDUCTION_PARTITION_NUM) == 0);
      size_t REDUCTION_PARTITION_SIZE = PARTITION_SIZE * PARTITION_NUM / REDUCTION_PARTITION_NUM;

      std::vector<job_id_t> reduction_l1_job_id;
      GetNewJobID(&reduction_l1_job_id, REDUCTION_PARTITION_NUM);
      for (size_t i = 0; i < REDUCTION_PARTITION_NUM; ++i) {
        GeometricRegion r(i * REDUCTION_PARTITION_SIZE, 0, 0, REDUCTION_PARTITION_SIZE, 1, 1);
        read.clear();
        LoadLdoIdsInSet(&read, r, MEANS_DATA_NAME, NULL);
        write.clear();
        LoadLdoIdsInSet(&write, r, SCRATCH_MEANS_DATA_NAME, NULL);
        before.clear();
        StageJobAndLoadBeforeSet(&before, REDUCE_L1_JOB_NAME, reduction_l1_job_id[i], read, write);
        after.clear();
        SpawnComputeJob(REDUCE_L1_JOB_NAME, reduction_l1_job_id[i], read, write, before, after, par, true, r); // NOLINT
      }

      MarkEndOfStage();

      std::vector<job_id_t> reduction_l2_job_id;
      GetNewJobID(&reduction_l2_job_id, 1);
      {
        GeometricRegion r(0, 0, 0, PARTITION_NUM * PARTITION_SIZE, 1, 1);
        read.clear();
        LoadLdoIdsInSet(&read, r, SCRATCH_MEANS_DATA_NAME, NULL);
        write.clear();
        LoadLdoIdsInSet(&write, r, SCRATCH_MEANS_DATA_NAME, NULL);
        before.clear();
        StageJobAndLoadBeforeSet(&before, REDUCE_L2_JOB_NAME, reduction_l2_job_id[0], read, write);
        after.clear();
        SpawnComputeJob(REDUCE_L2_JOB_NAME, reduction_l2_job_id[0], read, write, before, after, par, true, r); // NOLINT
      }

      MarkEndOfStage();

      std::vector<job_id_t> synch_job_id;
      GetNewJobID(&synch_job_id, REDUCTION_PARTITION_NUM);
      for (size_t i = 0; i < REDUCTION_PARTITION_NUM; ++i) {
        GeometricRegion r(i * REDUCTION_PARTITION_SIZE, 0, 0, REDUCTION_PARTITION_SIZE, 1, 1);
        read.clear();
        LoadLdoIdsInSet(&read, r, SCRATCH_MEANS_DATA_NAME, NULL);
        write.clear();
        LoadLdoIdsInSet(&write, r, MEANS_DATA_NAME, NULL);
        before.clear();
        StageJobAndLoadBeforeSet(&before, SYNCH_JOB_NAME, synch_job_id[i], read, write);
        after.clear();
        SerializeParameter(&par, i * ITERATION_NUM + loop_counter - 1);
        SpawnComputeJob(SYNCH_JOB_NAME, synch_job_id[i], read, write, before, after, par, true, r); // NOLINT
      }

      MarkEndOfStage();
    }

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

    MarkEndOfStage();

    EndTemplate("__MARK_STAT_for_loop");
  } else {
    // StartTemplate("for_loop_end");

    // EndTemplate("for_loop_end");

    TerminateApplication();
  }
};

Cluster::Cluster(Application* app) {
  set_application(app);
};

Job * Cluster::Clone() {
  dbg(DBG_APP, "Cloning cluster job!\n");
  return new Cluster(application());
};

void Cluster::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the cluster job: %lu\n", id().elem());
  assert(da.size() == 3);
  Means *m = NULL;
  SampleBatch *sb = NULL;
  if (da[0]->name() == SAMPLE_BATCH_DATA_NAME) {
    assert(da[1]->name() == MEANS_DATA_NAME);
    m = static_cast<Means*>(da[1]);
    sb = static_cast<SampleBatch*>(da[0]);
  } else {
    assert(da[0]->name() == MEANS_DATA_NAME);
    assert(da[1]->name() == SAMPLE_BATCH_DATA_NAME);
    m = static_cast<Means*>(da[0]);
    sb = static_cast<SampleBatch*>(da[1]);
  }
  assert(da[2]->name() == MEANS_DATA_NAME);
  assert(da[2]->physical_id() == m->physical_id());

  assert(m->means()->size() == CLUSTER_NUM);
  size_t dimension = m->dimension();

  for (size_t i = 0; i < CLUSTER_NUM; ++i) {
    m->means()->operator[](i).set_scratch_weight(0);
    for (size_t k = 0; k < dimension; ++k) {
      m->means()->operator[](i).scratch()->operator[](k) = 0;
    }
  }

  std::vector<Sample>::iterator iter = sb->samples()->begin();
  for (; iter != sb->samples()->end(); ++iter) {
    double distance = -1;
    size_t cluster_index = 0;
    size_t closest_cluster = 0;
    std::vector<Mean>::iterator it = m->means()->begin();
    for (; it != m->means()->end(); ++it, ++cluster_index) {
      double temp = VectorDistance(it->vector(), iter->vector());
      if ((temp < distance) || (distance == -1)) {
        distance = temp;
        closest_cluster = cluster_index;
      }
    }
    VectorAddWithScale(m->means()->operator[](closest_cluster).scratch(), iter->vector(), 1);
    size_t weight = m->means()->operator[](closest_cluster).scratch_weight();
    m->means()->operator[](closest_cluster).set_scratch_weight(++weight);
  }
};



Reduce::Reduce(Application *app) {
  set_application(app);
};

Job * Reduce::Clone() {
  dbg(DBG_APP, "Cloning reduce job!\n");
  return new Reduce(application());
};


void Reduce::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the reduce job: %lu\n", id().elem());
  assert(da.size() >= 1);
  Means *m = static_cast<Means*>(da[0]);
  assert(m->name() == MEANS_DATA_NAME);
  assert(m->cluster_num() == CLUSTER_NUM);
  size_t dimension = m->dimension();

  Means* reduced = static_cast<Means*>(m->Clone());
  reduced->Create();

  for (size_t i = 0; i < CLUSTER_NUM; ++i) {
    size_t sw = 0;
    DataArray::const_iterator iter = da.begin();
    for (size_t j = 0; j < (da.size() - 1); ++j, ++iter) {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == MEANS_DATA_NAME);
      VectorAddWithScale(reduced->means()->operator[](i).vector(),
                        mp->means()->operator[](i).scratch(), 1);
      sw += mp->means()->operator[](i).scratch_weight();
    }

    if (sw > 0) {
      VectorScale(reduced->means()->operator[](i).vector(), 1 / static_cast<double>(sw));
    }

    {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == MEANS_DATA_NAME);
      for (size_t k = 0; k < dimension; ++k) {
        mp->means()->operator[](i).vector()->operator[](k) =
          reduced->means()->operator[](i).vector()->operator[](k);
      }
    }
    iter++;
    assert(iter == da.end());
  }

  {
    DataArray::const_reverse_iterator r_iter = da.rbegin();
    Means *mp = static_cast<Means*>(*r_iter);
    size_t loop_counter;
    LoadParameter(&params, &loop_counter);
    PrintMeans(mp, loop_counter, ITERATION_NUM);
  }
};

Combine::Combine(Application *app) {
  set_application(app);
};

Job * Combine::Clone() {
  dbg(DBG_APP, "Cloning combine job!\n");
  return new Combine(application());
};

void Combine::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the combine job: %lu\n", id().elem());
  assert(da.size() >= (1));

  Means *m = static_cast<Means*>(da[0]);
  assert(m->name() == MEANS_DATA_NAME);
  assert(m->cluster_num() == CLUSTER_NUM);
  size_t dimension = m->dimension();

  Means* reduced = static_cast<Means*>(m->Clone());
  reduced->Create();

  for (size_t i = 0; i < CLUSTER_NUM; ++i) {
    size_t sw = 0;
    DataArray::const_iterator iter = da.begin();
    for (size_t j = 0; j < (da.size() - 1); ++j, ++iter) {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == MEANS_DATA_NAME);
      VectorAddWithScale(reduced->means()->operator[](i).vector(),
                        mp->means()->operator[](i).scratch(), 1);
      sw += mp->means()->operator[](i).scratch_weight();
    }

    {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == MEANS_DATA_NAME);
      mp->means()->operator[](i).set_scratch_weight(sw);
      for (size_t k = 0; k < dimension; ++k) {
        mp->means()->operator[](i).scratch()->operator[](k) =
          reduced->means()->operator[](i).vector()->operator[](k);
      }
    }
    iter++;
    assert(iter == da.end());
  }
};


ReduceL1::ReduceL1(Application *app) {
  set_application(app);
};

Job * ReduceL1::Clone() {
  dbg(DBG_APP, "Cloning reduce l1 job!\n");
  return new ReduceL1(application());
};

void ReduceL1::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the reduce l1 job: %lu\n", id().elem());
  assert((PARTITION_NUM % REDUCTION_PARTITION_NUM) == 0);
  size_t REDUCTION_SIZE = PARTITION_NUM / REDUCTION_PARTITION_NUM;
  assert(da.size() == (REDUCTION_SIZE + 1));

  Means *m = static_cast<Means*>(da[0]);
  assert(m->name() == MEANS_DATA_NAME);
  assert(m->cluster_num() == CLUSTER_NUM);
  size_t dimension = m->dimension();

  Means* reduced = static_cast<Means*>(m->Clone());
  reduced->Create();

  for (size_t i = 0; i < CLUSTER_NUM; ++i) {
    size_t sw = 0;
    DataArray::const_iterator iter = da.begin();
    for (size_t j = 0; j < (da.size() - 1); ++j, ++iter) {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == MEANS_DATA_NAME);
      VectorAddWithScale(reduced->means()->operator[](i).vector(),
                        mp->means()->operator[](i).scratch(), 1);
      sw += mp->means()->operator[](i).scratch_weight();
    }

    {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == SCRATCH_MEANS_DATA_NAME);
      mp->means()->operator[](i).set_scratch_weight(sw);
      for (size_t k = 0; k < dimension; ++k) {
        mp->means()->operator[](i).scratch()->operator[](k) =
          reduced->means()->operator[](i).vector()->operator[](k);
      }
    }
    iter++;
    assert(iter == da.end());
  }
};


ReduceL2::ReduceL2(Application *app) {
  set_application(app);
};

Job * ReduceL2::Clone() {
  dbg(DBG_APP, "Cloning reduce l2 job!\n");
  return new ReduceL2(application());
};

void ReduceL2::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the reduce l2 job: %lu\n", id().elem());
  assert(da.size() == 2 * REDUCTION_PARTITION_NUM);

  Means *m = static_cast<Means*>(da[0]);
  assert(m->name() == SCRATCH_MEANS_DATA_NAME);
  assert(m->cluster_num() == CLUSTER_NUM);
  size_t dimension = m->dimension();

  Means* reduced = static_cast<Means*>(m->Clone());
  reduced->Create();

  for (size_t i = 0; i < CLUSTER_NUM; ++i) {
    size_t sw = 0;
    DataArray::const_iterator iter = da.begin();
    for (size_t j = 0; j < REDUCTION_PARTITION_NUM; ++j, ++iter) {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == SCRATCH_MEANS_DATA_NAME);
      VectorAddWithScale(reduced->means()->operator[](i).vector(),
                        mp->means()->operator[](i).scratch(), 1);
      sw += mp->means()->operator[](i).scratch_weight();
    }

    for (size_t j = 0; j < REDUCTION_PARTITION_NUM; ++j, ++iter) {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == SCRATCH_MEANS_DATA_NAME);
      mp->means()->operator[](i).set_scratch_weight(sw);
      for (size_t k = 0; k < dimension; ++k) {
        mp->means()->operator[](i).scratch()->operator[](k) =
          reduced->means()->operator[](i).vector()->operator[](k);
      }
    }
    assert(iter == da.end());
  }
};


Synch::Synch(Application *app) {
  set_application(app);
};

Job * Synch::Clone() {
  dbg(DBG_APP, "Cloning reduce l3 job!\n");
  return new Synch(application());
};

void Synch::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the synch job: %lu\n", id().elem());
  assert((PARTITION_NUM % REDUCTION_PARTITION_NUM) == 0);
  size_t REDUCTION_SIZE = PARTITION_NUM / REDUCTION_PARTITION_NUM;
  assert(da.size() == (REDUCTION_SIZE + 1));

  for (size_t i = 0; i < CLUSTER_NUM; ++i) {
    DataArray::const_iterator iter = da.begin();
    Means *m = NULL;
    {
      m = static_cast<Means*>(*iter);
      assert(m->name() == SCRATCH_MEANS_DATA_NAME);
    }
    size_t dimension = m->dimension();
    iter++;
    for (size_t j = 0; j < REDUCTION_SIZE; ++j, ++iter) {
      Means *mp = static_cast<Means*>(*iter);
      assert(mp->name() == MEANS_DATA_NAME);

      for (size_t k = 0; k < dimension; ++k) {
        mp->means()->operator[](i).vector()->operator[](k) =
          m->means()->operator[](i).scratch()->operator[](k);
      }
      size_t sw = m->means()->operator[](i).scratch_weight();
      if (sw > 0) {
        VectorScale(mp->means()->operator[](i).vector(), 1 / static_cast<double>(sw));
      }
    }
    assert(iter == da.end());
  }

  {
    DataArray::const_reverse_iterator r_iter = da.rbegin();
    Means *mp = static_cast<Means*>(*r_iter);
    size_t loop_counter;
    LoadParameter(&params, &loop_counter);
    PrintMeans(mp, loop_counter, ITERATION_NUM);
  }
};













