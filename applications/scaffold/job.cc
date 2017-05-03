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
 * Job classes are implemented here.
 *
 */

#include "./job.h"
#include "./data.h"
#include "./utils.h"


#define PART_NUM static_cast<ScaffoldApp*>(application())->part_num_

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  dbg(DBG_APP, "Cloning main job!\n");
  return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the main job\n");
  printf("Main message is %s.\n",
      static_cast<ScaffoldApp*>(application())->message_.c_str());

  /*
   * Defining partitions and data.
   */
  std::vector<logical_data_id_t> data_ids;
  GetNewLogicalDataID(&data_ids, PART_NUM);
  IDSet<partition_id_t> neighbor_partitions;
  for (size_t i = 0; i < PART_NUM; ++i) {
    GeometricRegion region(i, 0, 0, 1, 1, 1);
    DefinePartition(ID<partition_id_t>(i), region);
    DefineData(DATA_NAME, data_ids[i], i, neighbor_partitions);
  }

  /*
   * Spawning jobs
   */
  std::vector<job_id_t> job_ids;
  GetNewJobID(&job_ids, PART_NUM);
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;
  for (size_t i = 0; i < PART_NUM; ++i) {
    GeometricRegion region(i, 0, 0, 1, 1, 1);
    read.clear();
    LoadLdoIdsInSet(&read, region, DATA_NAME, NULL);
    write.clear();
    LoadLdoIdsInSet(&write, region, DATA_NAME, NULL);
    before.clear();
    StageJobAndLoadBeforeSet(&before, PRINT_JOB_NAME, job_ids[i], read, write);
    SerializeParameter(&par, i);
    SpawnComputeJob(PRINT_JOB_NAME, job_ids[i], read, write, before, after, par, true, region);
  }

  MarkEndOfStage();

  TerminateApplication();
};

Print::Print() {
};

Job * Print::Clone() {
  dbg(DBG_APP, "Cloning print job!\n");
  return new Print();
};

void Print::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the print job\n");
  size_t part_num;
  LoadParameter(&params, &part_num);
  printf("PRINT: PART %lu\n", part_num);
};


