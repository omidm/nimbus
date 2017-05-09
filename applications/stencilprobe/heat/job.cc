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
#include "src/application_utils/data_definer.h"


#define NX static_cast<Heat*>(application())->nx_
#define NY static_cast<Heat*>(application())->ny_
#define NZ static_cast<Heat*>(application())->nz_
#define PNX static_cast<Heat*>(application())->pnx_
#define PNY static_cast<Heat*>(application())->pny_
#define PNZ static_cast<Heat*>(application())->pnz_
#define PART_NUM PNX*PNY*PNZ
#define ITER_NUM static_cast<Heat*>(application())->iter_num_

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  dbg(DBG_APP, "Cloning main job!\n");
  return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the main job\n");

  /*
   * Defining partitions and data.
   */
  nimbus::DataDefiner df(this);

  df.DefineData(DATA_NAME,
                NX, NY, NZ,
                1, 1, 1,
                PNX, PNY, PNZ,
                false);  // No global boundary.

  /*
   * Spawning jobs
   */
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;

  // // Spawning init jobs
  // std::vector<job_id_t> job_ids;
  // GetNewJobID(&job_ids, 1);
  // for (size_t i = 0; i < PART_NUM; ++i) {
  //   GeometricRegion region(i, 0, 0, 1, 1, 1);
  //   read.clear();
  //   LoadLdoIdsInSet(&read, region, DATA_NAME, NULL);
  //   write.clear();
  //   LoadLdoIdsInSet(&write, region, DATA_NAME, NULL);
  //   before.clear();
  //   StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, job_ids[i], read, write);
  //   SerializeParameter(&par, i);
  //   SpawnComputeJob(LOOP_JOB_NAME, job_ids[i], read, write, before, after, par, true, region);
  // }

  // MarkEndOfStage();

  // Spawning loop job
  std::vector<job_id_t> loop_job_id;
  GetNewJobID(&loop_job_id, 1);
  {
    read.clear();
    write.clear();
    before.clear();
    StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, loop_job_id[0], read, write, true);
    SerializeParameter(&par, ITER_NUM);
    SpawnComputeJob(LOOP_JOB_NAME, loop_job_id[0], read, write, before, after, par);
  }
};


Loop::Loop(Application* app) {
  set_application(app);
};

Job * Loop::Clone() {
  dbg(DBG_APP, "Cloning loop job!\n");
  return new Loop(application());
};

void Loop::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the loop job\n");
  size_t iter_num;
  LoadParameter(&params, &iter_num);
  printf("Iteration number remaining: %lu\n", iter_num);

  if (iter_num > 0) {
    StartTemplate("__MARK_STAT_for_loop");
    /*
     * Spawning stencil jobs
     */
    std::vector<job_id_t> stencil_job_ids;
    GetNewJobID(&stencil_job_ids, PART_NUM);
    IDSet<logical_data_id_t> read, write;
    IDSet<job_id_t> before, after;
    Parameter par;
    for (int i = 0; i < PART_NUM; ++i) {
      read.clear();
      LoadLdoIdsInSet(&read, ph.map()["kRegionsOuter"][i], DATA_NAME, NULL);
      write.clear();
      LoadLdoIdsInSet(&write, ph.map()["kRegionsCentral"][i], DATA_NAME, NULL);
      before.clear();
      StageJobAndLoadBeforeSet(&before, STENCIL_JOB_NAME, stencil_job_ids[i], read, write);
      SerializeParameter(&par, i);
      SpawnComputeJob(STENCIL_JOB_NAME, stencil_job_ids[i], read, write, before, after, par, true, ph.map()["kRegionsCentral"][i]);
      printf("SPAWNED STENCIL\n");
    }

    MarkEndOfStage();

    /*
     * Spawning loop job
     */
    std::vector<job_id_t> loop_job_id;
    GetNewJobID(&loop_job_id, 1);
    {
      read.clear();
      write.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, LOOP_JOB_NAME, loop_job_id[0], read, write, true);
      SerializeParameter(&par, iter_num - 1);
      SpawnComputeJob(LOOP_JOB_NAME, loop_job_id[0], read, write, before, after, par);
    }
    EndTemplate("__MARK_STAT_for_loop");
  } else {  // if (iter_num > 0)

    TerminateApplication();
  }
};


Stencil::Stencil(Application *app) {
  set_application(app);
};

Job * Stencil::Clone() {
  dbg(DBG_APP, "Cloning stencil job!\n");
  return new Stencil(application());
};

void Stencil::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Executing the stencil job\n");
  size_t part_num;
  LoadParameter(&params, &part_num);
  printf("PRINT: Stencil %lu\n", part_num);

  // DataArray::const_iterator iter = da.begin();
  // std::cout << "There are " << da.size() << std::endl;
  // for (; iter != da.end(); ++ iter) {
  //   std::cout << (*iter)->region().ToNetworkData() << std::endl;
  // }
  

  nimbus::DataArray read, write;
  GetDataAccess(this, da, READ, &read);
  GetDataAccess(this, da, WRITE, &write);
  printf("Read size: %lu, Write size: %lu\n", read.size(), write.size());

  nimbus::AppDataManager *cm = this->GetAppDataManager();
  nimbus::AppVar *app_var =
    cm->GetAppVar(
        read, ph.map()["kRegionsOuter"][part_num],
        write, ph.map()["kRegionsCentral"][part_num],
        AppDataVecPrototype, ph.map()["kRegionsCentral"][part_num],
        nimbus::app_data::SHARED);
  double *data = dynamic_cast<AppDataVec*>(app_var)->data();
  // Compute
  cm->ReleaseAccess(app_var);


};


