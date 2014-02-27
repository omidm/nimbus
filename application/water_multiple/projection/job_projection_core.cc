/* Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 vd* are met:
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
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <sstream>
#include <string>

#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/projection/nimbus_pcg_sparse_mpi.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "data/scalar_data.h"
#include "application/water_multiple/data_include.h"
#include "application/water_multiple/projection/job_projection_core.h"

namespace application {

JobProjectionCore::JobProjectionCore(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobProjectionCore::Clone() {
  return new JobProjectionCore(application());
}

void JobProjectionCore::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing PROJECTION_CORE job.\n");

  InitConfig init_config;
  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                &init_config.global_region, &init_config.local_region);

  // Assume time, dt, frame is ready from here.
  dbg(APP_LOG,
      "In PROJECTION: Initialize WATER_DRIVER/WATER_EXAMPLE"
      "(Frame=%d, Time=%f).\n",
      init_config.frame, init_config.time);

  DataConfig data_config;
  data_config.SetFlag(DataConfig::MATRIX_A);
  data_config.SetFlag(DataConfig::VECTOR_X);
  data_config.SetFlag(DataConfig::VECTOR_B);
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE);
  data_config.SetFlag(DataConfig::INDEX_M2C);
  data_config.SetFlag(DataConfig::INDEX_C2M);

  // MPI reference version:
  // laplace_mpi->Solve(A, x, b, q, s, r, k, z, tolerance, color);
  // color only used for MPI version.
  // laplace->pcg.Solve(A, x, b, q, s, r, k, z, laplace->tolerance);
  PhysBAM::PCG_SPARSE<float> pcg_temp;
  pcg_temp.Set_Maximum_Iterations(40);
  pcg_temp.evolution_solver_type = PhysBAM::krylov_solver_cg;
  pcg_temp.cg_restart_iterations = 0;
  pcg_temp.Show_Results();

  PhysBAM::NIMBUS_PCG_SPARSE_MPI pcg_mpi(pcg_temp);
  dbg(APP_LOG, "Job PROJECTION_CORE starts (dt=%f).\n", dt);

  {
    nimbus::Job* job = this;
    typedef nimbus::Data Data;
    if (data_config.GetFlag(DataConfig::MATRIX_A)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_MATRIX_A), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataSparseMatrix* data_real =
            dynamic_cast<application::DataSparseMatrix*>(data_temp);
        data_real->LoadFromNimbus(&pcg_mpi.projection_data.matrix_a);
        dbg(APP_LOG, "Finish reading MATRIX_A\n");
      }
    }
    if (data_config.GetFlag(DataConfig::VECTOR_X)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_VECTOR_X), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawVectorNd* data_real =
            dynamic_cast<application::DataRawVectorNd*>(data_temp);
        data_real->LoadFromNimbus(&pcg_mpi.projection_data.vector_x);
        dbg(APP_LOG, "Finish reading VECTOR_X\n");
      }
    }
    if (data_config.GetFlag(DataConfig::VECTOR_B)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_VECTOR_B), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawVectorNd* data_real =
            dynamic_cast<application::DataRawVectorNd*>(data_temp);
        data_real->LoadFromNimbus(&pcg_mpi.projection_data.vector_b);
        dbg(APP_LOG, "Finish reading VECTOR_B\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_C2M), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawGridArray* data_real =
            dynamic_cast<application::DataRawGridArray*>(data_temp);
        pcg_mpi.projection_data.cell_index_to_matrix_index.Resize(
            PhysBAM::RANGE<TV_INT>(TV_INT(0, 0, 0),
                                   TV_INT(init_config.local_region.dx()+1,
                                          init_config.local_region.dy()+1,
                                          init_config.local_region.dz()+1)));

        data_real->LoadFromNimbus(
            &pcg_mpi.projection_data.cell_index_to_matrix_index);
        dbg(APP_LOG, "Finish reading INDEX_C2M\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_M2C), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawArrayM2C* data_real =
            dynamic_cast<application::DataRawArrayM2C*>(data_temp);
        data_real->LoadFromNimbus(
            &pcg_mpi.projection_data.matrix_index_to_cell_index);
        dbg(APP_LOG, "Finish reading INDEX_M2C\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_LOCAL_TOLERANCE),
          da, application::READ_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
        pcg_mpi.projection_data.local_tolerance = data_real->scalar();
        dbg(APP_LOG, "Finish reading tolerance\n");
      }
    }
  }
  /*
  pcg_mpi.projection_data.matrix_index_to_cell_index =
      &laplace_solver_wrapper.matrix_index_to_cell_index_array(1);
  pcg_mpi.projection_data.cell_index_to_matrix_index =
      &laplace_solver_wrapper.cell_index_to_matrix_index;
  pcg_mpi.projection_data.matrix_a = &laplace_solver_wrapper.A_array(1);
  pcg_mpi.projection_data.vector_b = &laplace_solver_wrapper.b_array(1);
  pcg_mpi.projection_data.vector_x = &laplace_solver_wrapper.x;
  pcg_mpi.projection_data.local_tolerance = laplace_solver.tolerance;
  */
  pcg_mpi.Initialize();
  pcg_mpi.CommunicateConfig();
  pcg_mpi.Parallel_Solve();

  if (data_config.GetFlag(DataConfig::VECTOR_X)) {
    Data* data_temp = application::GetTheOnlyData(
        this, std::string(APP_VECTOR_X), da, application::WRITE_ACCESS);
    if (data_temp) {
      application::DataRawVectorNd* data_real =
          dynamic_cast<application::DataRawVectorNd*>(data_temp);
      data_real->SaveToNimbus(pcg_mpi.projection_data.vector_x);
    }
  }

  dbg(APP_LOG, "Completed executing PROJECTION_CORE job\n");
}

}  // namespace application
