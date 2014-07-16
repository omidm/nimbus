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
 * Iterative projection solving driver. Still in progress.
 *
 * Author: Hang Qu<quhang@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_SMOKE_PROJECTION_PROJECTION_DRIVER_H_
#define NIMBUS_APPLICATION_SMOKE_PROJECTION_PROJECTION_DRIVER_H_

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>

#include "data/physbam/translator_physbam_old.h"
#include "application/smoke/app_utils.h"
#include "application/smoke/options.h"

namespace PhysBAM {

class SPARSE_MATRIX_PARTITION;
typedef GRID<application::TV> T_GRID;

class ProjectionDriver {
 public:
  typedef typename T_GRID::VECTOR_T TV;
  typedef typename TV::SCALAR T;
  typedef typename T_GRID::VECTOR_INT TV_INT;
  typedef application::DataConfig DataConfig;
  typedef application::InitConfig InitConfig;

  ProjectionDriver(
      PCG_SPARSE<T>& pcg_input,
      application::InitConfig& init_config_input,
      application::DataConfig& data_config_input)
          : pcg(pcg_input),
            init_config(init_config_input),
            data_config(data_config_input) {
    projection_data.local_n = 0;
    projection_data.interior_n = 0;

    projection_data.local_tolerance = 0;
    projection_data.global_tolerance = 0;
    projection_data.global_n = 0;
    projection_data.desired_iterations = 0;

    projection_data.local_rho = 0;
    projection_data.rho = 0;
    projection_data.rho_last = 0;
    projection_data.local_dot_product_for_alpha = 0;
    projection_data.alpha = 0;
    projection_data.beta = 0;

    projection_data.local_residual = 0;
    projection_data.residual = 0;
    projection_data.iteration = 0;
    projection_data.matrix_a = NULL;
    projection_data.matrix_index_to_cell_index = NULL;
    cache_pressure = NULL;
    cache_vector_p = NULL;
    cache_matrix_a = NULL;
    cache_matrix_c = NULL;
    cache_index_m2c = NULL;
  }

  virtual ~ProjectionDriver() {
    if (projection_data.matrix_a) {
      delete projection_data.matrix_a;
    }
    if (projection_data.matrix_index_to_cell_index) {
      delete projection_data.matrix_index_to_cell_index;
    }
  }

  class ProjectionData {
   public:
    ARRAY<float, VECTOR<int,3> > pressure;
    ARRAY<float, VECTOR<int,3> > grid_format_vector_p;
    ARRAY<TV_INT>* matrix_index_to_cell_index;
    ARRAY<int, TV_INT> cell_index_to_matrix_index;
    SPARSE_MATRIX_FLAT_NXN<T>* matrix_a;
    SPARSE_MATRIX_FLAT_NXN<T>* matrix_c;
    VECTOR_ND<T> vector_b;
    VECTOR_ND<T> vector_x;
    // VECTOR_ND<T> x_interior;
    VECTOR_ND<T> temp, temp_interior;
    VECTOR_ND<T> p, p_interior;
    VECTOR_ND<T> z_interior;
    VECTOR_ND<T> b_interior;

    int local_n;
    int interior_n;  // SUM.

    T local_tolerance;  // MAX.
    T global_tolerance;
    int global_n;
    int desired_iterations;

    double local_residual;  // SUM.
    double local_rho, rho, rho_last;  // SUM.
    double local_dot_product_for_alpha;  // SUM.
    T alpha;
    T beta;

    // Not a Nimbus data type.
    T residual;

    // Passed by parameter.
    int iteration;
  } projection_data;

  nimbus::TranslatorPhysBAMOld<TV> translator;
  PCG_SPARSE<T>& pcg;
  InitConfig& init_config;
  DataConfig& data_config;

  SPARSE_MATRIX_PARTITION partition;
  typedef typename application::CacheScalarArray<T> TCacheScalarArray;
  TCacheScalarArray *cache_pressure;
  TCacheScalarArray *cache_vector_p;
  application::CacheSparseMatrix *cache_matrix_a;
  application::CacheSparseMatrix *cache_matrix_c;
  application::CacheArrayM2C * cache_index_m2c;

  template<class TYPE> TYPE Global_Sum(const TYPE& input) {
    return input;
  }
  template<class TYPE> TYPE Global_Max(const TYPE& input) {
    return input;
  }

  void Initialize(int local_n, int interior_n);
  void Cache_Initialize(int local_n, int interior_n);
  void LocalInitialize();
  void GlobalInitialize();
  void InitializeResidual();
  bool SpawnFirstIteration();
  void DoPrecondition();
  void CalculateLocalRho();
  void ReduceRho();
  void UpdateSearchVector();
  void UpdateTempVector();
  void CalculateLocalAlpha();
  void ReduceAlpha();
  void UpdateOtherVectors();
  void CalculateLocalResidual();
  bool DecideToSpawnNextIteration();
  void LoadFromNimbus(const nimbus::Job* job, const nimbus::DataArray& da);
  void Cache_LoadFromNimbus(const nimbus::Job* job, const nimbus::DataArray& da);
  void SaveToNimbus(const nimbus::Job* job, const nimbus::DataArray& da);
  void Cache_SaveToNimbus(const nimbus::Job* job, const nimbus::DataArray& da);
  template<typename TYPE_NAME> void ReadScalarData(
      const nimbus::Job* job, const nimbus::DataArray& da,
      const char* variable_name, TYPE_NAME& value);
  void ReadVectorData(
      const nimbus::Job* job, const nimbus::DataArray& da,
      const char* variable_name, VECTOR_ND<float>& value);
  template<typename TYPE_NAME> void WriteScalarData(
      const nimbus::Job* job, const nimbus::DataArray& da,
      const char* variable_name, const TYPE_NAME& value);
  void WriteVectorData(
      const nimbus::Job* job, const nimbus::DataArray& da,
      const char* variable_name, const VECTOR_ND<float>& value);
};

}  // namespace PhysBAM

#endif  // NIMBUS_APPLICATION_SMOKE_PROJECTION_PROJECTION_DRIVER_H_
