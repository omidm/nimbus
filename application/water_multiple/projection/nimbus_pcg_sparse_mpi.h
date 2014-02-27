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
 * Helper function for projection job. Still in progress.
 * Author: Hang Qu<quhang@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_NIMBUS_PCG_SPARSE_MPI_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_NIMBUS_PCG_SPARSE_MPI_H_

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>

#include "application/water_multiple/app_utils.h"

namespace PhysBAM {

class SPARSE_MATRIX_PARTITION;
typedef GRID<application::TV> T_GRID;

class NIMBUS_PCG_SPARSE_MPI {
 public:
  typedef typename T_GRID::VECTOR_T TV;
  typedef typename TV::SCALAR T;
  typedef typename T_GRID::VECTOR_INT TV_INT;

  NIMBUS_PCG_SPARSE_MPI(PCG_SPARSE<T>& pcg_input) : pcg(pcg_input) {
    projection_data.x_interior = NULL;
    projection_data.temp = NULL;
    projection_data.temp_interior = NULL;
    projection_data.p = NULL;
    projection_data.p_interior = NULL;
    projection_data.z_interior = NULL;
    projection_data.b_interior = NULL;
    projection_data.global_n = 0;
    projection_data.local_tolerance = 0;
    projection_data.global_tolerance = 0;
    projection_data.global_desired_iterations = 0;
    projection_data.alpha = 0;
    projection_data.beta = 0;
    projection_data.rho = 0;
    projection_data.rho_last = 0;
    projection_data.residual = 0;
    projection_data.iteration = 0;
  }

  virtual ~NIMBUS_PCG_SPARSE_MPI() {
    delete projection_data.x_interior;
    delete projection_data.temp;
    delete projection_data.temp_interior;
    delete projection_data.p;
    delete projection_data.p_interior;
    delete projection_data.z_interior;
    delete projection_data.b_interior;
  }

  class ProjectionData {
   public:
    typedef VECTOR_ND<T>* P_VECTOR;
    ARRAY<TV_INT> matrix_index_to_cell_index;
    ARRAY<int, TV_INT> cell_index_to_matrix_index;
    SPARSE_MATRIX_FLAT_NXN<T> matrix_a;
    VECTOR_ND<T> vector_b;
    VECTOR_ND<T> vector_x;
    P_VECTOR x_interior;
    P_VECTOR temp, temp_interior;
    P_VECTOR p, p_interior;
    P_VECTOR z_interior;
    P_VECTOR b_interior;

    // Static config.
    int global_n;
    T local_tolerance;
    T global_tolerance;
    int global_desired_iterations;

    double rho, rho_last;
    T beta, alpha;
    T residual;

    int iteration;
  } projection_data;

  PCG_SPARSE<T>& pcg;
  // [TODO]
  SPARSE_MATRIX_PARTITION partition;

  // [TODO] Global sum should not use MPI.
  template<class TYPE> TYPE Global_Sum(const TYPE& input) {
    return input;
  }

  // [TODO] Global max should not use MPI.
  template<class TYPE> TYPE Global_Max(const TYPE& input) {
    return input;
  }

  // [TODO] Ghost cell transmission should not use MPI.
  virtual void Fill_Ghost_Cells(VECTOR_ND<T>& v) {
    return;
  }

  void Initialize();
  void CommunicateConfig();
  void Parallel_Solve();
  void ExchangePressure();
  void InitializeResidual();
  bool SpawnFirstIteration();
  void DoPrecondition();
  void CalculateBeta();
  void UpdateSearchVector();
  void ExchangeSearchVector();
  void UpdateTempVector();
  void CalculateAlpha();
  void UpdateOtherVectors();
  void CalculateResidual();
  bool DecideToSpawnNextIteration();
};
}  // namespace PhysBAM

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_NIMBUS_PCG_SPARSE_MPI_H_
