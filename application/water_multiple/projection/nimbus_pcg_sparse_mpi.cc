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

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>

#include "application/water_multiple/projection/nimbus_pcg_sparse_mpi.h"

namespace PhysBAM {

void NIMBUS_PCG_SPARSE_MPI::Initialize() {
  int local_n = (*projection_data.matrix_a).n;
  // [TODO] Not accurate.
  partition.interior_indices.min_corner = 1;
  partition.interior_indices.max_corner = local_n;
  int interior_n = partition.interior_indices.Size()+1;

  projection_data.x_interior = new VECTOR_ND<T>;
  projection_data.x_interior->Set_Subvector_View(
      *projection_data.vector_x,
      partition.interior_indices);

  projection_data.temp = new VECTOR_ND<T>(local_n, false);
  projection_data.temp_interior = new VECTOR_ND<T>;
  projection_data.temp_interior->Set_Subvector_View(
      *projection_data.temp,
      partition.interior_indices);

  projection_data.p = new VECTOR_ND<T>(local_n, false);
  projection_data.p_interior = new VECTOR_ND<T>;
  projection_data.p_interior->Set_Subvector_View(
      *projection_data.p,
      partition.interior_indices);

  projection_data.z_interior = new VECTOR_ND<T>(interior_n, false);

  projection_data.b_interior = new VECTOR_ND<T>;
  projection_data.b_interior->Set_Subvector_View(
      *projection_data.vector_b,
      partition.interior_indices);

  projection_data.alpha = 0;
  projection_data.beta = 0;
  projection_data.rho = 0;
  projection_data.rho_last = 0;
  projection_data.residual = 0;

  projection_data.iteration = 0;
}

void NIMBUS_PCG_SPARSE_MPI::CommunicateConfig() {
  int interior_n = partition.interior_indices.Size()+1;

  int global_n = Global_Sum(interior_n);
  projection_data.global_n = global_n;

  T global_tolerance = Global_Max(projection_data.local_tolerance);
  projection_data.global_tolerance = global_tolerance;

  int global_desired_iterations = projection_data.global_n;
  if (pcg.maximum_iterations) {
    global_desired_iterations =
        min(global_desired_iterations, pcg.maximum_iterations);
  }

  projection_data.global_desired_iterations = global_desired_iterations;
}

void NIMBUS_PCG_SPARSE_MPI::Parallel_Solve() {
  ExchangePressure();
  InitializeResidual();
  if (SpawnFirstIteration()) {
    do {
      projection_data.iteration++;
      DoPrecondition();
      CalculateBeta();
      UpdateSearchVector();
      ExchangeSearchVector();
      UpdateTempVector();
      CalculateAlpha();
      UpdateOtherVectors();
      CalculateResidual();
    } while (DecideToSpawnNextIteration());
    ExchangePressure();
  }
}

// Projection is broken to "smallest" code piece to allow future changes.

void NIMBUS_PCG_SPARSE_MPI::ExchangePressure() {
  VECTOR_ND<T>& x = (*projection_data.vector_x);
  Fill_Ghost_Cells(x);
}

void NIMBUS_PCG_SPARSE_MPI::InitializeResidual() {
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data.matrix_a);
  VECTOR_ND<T>& x = (*projection_data.vector_x);
  VECTOR_ND<T>& temp = (*projection_data.temp);
  VECTOR_ND<T>& b_interior = (*projection_data.b_interior);
  VECTOR_ND<T>& temp_interior = (*projection_data.temp_interior);
  A.Times(x, temp);
  b_interior -= temp_interior;
}

bool NIMBUS_PCG_SPARSE_MPI::SpawnFirstIteration() {
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data.matrix_a);
  const bool recompute_preconditioner = true;
  VECTOR_ND<T>& b_interior = (*projection_data.b_interior);

  double local_norm = b_interior.Max_Abs();
  if (Global_Max(local_norm) <= projection_data.global_tolerance) {
    return false;
  }

  if (pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)) {
      delete A.C;
      A.C = A.Create_Submatrix(partition.interior_indices);
      A.C->In_Place_Incomplete_Cholesky_Factorization(
          pcg.modified_incomplete_cholesky,
          pcg.modified_incomplete_cholesky_coefficient,
          pcg.preconditioner_zero_tolerance,
          pcg.preconditioner_zero_replacement);
  }
  return true;
}

void NIMBUS_PCG_SPARSE_MPI::DoPrecondition() {
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data.matrix_a);
  VECTOR_ND<T>& z_interior = (*projection_data.z_interior);
  VECTOR_ND<T>& b_interior = (*projection_data.b_interior);
  VECTOR_ND<T>& temp_interior = (*projection_data.temp_interior);
  A.C->Solve_Forward_Substitution(b_interior, temp_interior, true);
  A.C->Solve_Backward_Substitution(temp_interior, z_interior, false, true);
}

void NIMBUS_PCG_SPARSE_MPI::CalculateBeta() {
  VECTOR_ND<T>& z_interior = (*projection_data.z_interior);
  VECTOR_ND<T>& b_interior = (*projection_data.b_interior);
  projection_data.rho_last = projection_data.rho;
  projection_data.rho = Global_Sum(
      VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior, b_interior));
  projection_data.beta = (T)(projection_data.rho / projection_data.rho_last);
}

void NIMBUS_PCG_SPARSE_MPI::UpdateSearchVector() {
  int interior_n = partition.interior_indices.Size() + 1;
  VECTOR_ND<T>& z_interior = (*projection_data.z_interior);
  VECTOR_ND<T>& p_interior = (*projection_data.p_interior);
  if (projection_data.iteration == 1) {
    p_interior = z_interior;
  } else {
    for (int i = 1; i <= interior_n; i++)
      p_interior(i) = z_interior(i) + projection_data.beta * p_interior(i);
  }
}

void NIMBUS_PCG_SPARSE_MPI::ExchangeSearchVector() {
  VECTOR_ND<T>& p = (*projection_data.p);
  Fill_Ghost_Cells(p);
}

void NIMBUS_PCG_SPARSE_MPI::UpdateTempVector() {
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data.matrix_a);
  VECTOR_ND<T>& temp = (*projection_data.temp);
  VECTOR_ND<T>& p = (*projection_data.p);
  A.Times(p, temp);
}

void NIMBUS_PCG_SPARSE_MPI::CalculateAlpha() {
  VECTOR_ND<T>& p_interior = (*projection_data.p_interior);
  VECTOR_ND<T>& temp_interior = (*projection_data.temp_interior);
  projection_data.alpha =
      (T) (projection_data.rho /
           Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(
                   p_interior, temp_interior)));
}

void NIMBUS_PCG_SPARSE_MPI::UpdateOtherVectors() {
  int interior_n = partition.interior_indices.Size()+1;
  VECTOR_ND<T>& x_interior = (*projection_data.x_interior);
  VECTOR_ND<T>& p_interior = (*projection_data.p_interior);
  VECTOR_ND<T>& b_interior = (*projection_data.b_interior);
  VECTOR_ND<T>& temp_interior = (*projection_data.temp_interior);
  for (int i = 1; i <= interior_n; i++) {
    x_interior(i) += projection_data.alpha * p_interior(i);
    b_interior(i) -= projection_data.alpha * temp_interior(i);
  }
}

void NIMBUS_PCG_SPARSE_MPI::CalculateResidual() {
  VECTOR_ND<T>& b_interior = (*projection_data.b_interior);
  double local_norm = b_interior.Max_Abs();
  projection_data.residual = Global_Max(local_norm);
}

bool NIMBUS_PCG_SPARSE_MPI::DecideToSpawnNextIteration() {
  if (projection_data.residual <= projection_data.global_tolerance) {
    return false;
  }
  if (projection_data.iteration == projection_data.global_desired_iterations) {
    return false;
  }
  return true;
}

}  // namespace PhysBAM
