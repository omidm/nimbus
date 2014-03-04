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

#include <cassert>

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>

#include "application/water_multiple/data_include.h"
#include "data/scalar_data.h"
#include "application/water_multiple/projection/projection_driver.h"

namespace PhysBAM {

// local_n is the dimentsion of matrixa_a, interior_n is the size of the
// internal.
void ProjectionDriver::Initialize(int local_n, int interior_n) {
  partition.interior_indices.min_corner = 1;
  partition.interior_indices.max_corner = interior_n;
  // int interior_n = partition.interior_indices.Size()+1;

  if (projection_data.temp.Size() == 0 &&
      data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    projection_data.temp.Resize(local_n, false);
  }
  if (projection_data.p.Size() == 0 &&
      data_config.GetFlag(DataConfig::VECTOR_P)) {
    projection_data.p.Resize(local_n, false);
  }
  if (projection_data.z_interior.Size() == 0 &&
      data_config.GetFlag(DataConfig::VECTOR_Z)) {
    projection_data.z_interior.Resize(interior_n, false);
  }

  if (data_config.GetFlag(DataConfig::VECTOR_X)) {
    projection_data.x_interior.Set_Subvector_View(
        projection_data.vector_x,
        partition.interior_indices);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    projection_data.temp_interior.Set_Subvector_View(
        projection_data.temp,
        partition.interior_indices);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_P)) {
    projection_data.p_interior.Set_Subvector_View(
        projection_data.p,
        partition.interior_indices);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    projection_data.b_interior.Set_Subvector_View(
        projection_data.vector_b,
        partition.interior_indices);
  }
}

// Projection is broken to "smallest" code piece to allow future changes.
void ProjectionDriver::LocalInitialize() {
  SPARSE_MATRIX_FLAT_NXN<T>& A = projection_data.matrix_a;
  VECTOR_ND<T>& x = projection_data.vector_x;
  VECTOR_ND<T>& temp = projection_data.temp;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  A.Times(x, temp);
  b_interior -= temp_interior;
  projection_data.local_residual = b_interior.Max_Abs();
  const bool recompute_preconditioner = true;
  if (pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)) {
      delete A.C;
      A.C = A.Create_Submatrix(partition.interior_indices);
      A.C->In_Place_Incomplete_Cholesky_Factorization(
          pcg.modified_incomplete_cholesky,
          pcg.modified_incomplete_cholesky_coefficient,
          pcg.preconditioner_zero_tolerance,
          pcg.preconditioner_zero_replacement);
  }
  projection_data.temp.Resize(projection_data.local_n, false);
  projection_data.temp.Fill(0);
  projection_data.p.Resize(projection_data.local_n, false);
  projection_data.p.Fill(0);
  projection_data.z_interior.Resize(projection_data.interior_n, false);
  projection_data.z_interior.Fill(0);
}

void ProjectionDriver::GlobalInitialize() {
  // projection_data.interior_n = partition.interior_indices.Size()+1;
  projection_data.global_n = Global_Sum(projection_data.interior_n);
  projection_data.global_tolerance =
      Global_Max(projection_data.local_tolerance);

  projection_data.desired_iterations = projection_data.global_n;
  if (pcg.maximum_iterations) {
    projection_data.desired_iterations =
        min(projection_data.desired_iterations, pcg.maximum_iterations);
  }
}

void ProjectionDriver::ExchangePressure() {
  VECTOR_ND<T>& x = projection_data.vector_x;
  Fill_Ghost_Cells(x);
}

void ProjectionDriver::DoPrecondition() {
  SPARSE_MATRIX_FLAT_NXN<T>& A = projection_data.matrix_a;
  VECTOR_ND<T>& z_interior = projection_data.z_interior;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  A.C->Solve_Forward_Substitution(b_interior, temp_interior, true);
  A.C->Solve_Backward_Substitution(temp_interior, z_interior, false, true);
}

void ProjectionDriver::CalculateLocalRho() {
  VECTOR_ND<T>& z_interior = projection_data.z_interior;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  projection_data.local_rho =
      VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior, b_interior);
}

void ProjectionDriver::ReduceRho() {
  projection_data.rho_last = projection_data.rho;
  // TODO(quhang), change.
  projection_data.rho = Global_Sum(projection_data.local_rho);
  projection_data.beta = (T)(projection_data.rho / projection_data.rho_last);
}

void ProjectionDriver::UpdateSearchVector() {
  int interior_n = partition.interior_indices.Size() + 1;
  VECTOR_ND<T>& z_interior = projection_data.z_interior;
  VECTOR_ND<T>& p_interior = projection_data.p_interior;
  if (projection_data.iteration == 1) {
    p_interior = z_interior;
  } else {
    for (int i = 1; i <= interior_n; i++)
      p_interior(i) = z_interior(i) + projection_data.beta * p_interior(i);
  }
}

void ProjectionDriver::ExchangeSearchVector() {
  VECTOR_ND<T>& p = projection_data.p;
  Fill_Ghost_Cells(p);
}

void ProjectionDriver::UpdateTempVector() {
  SPARSE_MATRIX_FLAT_NXN<T>& A = projection_data.matrix_a;
  VECTOR_ND<T>& temp = projection_data.temp;
  VECTOR_ND<T>& p = projection_data.p;
  A.Times(p, temp);
}

void ProjectionDriver::CalculateLocalAlpha() {
  VECTOR_ND<T>& p_interior = projection_data.p_interior;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  projection_data.local_dot_product_for_alpha =
      VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior, temp_interior);
}

void ProjectionDriver::ReduceAlpha() {
  projection_data.alpha =
      (T) (projection_data.rho /
           Global_Sum(projection_data.local_dot_product_for_alpha));
}

void ProjectionDriver::UpdateOtherVectors() {
  int interior_n = partition.interior_indices.Size()+1;
  VECTOR_ND<T>& x_interior = projection_data.x_interior;
  VECTOR_ND<T>& p_interior = projection_data.p_interior;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  for (int i = 1; i <= interior_n; i++) {
    x_interior(i) += projection_data.alpha * p_interior(i);
    b_interior(i) -= projection_data.alpha * temp_interior(i);
  }
}

void ProjectionDriver::CalculateLocalResidual() {
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  projection_data.local_residual = b_interior.Max_Abs();
}

bool ProjectionDriver::DecideToSpawnNextIteration() {
  if (projection_data.residual <= projection_data.global_tolerance) {
    return false;
  }
  if (projection_data.iteration == projection_data.desired_iterations) {
    return false;
  }
  return true;
}

void ProjectionDriver::LoadFromNimbus(
    const nimbus::Job* job, const nimbus::DataArray& da) {
  // MATRIX_A. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::MATRIX_A)) {
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_MATRIX_A), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataSparseMatrix* data_real =
          dynamic_cast<application::DataSparseMatrix*>(data_temp);
      data_real->LoadFromNimbus(&projection_data.matrix_a);
      dbg(APP_LOG, "Finish reading MATRIX_A.\n");
    }
  }
  // VECTOR_B. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    ReadVectorData(job, da, APP_VECTOR_B, projection_data.vector_b);
  }
  // VECTOR_X. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_X)) {
    ReadVectorData(job, da, APP_VECTOR_X, projection_data.vector_x);
  }
  // INDEX_C2M. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_INDEX_C2M), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataRawGridArray* data_real =
          dynamic_cast<application::DataRawGridArray*>(data_temp);
      projection_data.cell_index_to_matrix_index.Resize(
          PhysBAM::RANGE<TV_INT>(TV_INT(0, 0, 0),
                                 TV_INT(init_config.local_region.dx()+1,
                                        init_config.local_region.dy()+1,
                                        init_config.local_region.dz()+1)));

      data_real->LoadFromNimbus(
          &projection_data.cell_index_to_matrix_index);
      dbg(APP_LOG, "Finish reading INDEX_C2M.\n");
    }
  }
  // INDEX_M2C. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_INDEX_M2C), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataRawArrayM2C* data_real =
          dynamic_cast<application::DataRawArrayM2C*>(data_temp);
      data_real->LoadFromNimbus(
          &projection_data.matrix_index_to_cell_index);
      dbg(APP_LOG, "Finish reading INDEX_M2C.\n");
    }
  }
  // LOCAL_N. Reduction cannot take this branch.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N)) {
    ReadScalarData<int>(job, da, APP_PROJECTION_LOCAL_N,
                        projection_data.local_n);
  }
  // INTERIOR_N. Reduction cannot take this branch.
  if (data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N)) {
    ReadScalarData<int>(job, da, APP_PROJECTION_INTERIOR_N,
                        projection_data.interior_n);
  }
  // Groud III.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
    ReadScalarData<float>(job, da, APP_PROJECTION_LOCAL_TOLERANCE,
                          projection_data.local_tolerance);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_TOLERANCE)) {
    ReadScalarData<float>(job, da, APP_PROJECTION_GLOBAL_TOLERANCE,
                          projection_data.global_tolerance);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_N)) {
    ReadScalarData<int>(job, da, APP_PROJECTION_GLOBAL_N,
                        projection_data.global_n);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_DESIRED_ITERATIONS)) {
    ReadScalarData<int>(job, da, APP_PROJECTION_DESIRED_ITERATIONS,
                        projection_data.desired_iterations);
  }
  // Group IV.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_RESIDUAL)) {
    ReadScalarData<double>(job, da, APP_PROJECTION_LOCAL_RESIDUAL,
                           projection_data.local_residual);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_RHO)) {
    ReadScalarData<double>(job, da, APP_PROJECTION_LOCAL_RHO,
                           projection_data.local_rho);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_RHO)) {
    ReadScalarData<double>(job, da, APP_PROJECTION_GLOBAL_RHO,
                           projection_data.rho);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_RHO_OLD)) {
    ReadScalarData<double>(job, da, APP_PROJECTION_GLOBAL_RHO_OLD,
                           projection_data.rho_last);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA)) {
    ReadScalarData<double>(job, da, APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA,
                           projection_data.local_dot_product_for_alpha);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_ALPHA)) {
    ReadScalarData<float>(job, da, APP_PROJECTION_ALPHA, projection_data.alpha);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_BETA)) {
    ReadScalarData<float>(job, da, APP_PROJECTION_BETA, projection_data.beta);
  }
  // MATRIX_C. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::MATRIX_C)) {
    if (projection_data.matrix_a.C == NULL) {
      projection_data.matrix_a.C = new SPARSE_MATRIX_FLAT_NXN<float>;
    }
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_MATRIX_C), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataSparseMatrix* data_real =
          dynamic_cast<application::DataSparseMatrix*>(data_temp);
      data_real->LoadFromNimbus(projection_data.matrix_a.C);
      dbg(APP_LOG, "Finish reading MATRIX_A.\n");
    }
  }
  // VECTOR_Z. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_Z)) {
    ReadVectorData(job, da, APP_VECTOR_Z, projection_data.z_interior);
  }
  // VECTOR_P. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_P)) {
    ReadVectorData(job, da, APP_VECTOR_P, projection_data.p);
  }
  // VECTOR_TEMP. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    ReadVectorData(job, da, APP_VECTOR_TEMP, projection_data.temp);
  }
  assert(data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N));
  assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
  Initialize(projection_data.local_n, projection_data.interior_n);
}

template<typename TYPE_NAME> void ProjectionDriver::ReadScalarData(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, TYPE_NAME& value) {
  Data* data_temp = application::GetTheOnlyData(
      job, std::string(variable_name), da, application::READ_ACCESS);
  if (data_temp) {
    nimbus::ScalarData<TYPE_NAME>* data_real =
        dynamic_cast<nimbus::ScalarData<TYPE_NAME>*>(data_temp);
    value = data_real->scalar();
    dbg(APP_LOG, "[Data Loading]%s: %0.1f", variable_name, (float)value);
    dbg(APP_LOG, "Finish reading %s.\n", variable_name);
  }
}

void ProjectionDriver::ReadVectorData(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, VECTOR_ND<float>& value) {
  Data* data_temp = application::GetTheOnlyData(
      job, std::string(variable_name), da, application::READ_ACCESS);
  if (data_temp) {
    application::DataRawVectorNd* data_real =
        dynamic_cast<application::DataRawVectorNd*>(data_temp);
    data_real->LoadFromNimbus(&value);
    dbg(APP_LOG, "Finish reading %s.\n", variable_name);
  }
}

void ProjectionDriver::SaveToNimbus(
    const nimbus::Job* job, const nimbus::DataArray& da) {
  // VECTOR_B. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    WriteVectorData(job, da, APP_VECTOR_B, projection_data.vector_b);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_X)) {
    WriteVectorData(job, da, APP_VECTOR_X, projection_data.vector_x);
  }
  // Groud III.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
    WriteScalarData<float>(job, da, APP_PROJECTION_LOCAL_TOLERANCE,
                          projection_data.local_tolerance);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_TOLERANCE)) {
    WriteScalarData<float>(job, da, APP_PROJECTION_GLOBAL_TOLERANCE,
                          projection_data.global_tolerance);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_N)) {
    WriteScalarData<int>(job, da, APP_PROJECTION_GLOBAL_N,
                         projection_data.global_n);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_DESIRED_ITERATIONS)) {
    WriteScalarData<int>(job, da, APP_PROJECTION_DESIRED_ITERATIONS,
                        projection_data.desired_iterations);
  }
  // Group IV.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_RESIDUAL)) {
    WriteScalarData<double>(job, da, APP_PROJECTION_LOCAL_RESIDUAL,
                           projection_data.local_residual);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_RHO)) {
    WriteScalarData<double>(job, da, APP_PROJECTION_LOCAL_RHO,
                           projection_data.local_rho);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_RHO)) {
    WriteScalarData<double>(job, da, APP_PROJECTION_GLOBAL_RHO,
                           projection_data.rho);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_RHO_OLD)) {
    WriteScalarData<double>(job, da, APP_PROJECTION_GLOBAL_RHO_OLD,
                           projection_data.rho_last);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA)) {
    WriteScalarData<double>(job, da, APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA,
                           projection_data.local_dot_product_for_alpha);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_ALPHA)) {
    WriteScalarData<float>(job, da, APP_PROJECTION_ALPHA,
                           projection_data.alpha);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_BETA)) {
    WriteScalarData<float>(job, da, APP_PROJECTION_BETA, projection_data.beta);
  }
  // MATRIX_C. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::MATRIX_C)) {
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_MATRIX_C), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataSparseMatrix* data_real =
          dynamic_cast<application::DataSparseMatrix*>(data_temp);
      data_real->SaveToNimbus(*projection_data.matrix_a.C);
      dbg(APP_LOG, "Finish reading MATRIX_A.\n");
    }
  }
  // VECTOR_Z. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_Z)) {
    WriteVectorData(job, da, APP_VECTOR_Z, projection_data.z_interior);
  }
  // VECTOR_P. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_P)) {
    WriteVectorData(job, da, APP_VECTOR_P, projection_data.p);
  }
  // VECTOR_TEMP. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    WriteVectorData(job, da, APP_VECTOR_TEMP, projection_data.temp);
  }
}

template<typename TYPE_NAME> void ProjectionDriver::WriteScalarData(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, const TYPE_NAME& value) {
  Data* data_temp = application::GetTheOnlyData(
      job, std::string(variable_name), da, application::WRITE_ACCESS);
  if (data_temp) {
    nimbus::ScalarData<TYPE_NAME>* data_real =
        dynamic_cast<nimbus::ScalarData<TYPE_NAME>*>(data_temp);
    data_real->set_scalar(value);
    dbg(APP_LOG, "[Data Saving]%s: %0.1f", variable_name, (float)value);
    dbg(APP_LOG, "Finish writing %s.\n", variable_name);
  }
}

void ProjectionDriver::WriteVectorData(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, const VECTOR_ND<float>& value) {
  Data* data_temp = application::GetTheOnlyData(
      job, std::string(variable_name), da, application::WRITE_ACCESS);
  if (data_temp) {
    application::DataRawVectorNd* data_real =
        dynamic_cast<application::DataRawVectorNd*>(data_temp);
    data_real->SaveToNimbus(value);
    dbg(APP_LOG, "Finish writing %s.\n", variable_name);
  }
}

}  // namespace PhysBAM
