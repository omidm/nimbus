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
 * PhysBAM projection simulation codes are here. It serves the similar
 * functionality as WATER_DRIVER, but is seperated out with code that is not
 * related to other parts of the simulation.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 * Edited for application app_data by Chinmayee Shah <chshah@stanford.edu>
 */

#include <cassert>

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>

#include "application/water_multiple/app_data_prototypes.h"
#include "application/water_multiple/data_include.h"
#include "application/water_multiple/physbam_utils.h"
#include "data/scalar_data.h"
#include "shared/nimbus.h"

#include "application/water_multiple/projection/projection_driver.h"

namespace PhysBAM {

void ProjectionDriver::Initialize(int local_n, int interior_n) {
  partition.interior_indices.min_corner = 1;
  partition.interior_indices.max_corner = interior_n;

  // Initializes the vector if it is not transmitted.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N));
    if (projection_data.temp.Size() != local_n) {
      projection_data.temp.Resize(local_n, false);
    }
  }
  // Sets subview if necessary.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
    projection_data.temp_interior.Set_Subvector_View(
        projection_data.temp,
        partition.interior_indices);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_P_META_FORMAT)) {
    if (projection_data.meta_p.Size() == 0) {
      assert(data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N));
      projection_data.meta_p.Resize(projection_data.local_n);
    }
    assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
    projection_data.p_interior.Set_Subvector_View(
        projection_data.meta_p,
        partition.interior_indices);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
    projection_data.b_interior.Set_Subvector_View(
        projection_data.vector_b,
        partition.interior_indices);
  }
}

// Projection is broken to "smallest" code piece to allow future changes.
void ProjectionDriver::LocalInitialize() {
  if (projection_data.interior_n == 0) {
    projection_data.local_residual = 0;
    return;
  }
  SPARSE_MATRIX_FLAT_NXN<T>& A = *projection_data.matrix_a;
  // Initializes VECTOR_X.
  // VECTOR_X is only used in this job.
  projection_data.vector_x.Resize(projection_data.local_n, false);
  for (int i = 1; i <= projection_data.local_n; ++i) {
    projection_data.vector_x(i) =
        projection_data.pressure((*projection_data.matrix_index_to_cell_index)(i));
  }
  VECTOR_ND<T>& x = projection_data.vector_x;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  VECTOR_ND<T>& temp = projection_data.temp;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  // Initializes vector B and local residual.
  A.thread_queue = thread_queue;
  A.Times(x, temp);
  A.thread_queue = NULL;
  b_interior -= temp_interior;
  projection_data.local_residual = b_interior.Max_Abs();
  // Calculate matrix C.
  assert(A.C != NULL);
  A.Nimbus_Create_Submatrix(partition.interior_indices, A.C);
  A.C->In_Place_Incomplete_Cholesky_Factorization(
      pcg.modified_incomplete_cholesky,
      pcg.modified_incomplete_cholesky_coefficient,
      pcg.preconditioner_zero_tolerance,
      pcg.preconditioner_zero_replacement);
  // Initializes vector pressure.
  projection_data.vector_pressure.Resize(projection_data.interior_n, false);
  for (int i = 1; i <= projection_data.interior_n; ++i) {
    projection_data.vector_pressure(i) =
        projection_data.pressure((*projection_data.matrix_index_to_cell_index)(i));
  }
  // Initializes other vectors.
  projection_data.temp.Resize(projection_data.local_n, false);
  projection_data.temp.Fill(0);
  projection_data.z_interior.Resize(projection_data.interior_n, false);
  projection_data.z_interior.Fill(0);
  projection_data.meta_p.Resize(projection_data.local_n, false);
  projection_data.meta_p.Fill(0);
}

void ProjectionDriver::GlobalInitialize() {
  projection_data.global_n = Global_Sum(projection_data.interior_n);
  projection_data.global_tolerance =
      Global_Max(projection_data.local_tolerance);

  projection_data.desired_iterations = projection_data.global_n;
  if (pcg.maximum_iterations) {
    projection_data.desired_iterations =
        min(projection_data.desired_iterations, pcg.maximum_iterations);
  }
}

void ProjectionDriver::TransformPressureResult() {
  for (int i = 1; i <= projection_data.interior_n; ++i) {
    projection_data.pressure((*projection_data.matrix_index_to_cell_index)(i))
        = projection_data.vector_pressure(i);
  }
}

// Step one starts.
void ProjectionDriver::DoPrecondition() {
  if (projection_data.interior_n == 0) {
    return;
  }
  SPARSE_MATRIX_FLAT_NXN<T>& A = *projection_data.matrix_a;
  VECTOR_ND<T>& z_interior = projection_data.z_interior;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  // Time consuming part.
  // Multi-threaded introduced.
  // A.C->thread_queue = thread_queue;
  A.C->Solve_Forward_Substitution(b_interior, temp_interior, true);
  A.C->Solve_Backward_Substitution(temp_interior, z_interior, false, true);
  // A.C->thread_queue = NULL;
}

void ProjectionDriver::CalculateLocalRho() {
  if (projection_data.interior_n == 0) {
    projection_data.local_rho = 0;
    return;
  }
  VECTOR_ND<T>& z_interior = projection_data.z_interior;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  projection_data.local_rho =
      VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior, b_interior);
}
// Step one finishes.

void ProjectionDriver::ReduceRho() {
  projection_data.rho_last = projection_data.rho;
  projection_data.rho = Global_Sum(projection_data.local_rho);
  if (projection_data.iteration == 1) {
    projection_data.beta = 0;
  } else {
    projection_data.beta = (T)(projection_data.rho / projection_data.rho_last);
  }
}

// Step two starts.
void ProjectionDriver::UpdateSearchVector() {
  if (projection_data.interior_n == 0) {
    return;
  }
  int interior_n = partition.interior_indices.Size() + 1;
  VECTOR_ND<T>& z_interior = projection_data.z_interior;
  VECTOR_ND<T>& p_interior = projection_data.p_interior;
  // Search vector p is updated here.
  if (projection_data.iteration == 1) {
    p_interior = z_interior;
  } else {
    for (int i = 1; i <= interior_n; i++)
      p_interior(i) = z_interior(i) + projection_data.beta * p_interior(i);
  }
}
// Step two finishes.

// Search-vector p is exchanged between the two steps.

// Step three starts.
void ProjectionDriver::UpdateTempVector() {
  if (projection_data.interior_n == 0) {
    return;
  }
  SPARSE_MATRIX_FLAT_NXN<T>& A = *projection_data.matrix_a;
  VECTOR_ND<T>& temp = projection_data.temp;
  VECTOR_ND<T>& p = projection_data.meta_p;
  // Search vector p is used here.
  // Time consuming part.
  A.thread_queue = thread_queue;
  A.Times(p, temp);
  A.thread_queue = NULL;
}

void ProjectionDriver::CalculateLocalAlpha() {
  if (projection_data.interior_n == 0) {
    projection_data.local_dot_product_for_alpha = 0;
    return;
  }
  VECTOR_ND<T>& p_interior = projection_data.p_interior;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  // Search vector p is used here.
  projection_data.local_dot_product_for_alpha =
      VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior, temp_interior);
}
// Step three finishes.

void ProjectionDriver::ReduceAlpha() {
  projection_data.alpha =
      (T) (projection_data.rho /
           Global_Sum(projection_data.local_dot_product_for_alpha));
}

// Step four starts.
void ProjectionDriver::UpdateOtherVectors() {
  if (projection_data.interior_n == 0) {
    return;
  }
  int interior_n = partition.interior_indices.Size()+1;
  VECTOR_ND<T>& p_interior = projection_data.p_interior;
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  VECTOR_ND<T>& temp_interior = projection_data.temp_interior;
  for (int i = 1; i <= interior_n; i++) {
    b_interior(i) -= projection_data.alpha * temp_interior(i);
  }
  for (int i = 1; i <= interior_n; i++) {
    projection_data.vector_pressure(i) += projection_data.alpha * p_interior(i);
  }
}

void ProjectionDriver::CalculateLocalResidual() {
  if (projection_data.interior_n == 0) {
    projection_data.local_residual = 0;
    return;
  }
  VECTOR_ND<T>& b_interior = projection_data.b_interior;
  projection_data.local_residual = b_interior.Max_Abs();
}

bool ProjectionDriver::DecideToSpawnNextIteration() {
  if (print_debug) dbg(APP_LOG, "[CONTROL_FLOW]");
  if (projection_data.residual <= projection_data.global_tolerance) {
    return false;
  }
  if (projection_data.iteration == projection_data.desired_iterations) {
    return false;
  }
  return true;
}
// Step four finishes.

void ProjectionDriver::LoadFromNimbus(
    const nimbus::Job* job, const nimbus::DataArray& da) {
  application::ScopeTimer scope_timer("projection_load");
  AppData_LoadFromNimbus(job, da);
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
    if (print_debug) dbg(APP_LOG, "Read %s=%0.9f.\n", variable_name, (float)value);
  } else {
    if (print_debug) dbg(APP_LOG, "Variable %s uninitialized.\n", variable_name);
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
    if (print_debug) dbg(APP_LOG, "Finish reading %s.\n", variable_name);
  } else {
    if (print_debug) dbg(APP_LOG, "Flag is set but data is not local:%s.\n", variable_name);
  }
}

void ProjectionDriver::SaveToNimbus(
    const nimbus::Job* job, const nimbus::DataArray& da) {
  application::ScopeTimer scope_timer("projection_save");
  AppData_SaveToNimbus(job, da);
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
    if (print_debug) dbg(APP_LOG, "Write %s=%0.9f.\n", variable_name, (float)value);
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
    if (print_debug) dbg(APP_LOG, "Finish writing %s.\n", variable_name);
  }
}

template void ProjectionDriver::WriteScalarData<float>(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, const float& value);
template void ProjectionDriver::WriteScalarData<double>(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, const double& value);
template void ProjectionDriver::WriteScalarData<int>(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, const int& value);

template void ProjectionDriver::ReadScalarData<float>(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, float& value);
template void ProjectionDriver::ReadScalarData<double>(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, double& value);
template void ProjectionDriver::ReadScalarData<int>(
    const nimbus::Job* job, const nimbus::DataArray& da,
    const char* variable_name, int& value);
}  // namespace PhysBAM
