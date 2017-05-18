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
 * Author: Hang Qu<quhang@stanford.edu>
 */

#include <cassert>

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>

#include "applications/physbam/smoke/app_data_prototypes.h"
#include "applications/physbam/smoke/data_include.h"
#include "applications/physbam/smoke/physbam_utils.h"
#include "src/data/scalar_data.h"
#include "src/shared/nimbus.h"

#include "applications/physbam/smoke/projection/projection_driver.h"

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
  A.Times(x, temp);
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
  A.C->Solve_Forward_Substitution(b_interior, temp_interior, true);
  A.C->Solve_Backward_Substitution(temp_interior, z_interior, false, true);
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
  A.Times(p, temp);
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
  dbg(APP_LOG, "[CONTROL_FLOW]");
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
  if (true) {
    AppData_LoadFromNimbus(job, da);
    return;
  }
  nimbus::int_dimension_t array_shift[3] = {
    init_config.local_region.x() - 1,
    init_config.local_region.y() - 1,
    init_config.local_region.z() - 1};
  nimbus::PdiVector pdv;
  GeometricRegion array_reg_central(init_config.local_region.x(),
                                    init_config.local_region.y(),
                                    init_config.local_region.z(),
                                    init_config.local_region.dx(),
                                    init_config.local_region.dy(),
                                    init_config.local_region.dz());
  GeometricRegion array_reg_thin_outer(init_config.local_region.x()-1,
                                       init_config.local_region.y()-1,
                                       init_config.local_region.z()-1,
                                       init_config.local_region.dx()+2,
                                       init_config.local_region.dy()+2,
                                       init_config.local_region.dz()+2);
  GRID<TV> grid;
  grid.Initialize(TV_INT(init_config.local_region.dx(),
                         init_config.local_region.dy(),
                         init_config.local_region.dz()),
                  application::GridToRange(init_config.global_region,
                                           init_config.local_region));

  Log log_timer;

  // LOCAL_N.
  // Reduction on LOCAL_N is never used and thus not supported.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N)) {
    ReadScalarData<int>(job, da, APP_PROJECTION_LOCAL_N,
                        projection_data.local_n);
  }
  // INTERIOR_N. Reducible.
  if (data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N)) {
    if (application::GetTranslatorData(
            job, std::string(APP_PROJECTION_INTERIOR_N), da, &pdv,
            application::READ_ACCESS)) {
      dbg(APP_LOG, "Reducing PROJECTION_INTERIOR_N sum(");
      projection_data.interior_n = 0;
      nimbus::PdiVector::const_iterator iter = pdv.begin();
      for (; iter != pdv.end(); ++iter) {
        const nimbus::PhysicalDataInstance* instance = *iter;
        nimbus::ScalarData<int>* data_real =
            dynamic_cast<nimbus::ScalarData<int>*>(instance->data());
        int value = data_real->scalar();
        dbg(APP_LOG, "%d ", value);
        projection_data.interior_n += value;
      }
      dbg(APP_LOG, ") = %d.\n", projection_data.interior_n);
    } else {
      dbg(APP_LOG, "PROJECTION_INTERIOR_N flag"
          " is set but data is not local.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }

  log_timer.StartTimer();
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    projection_data.pressure.Resize(grid.Domain_Indices(1));
    if (application::GetTranslatorData(job, std::string(APP_PRESSURE), da, &pdv,
                                       application::READ_ACCESS)) {
      translator.ReadScalarArrayFloat(
          &array_reg_thin_outer, array_shift, &pdv, &projection_data.pressure);
      dbg(APP_LOG, "Finish reading PRESSURE.\n");
    } else {
      dbg(APP_LOG, "PRESSURE flag is set but data is not local.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, pressure time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // MATRIX_A. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::MATRIX_A)) {
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_MATRIX_A), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataSparseMatrix* data_real =
          dynamic_cast<application::DataSparseMatrix*>(data_temp);
      // The memory will be allocated automatically.
      projection_data.matrix_a = new SPARSE_MATRIX_FLAT_NXN<T>;
      data_real->LoadFromNimbus(projection_data.matrix_a);
      dbg(APP_LOG, "Finish reading MATRIX_A.\n");
    } else {
      dbg(APP_LOG, "MATRIX_A flag is set but data is not local.\n");
    }
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, matrix_a time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // VECTOR_B. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    ReadVectorData(job, da, APP_VECTOR_B, projection_data.vector_b);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_b time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // INDEX_C2M. It cannot be splitted or merged.
  // index_c2m.
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

      data_real->LoadFromNimbus(&projection_data.cell_index_to_matrix_index);
      dbg(APP_LOG, "Finish reading INDEX_C2M.\n");
    } else {
      dbg(APP_LOG, "INDEX_C2M flag is set but data is not local.\n");
    }
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, index_c2m time:%f.\n",
      log_timer.timer());

  log_timer.StartTimer();
  if (data_config.GetFlag(DataConfig::VECTOR_P_META_FORMAT)) {
    assert(data_config.GetFlag(DataConfig::INDEX_C2M));
    if (application::GetTranslatorData(job,
        std::string(APP_VECTOR_P_META_FORMAT), da, &pdv, application::READ_ACCESS)) {
      translator.ReadCompressedScalarArray<float>(
          array_reg_thin_outer, array_shift, pdv,
          &projection_data.meta_p,
          projection_data.local_n,
          projection_data.cell_index_to_matrix_index);
      dbg(APP_LOG, "Finish reading VECTOR_P_META_FORMAT.\n");
    } else {
      dbg(APP_LOG, "VECTOR_P_META_FORMAT flag is set but data is not local.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_p_meta_format time:%f.\n",
      log_timer.timer());

  log_timer.StartTimer();
  // INDEX_M2C. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_INDEX_M2C), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataRawArrayM2C* data_real =
          dynamic_cast<application::DataRawArrayM2C*>(data_temp);
      // The memory will be allocated automatically.
      projection_data.matrix_index_to_cell_index = new ARRAY<TV_INT>;
      data_real->LoadFromNimbus(projection_data.matrix_index_to_cell_index);
      dbg(APP_LOG, "Finish reading INDEX_M2C.\n");
    } else {
      dbg(APP_LOG, "INDEX_M2C flag is set but data is not local.\n");
    }
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, index_m2c time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // Group III.
  // PROJECTION_LOCAL_TOLERANCE. Reducible.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
    if (application::GetTranslatorData(
            job, std::string(APP_PROJECTION_LOCAL_TOLERANCE), da, &pdv,
            application::READ_ACCESS)) {
      dbg(APP_LOG, "Reducing PROJECTION_LOCAL_TOLERANCE MAX(");
      projection_data.local_tolerance = 0;
      nimbus::PdiVector::const_iterator iter = pdv.begin();
      for (; iter != pdv.end(); ++iter) {
        const nimbus::PhysicalDataInstance* instance = *iter;
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(instance->data());
        float value = data_real->scalar();
        dbg(APP_LOG, "%f ", value);
        if (value > projection_data.local_tolerance) {
          projection_data.local_tolerance = value;
        }
      }
      dbg(APP_LOG, ") = %f.\n", projection_data.local_tolerance);
    } else {
      dbg(APP_LOG, "PROJECTION_LOCAL_TOLERANCE flag"
          " is set but data is not local.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
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
  // PROJECTION_LOCAL_RESIDUAL. Reducible.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_RESIDUAL)) {
    if (application::GetTranslatorData(
            job, std::string(APP_PROJECTION_LOCAL_RESIDUAL), da, &pdv,
            application::READ_ACCESS)) {
      dbg(APP_LOG, "Reducing PROJECTION_LOCAL_RESIDUAL max(\n");
      projection_data.local_residual = 0;
      nimbus::PdiVector::const_iterator iter = pdv.begin();
      for (; iter != pdv.end(); ++iter) {
        const nimbus::PhysicalDataInstance* instance = *iter;
        nimbus::ScalarData<double>* data_real =
            dynamic_cast<nimbus::ScalarData<double>*>(instance->data());
        double value = data_real->scalar();
        dbg(APP_LOG, "%f ", value);
        if (value > projection_data.local_residual) {
          projection_data.local_residual = value;
        }
      }
      dbg(APP_LOG, ") = %f.\n", projection_data.local_residual);
    } else {
      dbg(APP_LOG, "PROJECTION_LOCAL_RESIDUAL flag"
          "is set but data is not local.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }
  // PROJECTION_LOCAL_RHO. Reducible.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_RHO)) {
    if (application::GetTranslatorData(
            job, std::string(APP_PROJECTION_LOCAL_RHO), da, &pdv,
            application::READ_ACCESS)) {
      dbg(APP_LOG, "Reducing PROJECTION_LOCAL_RHO sum(:");
      projection_data.local_rho = 0;
      nimbus::PdiVector::const_iterator iter = pdv.begin();
      for (; iter != pdv.end(); ++iter) {
        const nimbus::PhysicalDataInstance* instance = *iter;
        nimbus::ScalarData<double>* data_real =
            dynamic_cast<nimbus::ScalarData<double>*>(instance->data());
        double value = data_real->scalar();
        dbg(APP_LOG, "%f ", value);
        projection_data.local_rho += value;
      }
      dbg(APP_LOG, ") = %f.\n", projection_data.local_rho);
    } else {
      dbg(APP_LOG, "PROJECTION_LOCAL_RHO flag"
          "is set but data is not local.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_RHO)) {
    ReadScalarData<double>(job, da, APP_PROJECTION_GLOBAL_RHO,
                           projection_data.rho);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_GLOBAL_RHO_OLD)) {
    ReadScalarData<double>(job, da, APP_PROJECTION_GLOBAL_RHO_OLD,
                           projection_data.rho_last);
  }
  // PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA. Reducible.
  if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA)) {
    if (application::GetTranslatorData(
            job, std::string(APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA),
            da, &pdv, application::READ_ACCESS)) {
      dbg(APP_LOG, "Reducing PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA sum(:");
      projection_data.local_dot_product_for_alpha = 0;
      nimbus::PdiVector::const_iterator iter = pdv.begin();
      for (; iter != pdv.end(); ++iter) {
        const nimbus::PhysicalDataInstance* instance = *iter;
        nimbus::ScalarData<double>* data_real =
            dynamic_cast<nimbus::ScalarData<double>*>(instance->data());
        double value = data_real->scalar();
        dbg(APP_LOG, "%f ", value);
        projection_data.local_dot_product_for_alpha += value;
      }
      dbg(APP_LOG, ") = %f.\n", projection_data.local_dot_product_for_alpha);
    } else {
      dbg(APP_LOG, "PROJECTION_LOCAL_PRODUCT_FOR_ALPHA flag"
          "is set but data is not local.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_ALPHA)) {
    ReadScalarData<float>(job, da, APP_PROJECTION_ALPHA, projection_data.alpha);
  }
  if (data_config.GetFlag(DataConfig::PROJECTION_BETA)) {
    ReadScalarData<float>(job, da, APP_PROJECTION_BETA, projection_data.beta);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, scalar time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // MATRIX_C. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::MATRIX_C)) {
    if (projection_data.matrix_a == NULL) {
      projection_data.matrix_a = new SPARSE_MATRIX_FLAT_NXN<float>;
    }
    if (projection_data.matrix_a->C == NULL) {
      projection_data.matrix_a->C = new SPARSE_MATRIX_FLAT_NXN<float>;
    }
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_MATRIX_C), da, application::READ_ACCESS);
    if (data_temp) {
      application::DataSparseMatrix* data_real =
          dynamic_cast<application::DataSparseMatrix*>(data_temp);
      // Memory allocation happens inside the call.
      data_real->LoadFromNimbus(projection_data.matrix_a->C);
      dbg(APP_LOG, "Finish reading MATRIX_C.\n");
    } else {
      dbg(APP_LOG, "MATRIX_C flag"
          "is set but data is not local.\n");
    }
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, matrix_c time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // VECTOR_Z. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_Z)) {
    ReadVectorData(job, da, APP_VECTOR_Z, projection_data.z_interior);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_z time:%f.\n",
      log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_TEMP. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    ReadVectorData(job, da, APP_VECTOR_TEMP, projection_data.temp);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_PRESSURE)) {
    ReadVectorData(job, da, APP_VECTOR_PRESSURE,
                   projection_data.vector_pressure);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_temp time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  Initialize(projection_data.local_n, projection_data.interior_n);
  dbg(APP_LOG, "[PROJECTION] LOAD, else time:%f.\n",
      log_timer.timer());
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
    dbg(APP_LOG, "[Data Loading]%s: %0.9f\n", variable_name, (float)value);
    dbg(APP_LOG, "Finish reading %s.\n", variable_name);
  } else {
    dbg(APP_LOG, "Flag is set but data is not local:%s.\n", variable_name);
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
  } else {
    dbg(APP_LOG, "Flag is set but data is not local:%s.\n", variable_name);
  }
}

void ProjectionDriver::SaveToNimbus(
    const nimbus::Job* job, const nimbus::DataArray& da) {
  application::ScopeTimer scope_timer("projection_save");
  if (true) {
    AppData_SaveToNimbus(job, da);
    return;
  }
  nimbus::int_dimension_t array_shift[3] = {
      init_config.local_region.x() - 1,
      init_config.local_region.y() - 1,
      init_config.local_region.z() - 1};
  nimbus::PdiVector pdv;
  GeometricRegion array_reg_central(init_config.local_region.x(),
                                    init_config.local_region.y(),
                                    init_config.local_region.z(),
                                    init_config.local_region.dx(),
                                    init_config.local_region.dy(),
                                    init_config.local_region.dz());
  GeometricRegion array_reg_thin_outer(init_config.local_region.x()-1,
                                       init_config.local_region.y()-1,
                                       init_config.local_region.z()-1,
                                       init_config.local_region.dx()+2,
                                       init_config.local_region.dy()+2,
                                       init_config.local_region.dz()+2);

  Log log_timer;
  log_timer.StartTimer();

  const std::string pressure_string = std::string(APP_PRESSURE);
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    if (application::GetTranslatorData(job, pressure_string, da, &pdv,
                                       application::WRITE_ACCESS)) {
      translator.WriteScalarArrayFloat(
          &array_reg_central, array_shift, &pdv, &projection_data.pressure);
      dbg(APP_LOG, "Finish writing pressure.\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, pressure time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // VECTOR_B. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    WriteVectorData(job, da, APP_VECTOR_B, projection_data.vector_b);
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, vector_b time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
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
  dbg(APP_LOG, "[PROJECTION] SAVE, scalar time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // MATRIX_C. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::MATRIX_C)) {
    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_MATRIX_C), da, application::WRITE_ACCESS);
    if (data_temp) {
      application::DataSparseMatrix* data_real =
          dynamic_cast<application::DataSparseMatrix*>(data_temp);
      data_real->SaveToNimbus(*projection_data.matrix_a->C);
      dbg(APP_LOG, "Finish writing MATRIX_C.\n");
    }
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, matrix_c time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();
  // VECTOR_Z. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_Z)) {
    WriteVectorData(job, da, APP_VECTOR_Z, projection_data.z_interior);
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, vector_z time:%f.\n",
      log_timer.timer());
  log_timer.StartTimer();

  log_timer.StartTimer();
  if (data_config.GetFlag(DataConfig::VECTOR_P_META_FORMAT)) {
    assert(data_config.GetFlag(DataConfig::INDEX_C2M));
    if (application::GetTranslatorData(job,
        std::string(APP_VECTOR_P_META_FORMAT), da, &pdv, application::WRITE_ACCESS)) {
      translator.WriteCompressedScalarArray<float>(
          array_reg_thin_outer, array_shift, pdv,
          projection_data.meta_p, projection_data.local_n,
          projection_data.cell_index_to_matrix_index);
      dbg(APP_LOG, "Finish reading VECTOR_P_META_FORMAT).\n");
    }
    application::DestroyTranslatorObjects(&pdv);
  }
  dbg(APP_LOG, "[PROJECTION] Save, vector_p_meta_format time:%f.\n",
      log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_TEMP. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    WriteVectorData(job, da, APP_VECTOR_TEMP, projection_data.temp);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_PRESSURE)) {
    WriteVectorData(job, da, APP_VECTOR_PRESSURE,
                    projection_data.vector_pressure);
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, vector_temp time:%f.\n",
      log_timer.timer());
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
    dbg(APP_LOG, "[Data Saving]%s: %0.9f\n", variable_name, (float)value);
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
