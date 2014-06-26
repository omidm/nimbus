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
 * All the cached version loading and saving function for projection loop
 * iteration.
 * Author: Hang Qu<quhang@stanford.edu>
 */

#include <cassert>

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>

#include "application/water_multiple/cache_prototypes.h"
#include "application/water_multiple/data_include.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/app_utils.h"
#include "data/scalar_data.h"
#include "shared/nimbus.h"

#include "application/water_multiple/projection/projection_driver.h"

namespace PhysBAM {

void ProjectionDriver::Cache_Initialize(int local_n, int interior_n) {
  partition.interior_indices.min_corner = 1;
  partition.interior_indices.max_corner = interior_n;

  // Initializes the vector if it is not transmitted.
  if (projection_data.temp.Size() == 0 &&
      data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N));
    projection_data.temp.Resize(local_n, false);
  }
  if (projection_data.p.Size() == 0 &&
      data_config.GetFlag(DataConfig::VECTOR_P)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N));
    projection_data.p.Resize(local_n, false);
    assert(data_config.GetFlag(DataConfig::INDEX_M2C));
    assert(data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N));
    // TODO(addcache): should be removed.
    // Translates from grid-format vector p to vector p, based on index.
    for (int i = 1; i <= local_n; ++i) {
      projection_data.p(i) =
          projection_data.grid_format_vector_p(
              (*projection_data.matrix_index_to_cell_index)(i));
    }
  }
  if (projection_data.z_interior.Size() == 0 &&
      data_config.GetFlag(DataConfig::VECTOR_Z)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
    projection_data.z_interior.Resize(interior_n, false);
  }

  // Sets subview if necessary.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
    projection_data.temp_interior.Set_Subvector_View(
        projection_data.temp,
        partition.interior_indices);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_P)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
    projection_data.p_interior.Set_Subvector_View(
        projection_data.p,
        partition.interior_indices);
  }
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    assert(data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N));
    projection_data.b_interior.Set_Subvector_View(
        projection_data.vector_b,
        partition.interior_indices);
  }
}

void ProjectionDriver::Cache_LoadFromNimbus(
    const nimbus::Job* job, const nimbus::DataArray& da) {
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
  nimbus::CacheManager *cm = job->GetCacheManager();
  Log log_timer;

  log_timer.StartTimer();
  // pressure.
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    nimbus::DataArray read, write;
    const std::string pressure_string = std::string(APP_PRESSURE);
    application::GetReadData(*job, pressure_string, da, &read);
    application::GetWriteData(*job, pressure_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_thin_outer,
            write, array_reg_thin_outer,
            application::kCachePressure, array_reg_thin_outer,
            nimbus::cache::EXCLUSIVE);
    cache_pressure = dynamic_cast<application::CacheScalarArray<T>*>(cache_var);
    assert(cache_pressure != NULL);
    typedef typename PhysBAM::ARRAY<T, TV_INT> T_SCALAR_ARRAY;
    T_SCALAR_ARRAY* pressure = cache_pressure->data();
    T_SCALAR_ARRAY::Exchange_Arrays(*pressure, projection_data.pressure);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, pressure time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // MATRIX_A.
  if (data_config.GetFlag(DataConfig::MATRIX_A)) {
    nimbus::DataArray read, write;
    const std::string matrix_a_string = std::string(APP_MATRIX_A);
    application::GetReadData(*job, matrix_a_string, da, &read);
    application::GetWriteData(*job, matrix_a_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_central,
            write, array_reg_central,
            application::kCacheSparseMatrixA, array_reg_central,
            nimbus::cache::EXCLUSIVE);
    cache_matrix_a = dynamic_cast<application::CacheSparseMatrix*>(cache_var);
    assert(cache_matrix_a != NULL);
    assert(projection_data.matrix_a == NULL);
    projection_data.matrix_a = cache_matrix_a->data();
    assert(projection_data.matrix_a != NULL);
    projection_data.matrix_a->C = NULL;
  } else {
    assert(projection_data.matrix_a == NULL);
    projection_data.matrix_a = new SPARSE_MATRIX_FLAT_NXN<T>;
    projection_data.matrix_a->C = NULL;
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, matrix_a time:%f.\n", log_timer.timer());

  // MATRIX_C.
  log_timer.StartTimer();
  if (data_config.GetFlag(DataConfig::MATRIX_C)) {
    assert(projection_data.matrix_a != NULL);
    assert(projection_data.matrix_a->C == NULL);
    nimbus::DataArray read, write;
    const std::string matrix_c_string = std::string(APP_MATRIX_C);
    application::GetReadData(*job, matrix_c_string, da, &read);
    application::GetWriteData(*job, matrix_c_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_central,
            write, array_reg_central,
            application::kCacheSparseMatrixC, array_reg_central,
            nimbus::cache::EXCLUSIVE);
    cache_matrix_c = dynamic_cast<application::CacheSparseMatrix*>(cache_var);
    assert(cache_matrix_c != NULL);
    projection_data.matrix_a->C = cache_matrix_c->data();
    assert(projection_data.matrix_a->C != NULL);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, matrix_c time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // INDEX_M2C. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
    nimbus::DataArray read, write;
    const std::string index_m2c_string = std::string(APP_INDEX_M2C);
    application::GetReadData(*job, index_m2c_string, da, &read);
    application::GetWriteData(*job, index_m2c_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_central,
            write, array_reg_central,
            application::kCacheArrayM2C, array_reg_central,
            nimbus::cache::EXCLUSIVE);
    cache_index_m2c = dynamic_cast<application::CacheArrayM2C*>(cache_var);
    assert(cache_index_m2c != NULL);
    projection_data.matrix_index_to_cell_index = cache_index_m2c->data();
    assert(projection_data.matrix_index_to_cell_index != NULL);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, index_m2c time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_B. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    ReadVectorData(job, da, APP_VECTOR_B, projection_data.vector_b);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_b time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
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

      data_real->LoadFromNimbus(&projection_data.cell_index_to_matrix_index);
      dbg(APP_LOG, "Finish reading INDEX_C2M.\n");
    } else {
      dbg(APP_LOG, "INDEX_C2M flag is set but data is not local.\n");
    }
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, index_c2m time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
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
  dbg(APP_LOG, "[PROJECTION] LOAD, scalar time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_Z. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_Z)) {
    ReadVectorData(job, da, APP_VECTOR_Z, projection_data.z_interior);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_z time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_P. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_P)) {
    nimbus::DataArray read, write;
    const std::string vector_p_string = std::string(APP_VECTOR_P);
    application::GetReadData(*job, vector_p_string, da, &read);
    application::GetWriteData(*job, vector_p_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_thin_outer,
            write, array_reg_central,
            application::kCacheVectorP, array_reg_thin_outer,
            nimbus::cache::EXCLUSIVE);
    cache_vector_p = dynamic_cast<application::CacheScalarArray<T>*>(cache_var);
    assert(cache_vector_p != NULL);
    typedef typename PhysBAM::ARRAY<T, TV_INT> T_SCALAR_ARRAY;
    T_SCALAR_ARRAY* vector_p = cache_vector_p->data();
    T_SCALAR_ARRAY::Exchange_Arrays(*vector_p,
                                    projection_data.grid_format_vector_p);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_p time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_TEMP. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    ReadVectorData(job, da, APP_VECTOR_TEMP, projection_data.temp);
  }
  dbg(APP_LOG, "[PROJECTION] LOAD, vector_temp time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  Initialize(projection_data.local_n, projection_data.interior_n);
  dbg(APP_LOG, "[PROJECTION] LOAD, else time:%f.\n", log_timer.timer());
}

void ProjectionDriver::Cache_SaveToNimbus(
    const nimbus::Job* job, const nimbus::DataArray& da) {
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

  nimbus::CacheManager *cm = job->GetCacheManager();

  Log log_timer;

  log_timer.StartTimer();
  if (cache_pressure) {
    typedef typename PhysBAM::ARRAY<T, TV_INT> T_SCALAR_ARRAY;
    T_SCALAR_ARRAY* pressure = cache_pressure->data();
    T_SCALAR_ARRAY::Exchange_Arrays(*pressure, projection_data.pressure);
    cm->ReleaseAccess(cache_pressure);
    cache_pressure = NULL;
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, pressure time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // MATRIX_C. It cannot be splitted or merged.
  if (cache_matrix_c) {
    cm->ReleaseAccess(cache_matrix_c);
    cache_matrix_c = NULL;
    projection_data.matrix_a->C = NULL;
  } else {
    assert(projection_data.matrix_a->C == NULL);
  }
  // MATRIX_A has to be deleted after MATRIX_C.
  if (cache_matrix_a) {
    cm->ReleaseAccess(cache_matrix_a);
    cache_matrix_a = NULL;
    projection_data.matrix_a = NULL;
  } else {
    assert(projection_data.matrix_a);
    delete projection_data.matrix_a;
    projection_data.matrix_a = NULL;
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, matrixes time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_B. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    WriteVectorData(job, da, APP_VECTOR_B, projection_data.vector_b);
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, vector_b time:%f.\n", log_timer.timer());

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
  dbg(APP_LOG, "[PROJECTION] SAVE, scalar time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_Z. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_Z)) {
    WriteVectorData(job, da, APP_VECTOR_Z, projection_data.z_interior);
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, vector_z time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_P.
  if (cache_vector_p) {
    typedef typename PhysBAM::ARRAY<T, TV_INT> T_SCALAR_ARRAY;
    T_SCALAR_ARRAY* vector_p = cache_vector_p->data();
    // TODO(addcache) to be replaced.
    for (int i = 1; i <= projection_data.local_n; ++i) {
      projection_data.grid_format_vector_p(
          (*projection_data.matrix_index_to_cell_index)(i)) =
          projection_data.p(i);
    }
    T_SCALAR_ARRAY::Exchange_Arrays(*vector_p,
                                    projection_data.grid_format_vector_p);
    cm->ReleaseAccess(cache_vector_p);
    cache_vector_p = NULL;
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, vector_p time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  // VECTOR_TEMP. It cannot be splitted or merged.
  if (data_config.GetFlag(DataConfig::VECTOR_TEMP)) {
    WriteVectorData(job, da, APP_VECTOR_TEMP, projection_data.temp);
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, vector_temp time:%f.\n", log_timer.timer());

  log_timer.StartTimer();
  if (cache_index_m2c) {
    cm->ReleaseAccess(cache_index_m2c);
    cache_index_m2c = NULL;
    projection_data.matrix_index_to_cell_index = NULL;
  }
  dbg(APP_LOG, "[PROJECTION] SAVE, index_m2c time:%f.\n", log_timer.timer());
}

}  // namespace PhysBAM
