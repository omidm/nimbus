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
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include "applications/physbam/water//app_data_prototypes.h"
#include "applications/physbam/water//parameters.h"

namespace application {


// face velocities
AppDataFaceArray<T> kAppDataFaceVel;
AppDataFaceArray<T> kAppDataFaceVelGhost;

// psi
AppDataFaceArray<bool> kAppDataPsiN;
AppDataScalarArray<bool> kAppDataPsiD;

//density
AppDataScalarArray<T> kAppDataDensity;
AppDataScalarArray<T> kAppDataDensityGhost;

// phi
AppDataScalarArray<T> kAppDataPhi3;
AppDataScalarArray<T> kAppDataPhi7;
AppDataScalarArray<T> kAppDataPhi8;

// particle levelset evolution
// AppDataParticleLevelsetEvolution<float> kAppDataPLE;

// pressure
AppDataScalarArray<T> kAppDataPressure;

// Varibales for projection.
AppDataScalarArray<int> kAppDataColors;
AppDataScalarArray<T> kAppDataDivergence;

// TODO(quhang): this app_data variable is questionable, because it cannot be
// deleted if meta_p is being used.
AppDataRawGridArray kAppDataIndexC2M;
AppDataArrayM2C kAppDataArrayM2C;
AppDataCompressedScalarArray<float> kAppDataMetaP;

AppDataSparseMatrix kAppDataSparseMatrixA;
AppDataSparseMatrix kAppDataSparseMatrixC;


AppDataVector kAppDataVectorB;
AppDataVector kAppDataVectorPressure;
AppDataVector kAppDataVectorZ;
AppDataVector kAppDataVectorTemp;


void InitializeAppDataPrototypes(nimbus::GeometricRegion kRegion) {
  // face velocities
  kAppDataFaceVel = AppDataFaceArray<T>(kRegion, 0, true, "face_vel");
  kAppDataFaceVelGhost = AppDataFaceArray<T>(kRegion, 3, true, "face_vel_ghost");
  
  kAppDataDensity = AppDataScalarArray<T>(kRegion, 0, true, "density");
  kAppDataDensityGhost = AppDataScalarArray<T>(kRegion, 3, true, "density_ghost");

  // psi
  kAppDataPsiN = AppDataFaceArray<bool>(kRegion, 1, true, "psi_n");
  kAppDataPsiD = AppDataScalarArray<bool>(kRegion, 1, true, "psi_d");
  
  // phi
  kAppDataPhi3 = AppDataScalarArray<T>(kRegion, 3, true, "phi_3");
  kAppDataPhi7 = AppDataScalarArray<T>(kRegion, 7, true, "phi_7");
  kAppDataPhi8 = AppDataScalarArray<T>(kRegion, 8, true, "phi_8");
  
  // particle levelset evolution
  // kAppDataPLE = AppDataParticleLevelsetEvolution<float>(kRegion, 3, true, "particle_container");
  
  // pressure
  kAppDataPressure = AppDataScalarArray<T>(kRegion, 1, true, "pressure");
  
  // Varibales for projection.
  kAppDataColors = AppDataScalarArray<int>(kRegion, 1, true, "filled_region_colors");
  kAppDataDivergence = AppDataScalarArray<T>(kRegion, 1, true, "divergence");
  
  // TODO(quhang): this app_data variable is questionable, because it cannot be
  // deleted if meta_p is being used.
  kAppDataIndexC2M = AppDataRawGridArray(kRegion, true, "index_c2m");
  kAppDataArrayM2C = AppDataArrayM2C(kRegion, true, "index_m2c");
  kAppDataMetaP = AppDataCompressedScalarArray<float>(kRegion, 1, true, "vector_p_meta_format");
  
  kAppDataSparseMatrixA = AppDataSparseMatrix(kRegion, true, "matrix_c");
  kAppDataSparseMatrixC = AppDataSparseMatrix(kRegion, true, "matrix_a");
  
  
  kAppDataVectorB = AppDataVector(kRegion, true, "vector_b");
  kAppDataVectorPressure = AppDataVector(kRegion, true, "vector_pressure");
  kAppDataVectorZ = AppDataVector(kRegion, true, "vector_z");
  kAppDataVectorTemp = AppDataVector(kRegion, true, "vector_temp");
}

} // namespace application

