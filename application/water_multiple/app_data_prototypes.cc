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

#include "application/water_multiple/app_data_prototypes.h"
#include "application/water_multiple/parameters.h"

namespace application {

// face velocities
AppDataFaceArray<T> kAppDataFaceVel(kDefaultRegion, 0, true, "face_vel");
AppDataFaceArray<T> kAppDataFaceVelGhost(kDefaultRegion, 3, true, "face_vel_ghost");

// psi
AppDataFaceArray<bool> kAppDataPsiN(kDefaultRegion, 1, true, "psi_n");
AppDataScalarArray<bool> kAppDataPsiD(kDefaultRegion, 1, true, "psi_d");

// phi
AppDataScalarArray<T> kAppDataPhi3(kDefaultRegion, 3, true, "phi_3");
AppDataScalarArray<T> kAppDataPhi7(kDefaultRegion, 7, true, "phi_7");
AppDataScalarArray<T> kAppDataPhi8(kDefaultRegion, 8, true, "phi_8");

// particle levelset evolution
AppDataParticleLevelsetEvolution<float> kAppDataPLE(kDefaultRegion, 3, true,
                                                "particle_container");

// pressure
AppDataScalarArray<T> kAppDataPressure(kDefaultRegion, 1, true, "pressure");

// Varibales for projection.
AppDataScalarArray<int> kAppDataColors(kDefaultRegion, 1, true,
                                   "filled_region_colors");
AppDataScalarArray<T> kAppDataDivergence(kDefaultRegion, 1, true, "divergence");

// TODO(quhang): this app_data variable is questionable, because it cannot be
// deleted if meta_p is being used.
AppDataRawGridArray kAppDataIndexC2M(kDefaultRegion, true, "index_c2m");
AppDataArrayM2C kAppDataArrayM2C(kDefaultRegion, true, "index_m2c");
AppDataCompressedScalarArray<float> kAppDataMetaP(kDefaultRegion, 1, true,
                                              "vector_p_meta_format");

AppDataSparseMatrix kAppDataSparseMatrixA(kDefaultRegion, true, "matrix_c");
AppDataSparseMatrix kAppDataSparseMatrixC(kDefaultRegion, true, "matrix_a");


AppDataVector kAppDataVectorB(kDefaultRegion, true, "vector_b");
AppDataVector kAppDataVectorPressure(kDefaultRegion, true, "vector_pressure");
AppDataVector kAppDataVectorZ(kDefaultRegion, true, "vector_z");
AppDataVector kAppDataVectorTemp(kDefaultRegion, true, "vector_temp");
} // namespace application

