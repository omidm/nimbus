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

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_PROTOTYPES_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_PROTOTYPES_H_

#include "applications/physbam/water//app_data_include.h"
#include "applications/physbam/water//parameters.h"

namespace application {

extern AppDataFaceArray<T> kAppDataFaceVel;
extern AppDataFaceArray<T> kAppDataFaceVelGhost;
extern AppDataFaceArray<bool> kAppDataPsiN;

extern AppDataScalarArray<T> kAppDataPhi3;
extern AppDataScalarArray<T> kAppDataPhi7;
extern AppDataScalarArray<T> kAppDataPhi8;
extern AppDataScalarArray<bool> kAppDataPsiD;

// Varibales for projection.
extern AppDataScalarArray<T> kAppDataPressure;
extern AppDataScalarArray<int> kAppDataColors;
extern AppDataScalarArray<T> kAppDataDivergence;
extern AppDataRawGridArray kAppDataIndexC2M;

extern AppDataParticleLevelsetEvolution<float> kAppDataPLE;

extern AppDataSparseMatrix kAppDataSparseMatrixA;
extern AppDataSparseMatrix kAppDataSparseMatrixC;

extern AppDataArrayM2C kAppDataArrayM2C;

extern AppDataCompressedScalarArray<float> kAppDataMetaP;

extern AppDataVector kAppDataVectorB;
extern AppDataVector kAppDataVectorPressure;
extern AppDataVector kAppDataVectorZ;
extern AppDataVector kAppDataVectorTemp;

void InitializeAppDataPrototypes(nimbus::GeometricRegion kRegion);

} // namespace application

#endif // NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_PROTOTYPES_H_
