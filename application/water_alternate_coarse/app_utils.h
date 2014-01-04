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
 * Definitions and typedef useful for application, data and jobs.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_ALTERNATE_COARSE_APP_UTILS_H_
#define NIMBUS_APPLICATION_WATER_ALTERNATE_COARSE_APP_UTILS_H_

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"
#include "shared/nimbus_types.h"
#include "worker/physical_data_instance.h"

#define APP_LOG DBG_TEMP
#define APP_LOG_STR "temp"
#define TRANSLATE_STR "translate"

#ifndef APP_FACE_VEL
#define APP_FACE_VEL "face_vel"
#endif

#ifndef APP_PHI
#define APP_PHI "phi"
#endif

#ifndef APP_PRESSURE
#define APP_PRESSURE "pressure"
#endif

#ifndef APP_POS_PARTICLES
#define APP_POS_PARTICLES "pos_particles"
#endif

#ifndef APP_NEG_PARTICLES
#define APP_NEG_PARTICLES "neg_particles"
#endif

#ifndef APP_POS_REM_PARTICLES
#define APP_POS_REM_PARTICLES "pos_rem_particles"
#endif

#ifndef APP_NEG_REM_PARTICLES
#define APP_NEG_REM_PARTICLES "neg_rem_particles"
#endif

namespace application {

    // simulation dimension
    const int kDimension = 3;

    // typedefs
    typedef float T;
    typedef float RW;
    typedef PhysBAM::VECTOR<T,   kDimension> TV;
    typedef PhysBAM::VECTOR<int, kDimension> TV_INT;

    // application specific parameters and constants
    const int kThreadsNum = 1;
    const int kScale = 30;
    const int kGhostNum = 3;
    const int kPressureGhostNum = 1;
    const int kLastFrame = 15;
    const std::string kOutputDir = "output";
    // follow physbam convenctions here, otherwise translator becomes messy
    const GeometricRegion kDomain(1, 1, 1, kScale, kScale, kScale);
    const GeometricRegion kDomainGhost(-kGhostNum + 1,
                                       -kGhostNum + 1,
                                       -kGhostNum + 1,
                                       kScale + kGhostNum*2,
                                       kScale + kGhostNum*2,
                                       kScale + kGhostNum*2);
    const GeometricRegion kDomainFaceVel = kDomain;
    const GeometricRegion kDomainPhi = kDomainGhost;
    const GeometricRegion kDomainParticles = kDomainGhost;
    const GeometricRegion kDomainPressure(-kPressureGhostNum + 1,
                                          -kPressureGhostNum + 1,
                                          -kPressureGhostNum + 1,
                                          kScale + kPressureGhostNum*2,
                                          kScale + kPressureGhostNum*2,
                                          kScale + kPressureGhostNum*2);

    const int_dimension_t kFaceVelBufSize = kScale *
                                            kScale * 
                                            (kScale+1) *
                                            kDimension * sizeof(T);
    const int_dimension_t kPhiBufSize = (kScale + 2*kGhostNum) *
                                        (kScale + 2*kGhostNum) *
                                        (kScale + 2*kGhostNum) * sizeof(T);
    const int_dimension_t kPressureBufSize = (kScale + 2*kPressureGhostNum) *
                                             (kScale + 2*kPressureGhostNum) *
                                             (kScale + 2*kPressureGhostNum) * sizeof(T);
    const int_dimension_t kParticlesBufSize = 0;

    // TODO: some hacks that need to be cleaned soon after a meeting/
    // discussion -- one option is to make region a part of data, and
    // let nimbus take care of initializing region correctly when creating
    // the data object
    bool GetTranslatorData(const nimbus::Job *job,
                           const std::string &name,
                           const nimbus::DataArray& da,
                           nimbus::PdiVector *vec);
    void DestroyTranslatorObjects(nimbus::PdiVector *vec);

   // TODO: lets make read/ write sets if possible, and also have separate
   // read/ write instead of one DataArray passed to a job/ a better indexing
    bool Contains(nimbus::IDSet<nimbus::logical_data_id_t> data_set,
                  nimbus::logical_data_id_t  id);

} // namespace application

#endif  // NIMBUS_APPLICATION_WATER_ALTERNATE_COARSE_APP_UTILS_H_
