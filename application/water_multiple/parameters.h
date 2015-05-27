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

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_PARAMETERS_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_PARAMETERS_H_

#include "application/water_multiple/physbam_include.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "application_utils/partition_handler.h"

#define APP_LOG DBG_TEMP
#define APP_LOG_STR "temp"
#define TRANSLATE_STR "translate"


// chnage the default variables here -omidm
#define DEFAULT_SCALE 40

#define DEFAULT_APP_PART_NUM_X 2
#define DEFAULT_APP_PART_NUM_Y 2
#define DEFAULT_APP_PART_NUM_Z 1

#define DEFAULT_APP_PROJ_PART_NUM_X 1
#define DEFAULT_APP_PROJ_PART_NUM_Y 2
#define DEFAULT_APP_PROJ_PART_NUM_Z 1

#define DEFAULT_LAST_FRAME 10

#define DEFAULT_MAX_ITERATIONS 40

#define DEFAULT_USE_GLOBAL_WRITE true



namespace application {

    // Partition Handler used by the jobs. 
    extern nimbus::PartitionHandler ph;

    // simulation dimension
    const int kDimension = 3;

    // typedefs
    typedef float T;
    typedef float RW;
    typedef PhysBAM::VECTOR<T,   kDimension> TV;
    typedef PhysBAM::VECTOR<int, kDimension> TV_INT;
    typedef typename PhysBAM::FACE_INDEX<TV::dimension> FaceIndex;
    typedef typename PhysBAM::ARRAY<T, FaceIndex> FaceArray;

    // application specific variables, set by WaterApp::Load() -omidm
    extern uint64_t kScale;
    extern uint64_t kAppPartNum;
    extern uint64_t kAppPartNumX;
    extern uint64_t kAppPartNumY;
    extern uint64_t kAppPartNumZ;
    extern uint64_t kProjAppPartNum;
    extern uint64_t kProjAppPartNumX;
    extern uint64_t kProjAppPartNumY;
    extern uint64_t kProjAppPartNumZ;
    extern nimbus::GeometricRegion kDefaultRegion;
    extern uint64_t kLastFrame;
    extern uint64_t kMaxIterations;
    extern bool kUseGlobalWrite;

    // constant application specific parameters.
    const int kThreadsNum = 1;
    const int kGhostNum = 3;
    const int kGhostW[3] = {kGhostNum, kGhostNum, kGhostNum};
    const int kPressureGhostNum = 1;
    const std::string kOutputDir = "output";
} // namespace application

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_PARAMETERS_H_
