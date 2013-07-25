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
 * Data types used by the application jobs and functions.
 */

#ifndef NIMBUS_APPLICATION_WATER_TEST_WATER_DATA_TYPES_H_
#define NIMBUS_APPLICATION_WATER_TEST_WATER_DATA_TYPES_H_

/* Include relevant PhysBAM files here.
 */
#include "./physbam_include.h"

using namespace PhysBAM;    // NOLINT

/* Face array for storing quantities like face velocities.
 */
template <class TV>
class FaceArray : public Data {
    typedef typename TV::SCALAR T;
    ARRAY<T, FACE_INDEX<TV::dimension> > *facedata;
};

/* Ghost face array for storing scalar quantities.
 */
template <class TV>
class FaceArrayGhost : public Data {
    typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS *faceghost;
};

/* Add all other data used by the water simulation here.  DO NOT add scalar
 * values. Scalar values can be passed around directly as parameters.
 */
template <class TV>
class NonVelAdvData : public Data {
};

#endif  // NIMBUS_APPLICATION_WATER_TEST_WATER_DATA_TYPES_H_
