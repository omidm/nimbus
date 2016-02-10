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

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_OPTIONS_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_OPTIONS_H_

#include "application/water_multiple/app_data_include.h"
#include "application/water_multiple/app_data_prototypes.h"

namespace application {

struct AppAppObjects {
  AppDataFaceArray<T> *fv;
  AppDataFaceArray<T> *fvg;
  AppDataFaceArray<bool> *psi_n;
  AppDataScalarArray<T> *phi3;
  AppDataScalarArray<T> *phi7;
  AppDataScalarArray<T> *phi8;
  AppDataScalarArray<bool> *psi_d;
  AppDataParticleLevelsetEvolution<T> *ple;
  AppDataScalarArray<T> *pressure;
  AppDataScalarArray<int> *color;
  AppDataScalarArray<T> *divergence;
  AppDataSparseMatrix *matrix_a;
  AppDataArrayM2C * index_m2c;
  AppDataRawGridArray *index_c2m;
  AppDataVector* vector_b;
  StaticConfigValidMask* static_config_valid_mask;
  StaticConfigUInterface* static_config_u_interface;
  StaticConfigForce* static_config_force;

  AppAppObjects() {
    fv    = NULL;
    fvg   = NULL;
    psi_n = NULL;
    phi3  = NULL;
    phi7  = NULL;
    phi8  = NULL;
    psi_d = NULL;
    ple   = NULL;
    pressure = NULL;
    color = NULL;
    divergence = NULL;
    matrix_a = NULL;
    index_m2c = NULL;
    index_c2m = NULL;
    vector_b = NULL;
    static_config_valid_mask = NULL;
    static_config_u_interface = NULL;
    static_config_force = NULL;
  }
};

}  // namespace application

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_OPTIONS_H_
