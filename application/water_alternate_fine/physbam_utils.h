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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 *
 * This file contains PhysBAM related utilities, e.g., the utilities for
 * initializing WATER_DRIVER and WATER_EXAMPLE.
 * Utilities that are heavily related to PhysBAM libraries should be put here,
 * while utilities that only depend on Nimbus should go to app_util.h.
 *
 * Modifier: Hang Qu <quhang@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_PHYSBAM_UTILS_H_
#define NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_PHYSBAM_UTILS_H_

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/physbam_include.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "shared/nimbus.h"

namespace application {

// Data structure to specify how to initialize WATER_EXAMPLE and WATER_DRIVER.
struct InitConfig {
  int frame;
  T time;
  bool init_phase;
  bool set_boundary_condition;
  TV_INT grid_size;

  InitConfig() {
    frame = 0;
    time = 0;
    init_phase = false;
    set_boundary_condition = true;
    grid_size = TV_INT::All_Ones_Vector() * kScale;
  }
};

// Initialzes WATER_EXAMPLE and WATER_DRIVER with the given "init_config"
// and the simulation variables in data array.
// Returns false if it fails.
bool InitializeExampleAndDriver(
    const InitConfig& init_config,
    const nimbus::Job* job,
    const nimbus::DataArray& da,
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver);

// Destroys WATER_EXAMPLE and WATER_DRIVER.
void DestroyExampleAndDriver(
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver);

}  // namespace application

#endif  // NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_PHYSBAM_UTILS_H_
