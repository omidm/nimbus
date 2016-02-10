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

#include <string>
#include <vector>

#include "application/physbam-app/translator_test/particle_test.h"
#include "application/water_multiple/physbam_include.h"
#include "data/physbam/translator_physbam.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace test {

void ParticleTest::
DeleteOutsideParticles(PhysBAMParticleContainer *particle_levelset,
                       bool positive) {
    Translator::DeleteParticles(shift, loc_region, enl_region,
        particle_levelset, scale, positive);
}

void ParticleTest::
ReadParticles(nimbus::GeometricRegion read_region,
              nimbus::DataArray &read_array,
              PhysBAMParticleContainer *particle_levelset, bool positive) {
    Translator::ReadParticles(read_region, shift, read_array, particle_levelset,
                              scale, positive, true);
}

void ParticleTest::
WriteParticles(nimbus::GeometricRegion write_region,
              nimbus::DataArray &write_array,
              PhysBAMParticleContainer *particle_levelset, bool positive) {
    Translator::WriteParticles(write_region, shift, write_array, particle_levelset,
                              scale, positive);
}

} // namespace test
