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
 */

#ifndef NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_PHYSBAM_UTILS_H_
#define NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_PHYSBAM_UTILS_H_

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/physbam_include.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "shared/nimbus.h"

namespace application {

    // Initialzes WATER_EXAMPLE and WATER_DRIVER with the given input parameters
    // and the simulation variables in data array. Notice, this function call
    // is not expected to be called at the start of the simulation when
    // WATER_EXAMPLE and WATER_DRIVER is not passed by Nimbus but initialized by
    // PhysBAM itself.  --quhang
    bool InitializeExampleAndDriver(
        const nimbus::DataArray& da,
        const int current_frame,
        const T time,
        const nimbus::Job* job,
        PhysBAM::WATER_EXAMPLE<TV>*& example,
        PhysBAM::WATER_DRIVER<TV>*& driver);

    // Destroys WATER_EXAMPLE and WATER_DRIVER.
    void DestroyExampleAndDriver(
        PhysBAM::WATER_EXAMPLE<TV>*& example,
        PhysBAM::WATER_DRIVER<TV>*& driver);

} // namespace application

#endif  // NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_PHYSBAM_UTILS_H_
