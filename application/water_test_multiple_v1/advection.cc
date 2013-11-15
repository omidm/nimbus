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
 * Methods used in advection in water application.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "advection.h"
#include "app_config.h"
#include "data_face_arrays.h"
#include "physbam_include.h"
#include "shared/geometric_region.h"
#include "water_app.h"
#include "water_data_driver.h"

void Advect_Velocities (
        ::nimbus::GeometricRegion region,
        T_FACE_ARRAY *face_velocities,
        T_FACE_ARRAY *face_vel_extended,
        WaterApp *water_app,
        T dt,
        T time) {

    typedef ::PhysBAM::GRID<TV> T_GRID;
    typedef ::PhysBAM::RANGE<TV> T_RANGE;

    TV_INT size(region.dx(), region.dy());
    T_GRID grid(size, T_RANGE::Unit_Box(), true);
    water_app->advection_scalar()->Update_Advection_Equation_Face(
            grid,
            *face_velocities,
            *face_vel_extended,
            *face_vel_extended,
            *water_app->boundary(),
            dt,
            time);
}
