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
#include "data_face_arrays.h"
#include "physbam_include.h"
#include "types.h"
#include "water_data_driver.h"
#include "water_driver.h"

void Advection (
        ::water_app_data::FaceArray<TV> *face_velocities,
        NonAdvData<TV, T> *sim_data) {

    typedef typename ::PhysBAM::GRID<TV> T_GRID;
    typedef typename ::PhysBAM::GRID_ARRAYS_POLICY<GRID<TV> >
        ::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename 
        ::PhysBAM::ARRAY<T, ::PhysBAM::FACE_INDEX<TV::dimension> >
        T_FACE_ARRAY;

    T_FACE_ARRAYS_SCALAR face_velocities_ghost;
    face_velocities_ghost.Resize(
            sim_data->incompressible->grid,
            sim_data->number_of_ghost_cells,
            false);
    sim_data->incompressible->boundary->Fill_Ghost_Cells_Face(
            sim_data->incompressible->grid,
            *face_velocities->data,
            face_velocities_ghost,
            sim_data->time + sim_data->dt,
            sim_data->number_of_ghost_cells);

    //TODO: serialize/ deserialize, advection needs:
    //sim_data->incompressible->advection (probably needed only for advect V)
    //sim_data->incompressible->boundary (needed elsewhere)
    //sim_data->dt (parameter)
    //face_velocities (needed elsewhere)
    //face_velocities_ghost (needed elsewhere)
    sim_data->incompressible->advection->Update_Advection_Equation_Face(
            *face_velocities->grid,
            *face_velocities->data,
            face_velocities_ghost,
            face_velocities_ghost,
            *sim_data->incompressible->boundary,
            sim_data->dt,
            sim_data->time);
}
