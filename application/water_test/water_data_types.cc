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
 * Helper functions in water_driver.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "lib/nimbus.h"
#include "./water_data_types.h"

template <class TV> FaceArray<TV>::
FaceArray()
{
    data = NULL;
}

template <class TV> bool FaceArray<TV>::
initialize()
{
    data = new ARRAY<T, FACE_INDEX<TV::dimension> >;
    if (data != NULL)
        return true;
    else
        return false;
}

template <class TV> FaceArrayGhost<TV>::
FaceArrayGhost()
{
    data = NULL;
}

template <class TV> bool FaceArrayGhost<TV>::
initialize()
{
    data = new typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS;
    if (data != NULL)
        return true;
    else
        return false;
}

template <class TV> Grid<TV>::
Grid()
{
    data = NULL;
}

template <class TV> bool Grid<TV>::
initialize(
        const TV_INT &counts,
        const RANGE<TV> &box,
        const bool MAC_grid
        )
{
    data = new GRID<TV>(counts, box, MAC_grid);
    if (data != NULL)
        return true;
    else
        return false;
}

template <class TV> MPIGrid<TV>::
MPIGrid()
{
    data = NULL;
}

template <class TV> bool MPIGrid<TV>::
initialize()
{
    data = new MPI_UNIFORM_GRID<GRID<TV> >;
    if (data != NULL)
        return true;
    else
        return false;
}

template <class TV, class T> NonAdvData<TV, T>::
NonAdvData()
{
}

template <class TV, class T> bool NonAdvData<TV, T>::
initialize()
{
    // TODO: to fill in
    // nothing for projection
    boundary_scalar = new  BOUNDARY_UNIFORM<GRID<TV>, T>;
    phi_boundary_water = new typename GEOMETRY_BOUNDARY_POLICY<GRID<TV> >::
        BOUNDARY_PHI_WATER;
    domain_boundary = new VECTOR<VECTOR<bool, 2>, TV::dimension>;
    sources = new ARRAY<IMPLICIT_OBJECT<TV>*>;
    particle_levelset_evolution = new 
        PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >;
    advection_scalar = new ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>, T>;
    rigid_geometry_collection = new RIGID_GEOMETRY_COLLECTION<TV>;
    collision_bodies_affecting_fluid = new typename
        COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::
        GRID_BASED_COLLISION_GEOMETRY;
    incompressible = new INCOMPRESSIBLE_UNIFORM<GRID<TV> >;
    kinematic_evolution = new KINEMATIC_EVOLUTION<TV>;
    if (boundary_scalar == NULL || phi_boundary_water == NULL ||
            domain_boundary == NULL || sources == NULL ||
            particle_levelset_evolution == NULL ||
            advection_scalar == NULL ||
            rigid_geometry_collection == NULL ||
            collision_bodies_affecting_fluid == NULL ||
            incompressible == NULL)
        return false;
    else
        return true;
}
