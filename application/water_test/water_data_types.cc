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

#include <iostream>
#include "assert.h"
#include "shared/nimbus.h"
#include "./water_data_types.h"

using namespace PhysBAM;
using nimbus::Data;

// TODO(someone): right now, create may/ may not allocate required amount of
// memory. Will have to dig deep into PhysBAM to allocate required amount of
// memory at the beginning itself.

template <class TV> FaceArray<TV>::
FaceArray(int size)
{
    this->size_ = size;
    grid = NULL;
    data = NULL;
}

template <class TV> void FaceArray<TV>::
create()
{
    std::cout << "Creating FaceArray\n";

    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    grid = new GRID<TV> (TV_INT::All_Ones_Vector()*size_,
            RANGE<TV>::Unit_Box(), true);
    assert(grid);

    data = new  ARRAY<T,FACE_INDEX<TV::dimension> >();
    assert(data);
}

template <class TV> Data* FaceArray<TV>::
clone()
{
    std::cout << "Cloning facearray\n";
    return new FaceArray<TV>(size_);
}

template <class TV> bool FaceArray<TV>::
initialize()
{
    return false;
}

template <class TV> FaceArrayGhost<TV>::
FaceArrayGhost(int size)
{
    this->size_ = size;
    grid = NULL;
    data = NULL;
}

template <class TV> void FaceArrayGhost<TV>::
create()
{
    std::cout << "Creating FaceArrayGhost\n";

    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    grid = new GRID<TV> (TV_INT::All_Ones_Vector()*size_,
            RANGE<TV>::Unit_Box(), true);
    assert(grid);

    data = new typename GRID_ARRAYS_POLICY<GRID<TV> >::
        FACE_ARRAYS();
    assert(data);
}

template <class TV> Data* FaceArrayGhost<TV>::
clone()
{
    std::cout << "Cloning facearrayghost\n";
    return new FaceArrayGhost<TV>(size_);
}

template <class TV> bool FaceArrayGhost<TV>::
initialize()
{
    return false;
}

template <class TV, class T> NonAdvData<TV, T>::
NonAdvData(int size)
{
    this->size_ = size;

    number_of_ghost_cells = 3;

    grid = NULL;

    boundary_scalar = NULL;
    boundary = NULL;
    phi_boundary = NULL;
    phi_boundary_water = NULL;
    domain_boundary = NULL;

    sources = NULL;

    particle_levelset_evolution = NULL;
    advection_scalar = NULL;

    collision_bodies_affecting_fluid = NULL;

    projection = NULL;
    incompressible = NULL;
}

template <class TV, class T> void NonAdvData<TV, T>::
create()
{
    std::cout << "Creating NonAdvData\n";

    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    grid = new GRID<TV> (TV_INT::All_Ones_Vector()*size_,
            RANGE<TV>::Unit_Box(), true);
    assert(grid);

    boundary_scalar = new BOUNDARY_UNIFORM<GRID<TV>, T>();
    phi_boundary_water = new
        typename GEOMETRY_BOUNDARY_POLICY<GRID<TV> >::
        BOUNDARY_PHI_WATER();
    domain_boundary = new VECTOR<VECTOR<bool,2>,TV::dimension>();

    particle_levelset_evolution = new
        PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >
        (*grid, number_of_ghost_cells);
    advection_scalar = new ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T>();

    collision_bodies_affecting_fluid = new
        typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::
        GRID_BASED_COLLISION_GEOMETRY(*grid);

    projection = new PROJECTION_DYNAMICS_UNIFORM< GRID<TV> >
        (*grid, false, false, false, false, NULL);
    incompressible = new INCOMPRESSIBLE_UNIFORM<GRID<TV> >(*grid, *projection);
}

template <class TV, class T> Data* NonAdvData<TV, T>::
clone()
{
    std::cout << "Cloning nonadvdata\n";
    return new NonAdvData<TV, T>(size_);
}

template <class TV, class T> bool NonAdvData<TV, T>::
initialize()
{
    return false;
}

#ifndef TEMPLATE_USE
#define TEMPLATE_USE
typedef VECTOR<float, 2> TVF2;
typedef float TF;
#endif  // TEMPLATE_USE

template class FaceArray<TVF2>;
template class FaceArrayGhost<TVF2>;
template class NonAdvData<TVF2, TF>;
