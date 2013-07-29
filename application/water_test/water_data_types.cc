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
    return false;
}
