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
 * Helper function for projection job. Still in progress.
 * Author: Hang Qu<quhang@stanford.edu>
 */

#include "application/water_multiple/app_utils.h"

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_PROJECTION_HELPER_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_PROJECTION_HELPER_H_

namespace PhysBAM {

typedef application::TV TV;
typedef application::TV_INT TV_INT;

typedef ARRAY<bool, VECTOR<int, TV::dimension> > T_PSI_D;
typedef ARRAY<bool, FACE_INDEX<TV::dimension> > T_PSI_N;
typedef ARRAY<int, VECTOR<int, TV::dimension> > T_COLOR;
typedef GRID<TV> T_GRID;

/* Fills in the color for projection. Assume only one colors. */
void FillUniformRegionColor(
    const T_GRID& grid,
    const T_PSI_D& psi_D,
    const T_PSI_N& psi_N,
    const bool solve_single_cell_neumann_regions,
    T_COLOR* filled_region_colors);

bool All_Cell_Faces_Neumann(
    const TV_INT& cell_index, const T_PSI_N& psi_N);

}  // namespace PhysBAM
#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_PROJECTION_HELPER_H_
