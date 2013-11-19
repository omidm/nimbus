/* Copyright 2013 Stanford University.
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
 * Author: Zhihao Jia <zhihao@stanford.edu>
 */

#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include "Projection_Example.h"
using namespace PhysBAM;

template<class TV> PROJECTION_EXAMPLE<TV>::PROJECTION_EXAMPLE():
	initial_time(0), 
	first_frame(0),
	last_frame(100), 
	frame_rate(24), 
	restart(0),
	mac_grid(TV_INT(), RANGE<TV>::Unit_Box(), true),	 
	projection(mac_grid, false, false),
	boundary(0) {
		for (int i=1; i<=TV::dimension; i++) {
			domain_boundary(i)(1)=true;
			domain_boundary(i)(2)=true;
		}
}

template<class TV> PROJECTION_EXAMPLE<TV>::~PROJECTION_EXAMPLE() {	
		delete boundary;
}

template<class TV> void PROJECTION_EXAMPLE<TV>::Set_Boundary_Conditions(const T time) {
	projection.elliptic_solver->psi_D.Fill(false);
	projection.elliptic_solver->psi_N.Fill(false);
	for (int axis=1; axis<=TV::dimension; axis++)
		for (int axis_side=1; axis_side<=2; axis_side++) {
			int side=2*(axis-1)+axis_side;
			if (domain_boundary(axis)(axis_side)) {
				TV_INT interior_cell_offset=axis_side==1 ? TV_INT()
						: -TV_INT::Axis_Vector(axis);
				for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid, 1,
						GRID<TV>::BOUNDARY_REGION, side); iterator.Valid(); iterator.Next()) {
					TV_INT cell=iterator.Face_Index()+interior_cell_offset;
					TV_INT boundary_face=axis_side==1 ? iterator.Face_Index()
							+TV_INT::Axis_Vector(axis) : iterator.Face_Index()
							-TV_INT::Axis_Vector(axis);
					projection.elliptic_solver->psi_D(cell)=true;
					projection.p(cell)=0;
				}
			}
		}
}
