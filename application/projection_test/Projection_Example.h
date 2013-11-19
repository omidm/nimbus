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

#ifndef NIMBUS_APPLICATION_PROJECTION_TEST_PROJECTION_EXAMPLE_H_
#define NIMBUS_APPLICATION_PROJECTION_TEST_PROJECTION_EXAMPLE_H_
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM {

template<class TV> class PROJECTION_EXAMPLE {
	typedef typename TV::SCALAR T;
	typedef typename TV::template REBIND<int>::TYPE TV_INT;
	enum workaround1 {d=TV::m};

public:
	STREAM_TYPE stream_type;
	T initial_time;
	int first_frame, last_frame;
	T frame_rate;
	int restart;
	std::string frame_title;

	GRID<TV> mac_grid;
	PROJECTION_UNIFORM<GRID<TV> > projection;
	ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
	BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
	BOUNDARY_UNIFORM<GRID<TV>,T> *boundary;
	VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;

	PROJECTION_EXAMPLE();
	virtual ~PROJECTION_EXAMPLE();

	T Time_At_Frame(const int frame) const {
		return initial_time+(frame-first_frame)/frame_rate;
	}

	void Initialize_Grid(TV_INT counts, RANGE<TV> domain) {
		mac_grid.Initialize(counts, domain, true);
	}

	void Initialize_Fields() {
		for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid); iterator.Valid(); iterator.Next())
			face_velocities(iterator.Full_Index())
					=(T)iterator.Face_Index()(iterator.Axis())
							/(T)mac_grid.counts(iterator.Axis());
	}

	virtual void Set_Boundary_Conditions(const T time);

	//#####################################################################
};
}

#endif /*NIMBUS_APPLICATION_PROJECTION_TEST_PROJECTION_EXAMPLE_H_*/
