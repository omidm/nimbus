//#####################################################################
// Copyright 2011 Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_THINSHELL_ADVECTION
//#####################################################################
#ifndef __EULER_THINSHELL_ADVECTION__
#define __EULER_THINSHELL_ADVECTION__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/CUT_CELL.h>

namespace PhysBAM {
template<class T,int d> void Advect_Near_Interface_Data(
        const GRID<VECTOR<T,d> >& grid,const T collision_thickness,const T dt,const ARRAY<bool,VECTOR<int,d> >& near_interface_mask,
        const ARRAY<bool,VECTOR<int,d> >& swept_cells,const ARRAY<VECTOR<T,d+2>,FACE_INDEX<d> >& flux_boundary_conditions, const ARRAY<T,VECTOR<int,d> >& cell_volumes_np1,
        HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >& visibility_graph,
        const ARRAY<bool,VECTOR<int,d> >& psi_n,  const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_n,  const ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >& U_n,
        const ARRAY<bool,VECTOR<int,d> >& psi_np1,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_np1,      ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >& U_np1);

template<class T,int d>
void Compute_Hybrid_Boundary_Fluxes(const GRID<VECTOR<T,d> >& grid,const T dt,const ARRAY<bool,VECTOR<int,d> >& near_interface_mask,const ARRAY<bool,VECTOR<int,d> >& psi,
                                    const ARRAY<VECTOR<T,d+2>,FACE_INDEX<d> >& flux_boundary_conditions,ARRAY<TRIPLE<VECTOR<int,d>,VECTOR<int,d>,VECTOR<T,d+2> > >& hybrid_flux_data);
};
#endif
