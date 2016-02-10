//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_INDEX_MAP_UNIFORM
//#####################################################################
#ifndef __LEVELSET_INDEX_MAP_UNIFORM__
#define __LEVELSET_INDEX_MAP_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;
template<class TV>
class LEVELSET_INDEX_MAP_UNIFORM:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
public:
    const GRID<TV>& grid;

    BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback;
    ARRAY<bool,FACE_INDEX<d> > inside;
    ARRAY<int,FACE_INDEX<d> > face_to_index;
    ARRAY<FACE_INDEX<d> > index_to_face;
    VECTOR<bool,d> periodic_boundary;

    LEVELSET_INDEX_MAP_UNIFORM(const GRID<TV>& grid_input,BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input);
    ~LEVELSET_INDEX_MAP_UNIFORM();

    void Compute(int axis,VECTOR<bool,d> periodic_boundary_input);
    void Gather(const ARRAY<T,FACE_INDEX<d> >& u,VECTOR_ND<T>& v) const;
    void Scatter(const VECTOR_ND<T>& u,ARRAY<T,FACE_INDEX<d> >& v) const;

};
}
#endif
