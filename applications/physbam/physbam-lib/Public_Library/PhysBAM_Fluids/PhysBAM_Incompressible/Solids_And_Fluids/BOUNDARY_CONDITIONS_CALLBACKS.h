//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONDITIONS_CALLBACKS
//##################################################################### 
#ifndef __BOUNDARY_CONDITIONS_CALLBACKS__
#define __BOUNDARY_CONDITIONS_CALLBACKS__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV>
class BOUNDARY_CONDITIONS_CALLBACKS
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    enum WORKAROUND{d=TV::m};
public:
    enum RAY_TYPE {unused,inside,outside,slip,noslip,free};

    BOUNDARY_CONDITIONS_CALLBACKS();
    virtual ~BOUNDARY_CONDITIONS_CALLBACKS();

//#####################################################################
    virtual void Mark_Outside(ARRAY<bool,FACE_INDEX<d> >& outside);
    virtual void Mark_Outside(ARRAY<bool,TV_INT>& outside);
    virtual RAY_TYPE Get_Boundary_Along_Ray(const FACE_INDEX<d>& f1,const FACE_INDEX<d>& f2,T& theta,T& value);
    virtual RAY_TYPE Get_Boundary_Along_Ray(const TV_INT& c1,const TV_INT& c2,T& theta,T& value);
//#####################################################################
};
}
#endif
