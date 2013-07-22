//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OBJECT_SOLID_VELOCITY  
//#####################################################################
#ifndef __BOUNDARY_OBJECT_SOLID_VELOCITY__
#define __BOUNDARY_OBJECT_SOLID_VELOCITY__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>

namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID>
class BOUNDARY_OBJECT_SOLID_VELOCITY:public BOUNDARY_OBJECT<T_GRID,VECTOR<typename T_GRID::SCALAR,T_GRID::dimension+2> >
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
    enum{d=T_GRID::dimension+2};
public:

    BOUNDARY_OBJECT_SOLID_VELOCITY()
    {}

    ~BOUNDARY_OBJECT_SOLID_VELOCITY()
    {}
//#####################################################################
    void Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis) PHYSBAM_OVERRIDE;
    void Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& object_velocity,const T unused) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
