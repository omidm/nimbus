//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OBJECT  
//#####################################################################
#ifndef __BOUNDARY_OBJECT__
#define __BOUNDARY_OBJECT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID,class TV_DIMENSION>
class BOUNDARY_OBJECT
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef GRID<VECTOR<T,TV::dimension-1> > T_GRID_LOWER_DIM;typedef typename T_GRID_LOWER_DIM::VECTOR_INT TV_INT_LOWER_DIM;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
public:
    BOUNDARY_OBJECT();
    virtual ~BOUNDARY_OBJECT();

//#####################################################################
    void Get_State_At_Location(const GRID<VECTOR<T,1> >& grid_1d,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_1d,const T location,const VECTOR<int,2>& region_boundaries,TV_DIMENSION& u_1d);
    void Fill_Ghost_Cells_Neumann(const GRID<VECTOR<T,1> >& grid_1d,ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_1d,const T_FACE_ARRAYS_SCALAR& face_velocities,const TV_INT_LOWER_DIM& node_lower_dimension,const int axis,
        const int ghost_cells,const bool use_exact_neumann_face_location,const VECTOR<int,2>& domain,const VECTOR<int,2>& region_boundaries,const VECTOR<bool,2>& psi_N,
        CONSERVATION_CALLBACKS<T_GRID,TV_DIMENSION>* callbacks);
    virtual void Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis)=0;
    virtual void Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& normal,const T object_velocity_normal_component)=0;
//#####################################################################
};
}
#endif
