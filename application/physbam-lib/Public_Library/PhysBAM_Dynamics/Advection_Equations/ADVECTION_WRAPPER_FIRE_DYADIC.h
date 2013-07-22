#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_FIRE_DYADIC
//#####################################################################
#ifndef __ADVECTION_WRAPPER_FIRE_DYADIC__
#define __ADVECTION_WRAPPER_FIRE_DYADIC__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Dynamics/Interpolation/FIRE_INTERPOLATION_FORWARD.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_NESTED_ADVECTION>
class ADVECTION_WRAPPER_FIRE_DYADIC:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::LINEAR_INTERPOLATION_DYADIC_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;typedef typename T_GRID::PROJECTION T_PROJECTION;
    typedef typename T_GRID::BOUNDARY_SCALAR T_BOUNDARY;typedef typename T_BOUNDARY::template REBIND<T2>::TYPE T_BOUNDARY_T2;
public:
    T_NESTED_ADVECTION& nested_advection;
    const T_PROJECTION& projection;
    const ARRAY<T>& phi_n_plus_one;

    ADVECTION_WRAPPER_FIRE_DYADIC(T_NESTED_ADVECTION& nested_advection_input,const T_PROJECTION& projection_input,const ARRAY<T>& phi_n_plus_one_input)
        :nested_advection(nested_advection_input),projection(projection_input),phi_n_plus_one(phi_n_plus_one_input)
    {}


    void Update_Advection_Equation_Node(const T_GRID& grid,ARRAY<T2>& Z,const ARRAY<T2>& Z_ghost,
        const ARRAY<TV>& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const ARRAY<T2>* Z_min_ghost=0,const ARRAY<T2>* Z_max_ghost=0,ARRAY<T2>* Z_min=0,ARRAY<T2>* Z_max=0)
    {nested_advection.Update_Advection_Equation_Node(grid,Z,Z_ghost,V,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2>& Z,const ARRAY<T2>& Z_ghost,
        const FACE_LOOKUP_DYADIC<T_GRID>& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const ARRAY<T2>* Z_min_ghost,const ARRAY<T2>* Z_max_ghost,ARRAY<T2>* Z_min,ARRAY<T2>* Z_max)
    {ARRAY<T> phi_n_plus_one_face(grid.number_of_faces,false);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Faces(grid,phi_n_plus_one,phi_n_plus_one_face);
    FACE_LOOKUP_FIRE_DYADIC<T_GRID> V_lookup(face_velocities.V_face,projection,&phi_n_plus_one,&phi_n_plus_one_face);
    nested_advection.Update_Advection_Equation_Cell_Lookup(grid,Z,Z_ghost,V_lookup,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,ARRAY<T>& Z,const FACE_LOOKUP_DYADIC<T_GRID>& Z_ghost,
        const FACE_LOOKUP_DYADIC<T_GRID>& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const FACE_LOOKUP_DYADIC<T_GRID>* Z_min_ghost,const FACE_LOOKUP_DYADIC<T_GRID>* Z_max_ghost,ARRAY<T>* Z_min,ARRAY<T>* Z_max)
    {ARRAY<T>& phi_n=projection.elliptic_solver->levelset->phi; //assumes poisson's internal levelset is at time n
    ARRAY<T> phi_n_face(grid.number_of_faces,false);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Faces(grid,phi_n,phi_n_face);
    FACE_LOOKUP_FIRE_DYADIC<T_GRID> Z_ghost_lookup(Z_ghost.V_face,projection,&phi_n,&phi_n_face);
    ARRAY<T> phi_n_plus_one_face(grid.number_of_faces,false);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Faces(grid,phi_n_plus_one,phi_n_plus_one_face);
    FACE_LOOKUP_FIRE_DYADIC<T_GRID> V_lookup(face_velocities.V_face,projection,&phi_n_plus_one,&phi_n_plus_one_face);
    nested_advection.Update_Advection_Equation_Face_Lookup(grid,Z,Z_ghost_lookup,V_lookup,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

//#####################################################################
};
}
#endif
#endif
