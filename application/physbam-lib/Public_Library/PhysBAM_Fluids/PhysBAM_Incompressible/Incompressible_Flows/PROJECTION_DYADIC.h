#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_DYADIC
//#####################################################################
//
// Solves V-V*+GRAD(p)=0 for a divergence free V, given V*. This is done by first solving DIV(GRAD(p))=DIV(V*).
// Input a LEVELSET_OCTREE class with phi as (1,m) by (1,n).
//
// In the flame case.
// Changes both the real and ghost values of u, v, u_fuel and v_fuel, although the new ghost values are wrong!
// u and v are real where phi is negative. u_fuel and v_fuel are real where where phi is positive.
//
//#####################################################################
#ifndef __PROJECTION_DYADIC__
#define __PROJECTION_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/PROJECTION.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_PDE_Linear/POISSON_COLLIDABLE_DYADIC.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_POLICY.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_DYADIC.h>
namespace PhysBAM{

template<class T_GRID>
class PROJECTION_DYADIC:public PROJECTION<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef DYADIC_GRID_ITERATOR_FACE<T_GRID> FACE_ITERATOR;typedef DYADIC_GRID_ITERATOR_NODE<T_GRID> NODE_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::INCOMPRESSIBLE INCOMPRESSIBLE;typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;
public:
    typedef PROJECTION<T> BASE;
    using BASE::use_non_zero_divergence;using BASE::flame;using BASE::density;

    T_GRID& grid;
    ARRAY<T> p,face_velocities;
    LAPLACE_COLLIDABLE_DYADIC<T_GRID> *elliptic_solver,*laplace;
    POISSON_COLLIDABLE_DYADIC<T_GRID>* poisson;
    ARRAY<T> divergence; // use this to set up a non-zero divergence
    ARRAY<T>* flame_speed;
    ARRAY<T> phi_face_for_flame;

    PROJECTION_DYADIC(T_GRID& grid_input,ARRAY<T>* flame_speed_input=0,const bool flame_input=false);
    ~PROJECTION_DYADIC();

    void Initialize_Grid()
    {if(poisson)poisson->Initialize_Grid();else laplace->Initialize_Grid();
    Use_Non_Zero_Divergence(use_non_zero_divergence); // call this since the grid changed
    p.Resize(grid.number_of_cells);face_velocities.Resize(grid.number_of_faces);
    if(flame)phi_face_for_flame.Resize(grid.number_of_faces);}

    void Use_Non_Zero_Divergence(const bool use_non_zero_divergence_input=true)
    {use_non_zero_divergence=use_non_zero_divergence_input;
    if(use_non_zero_divergence) divergence.Resize(grid.number_of_cells);else divergence.Resize(0);}

    T Face_Velocity_With_Ghost_Value(const ARRAY<T>& face_velocities_ghost,const int face_index,const T current_phi) const
    {assert(flame);assert(flame_speed);const T face_phi=phi_face_for_flame(face_index);
    if(!LEVELSET_UTILITIES<T>::Interface(current_phi,face_phi)) return face_velocities_ghost(face_index);
    if(face_phi>0)return face_velocities_ghost(face_index)-Face_Jump(face_index);else return face_velocities_ghost(face_index)+Face_Jump(face_index);}

    T Face_Jump(const int face_index) const
    {PHYSBAM_NOT_IMPLEMENTED();/*
    DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,face_index);
    CELL* deepest_cell=iterator.Deepest_Cell(),*other_cell=iterator.Other_Cell();
    if(other_cell&&other_cell->Depth_Of_This_Cell()==grid.maximum_depth){
        ARRAY<TV>& normals=*elliptic_solver->levelset->normals;
        T two_times_face_flame_speed=(*flame_speed)(deepest_cell->Cell())+(*flame_speed)(other_cell->Cell());
        T two_times_face_normal_component=normals(deepest_cell->Cell())[iterator.Axis()+1]+normals(other_cell->Cell())[iterator.Axis()+1];
        assert(Jump_Constant()==jump_constant);return (T).25*jump_constant*two_times_face_flame_speed*two_times_face_normal_component;}
        else return 0;*/} // TODO: revisit what to do in this case (ideally it shouldn't be called here, but currently advection semi lagrangian calls the fire lookup unconditionally)

//#####################################################################
    void Make_Divergence_Free(const T dt,const T time);
    void Enforce_Velocity_Compatibility();
    void Update_Phi_And_Move_Velocity_Discontinuity(const T_LEVELSET& levelset,const T time,const bool update_phi_only=false);
//#####################################################################
};
}
#endif
#endif
