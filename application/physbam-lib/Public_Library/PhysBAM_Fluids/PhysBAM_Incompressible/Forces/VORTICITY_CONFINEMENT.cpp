//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> VORTICITY_CONFINEMENT<T_GRID>::
VORTICITY_CONFINEMENT(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_body_list,T_FACE_ARRAYS_BOOL* valid_mask,const bool use_variable_vorticity_confinement,const T vorticity_confinement)
    :collision_body_list(collision_body_list),valid_mask(valid_mask),use_variable_vorticity_confinement(use_variable_vorticity_confinement),vorticity_confinement(vorticity_confinement)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> VORTICITY_CONFINEMENT<T_GRID>::
~VORTICITY_CONFINEMENT()
{
}
//#####################################################################
// Function Apply_Vorticity_Confinement_Force
//#####################################################################
template<class T_GRID,class T_ARRAYS_TV> static void
Apply_Vorticity_Confinement_Force_Helper(const T_GRID& grid,typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& face_velocities,T_ARRAYS_TV& F,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_body_list)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    // want cells to face averaging here
    if(collision_body_list){
        AVERAGING_COLLIDABLE_UNIFORM<T_GRID,FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> > vorticity_averaging_collidable(*collision_body_list,T());
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities.Component(axis)(iterator.Face_Index())+=vorticity_averaging_collidable.Cell_To_Face(grid,axis,iterator.Face_Index(),F);}}
    else
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
            face_velocities.Component(axis)(iterator.Face_Index())+=F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis];}
}
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<float,1> >&,ARRAY<float,FACE_INDEX<1> >&,ARRAY<VECTOR<float,1> ,VECTOR<int,1> >&,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<float,1> > >*)
{PHYSBAM_NOT_IMPLEMENTED();}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<double,1> >&,ARRAY<double,FACE_INDEX<1> >&,ARRAY<VECTOR<double,1> ,VECTOR<int,1> >&,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<double,1> > >*)
{PHYSBAM_NOT_IMPLEMENTED();}
#endif
template<class T_GRID> void VORTICITY_CONFINEMENT<T_GRID>::
Apply_Vorticity_Confinement_Force(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities,T_ARRAYS_VECTOR& F)
{
    Apply_Vorticity_Confinement_Force_Helper(grid,face_velocities,F,collision_body_list);
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force
//#####################################################################
template<class T_GRID,class T_FACE_ARRAYS,class T_ARRAYS_TV> static void 
Compute_Vorticity_Confinement_Force_Helper(const T_GRID& grid,const T_FACE_ARRAYS& face_velocities_ghost,T_ARRAYS_TV& F,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_body_list,const typename T_FACE_ARRAYS::template REBIND<bool>::TYPE* valid_mask)
{
    typedef typename T_GRID::VECTOR_T TV;
    if(TV::dimension==1){PHYSBAM_NOT_IMPLEMENTED();}
    typedef typename TV::SCALAR T;typedef typename T_ARRAYS_TV::template REBIND<T>::TYPE T_ARRAYS_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_ARRAYS_TV::template REBIND<typename TV::SPIN>::TYPE T_ARRAYS_SPIN;typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    T_ARRAYS_SPIN vorticity(grid.Cell_Indices(2),false);
    T_ARRAYS_SCALAR vorticity_magnitude(grid.Cell_Indices(2));
    if(collision_body_list){
        FACE_LOOKUP_UNIFORM<T_GRID> face_velocities_lookup_uniform(face_velocities_ghost);
        FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> face_velocities_lookup(face_velocities_lookup_uniform,*collision_body_list,valid_mask);
        VORTICITY_UNIFORM<TV>::Vorticity(grid,face_velocities_lookup,vorticity,vorticity_magnitude);}
    else VORTICITY_UNIFORM<TV>::Vorticity(grid,T_FACE_LOOKUP(face_velocities_ghost),vorticity,vorticity_magnitude);
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){ // do collision awareness when these are averaged to faces
        TV vortex_normal_vector=T_LEVELSET::Normal_At_Node(grid,vorticity_magnitude,iterator.Cell_Index());
        F(iterator.Cell_Index())=TV::Cross_Product(vortex_normal_vector,vorticity(iterator.Cell_Index()));}
}
template<class T_GRID> void VORTICITY_CONFINEMENT<T_GRID>::
Compute_Vorticity_Confinement_Force(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T_FACE_ARRAYS_BOOL* valid_mask,T_ARRAYS_VECTOR& F)
{
    Compute_Vorticity_Confinement_Force_Helper(grid,face_velocities_ghost,F,collision_body_list,valid_mask);
}
//#####################################################################
// Function Add_Explicit_Forces
//#####################################################################
template<class T_GRID> void VORTICITY_CONFINEMENT<T_GRID>::
Add_Explicit_Forces(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    T_ARRAYS_VECTOR F(grid.Cell_Indices(1),false);
    Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,valid_mask,F);
    if(collision_body_list){
        if(use_variable_vorticity_confinement){F*=dt;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement;}
    else{
        if(use_variable_vorticity_confinement){F*=dt*(T).5;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement*(T).5;}
    Apply_Vorticity_Confinement_Force(grid,face_velocities,F);
}
//#####################################################################
// Function Add_Implicit_Forces_Projection
//#####################################################################
template<class T_GRID> void VORTICITY_CONFINEMENT<T_GRID>::
Add_Implicit_Forces_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_GRID> void VORTICITY_CONFINEMENT<T_GRID>::
Initialize_Grids(const T_GRID& grid)
{
    if(use_variable_vorticity_confinement) variable_vorticity_confinement.Resize(grid.Cell_Indices(1));
    else variable_vorticity_confinement.Clean_Memory();
}
//#####################################################################
template class VORTICITY_CONFINEMENT<GRID<VECTOR<float,1> > >;
template class VORTICITY_CONFINEMENT<GRID<VECTOR<float,2> > >;
template class VORTICITY_CONFINEMENT<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTICITY_CONFINEMENT<GRID<VECTOR<double,1> > >;
template class VORTICITY_CONFINEMENT<GRID<VECTOR<double,2> > >;
template class VORTICITY_CONFINEMENT<GRID<VECTOR<double,3> > >;
#endif
