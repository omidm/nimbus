//#####################################################################
// Copyright 2009, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/AVERAGING_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
using namespace PhysBAM;
template<class T_GRID> LEVELSET_ADVECTION<T_GRID>::
LEVELSET_ADVECTION(T_LEVELSET* _levelset)
    :levelset(_levelset),advection(0),nested_semi_lagrangian_collidable(0),semi_lagrangian_collidable(0)
{
    Set_Reinitialization_Runge_Kutta_Order();
    Set_Reinitialization_CFL();
    Use_WENO_For_Reinitialization();
}
template<class T_GRID> LEVELSET_ADVECTION<T_GRID>::
~LEVELSET_ADVECTION()
{
    delete nested_semi_lagrangian_collidable;
    delete semi_lagrangian_collidable;
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Advection
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION<T_GRID>::
Use_Semi_Lagrangian_Collidable_Advection(const T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T phi_replacement_value,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
{
    assert(!nested_semi_lagrangian_collidable&&!semi_lagrangian_collidable);
    nested_semi_lagrangian_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL(body_list,levelset->valid_mask_current,levelset->valid_mask_next,phi_replacement_value,true);
    semi_lagrangian_collidable=new ADVECTION_WRAPPER_COLLIDABLE_CELL<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>(*nested_semi_lagrangian_collidable,body_list,
        face_velocities_valid_mask_input);
    Set_Custom_Advection(*semi_lagrangian_collidable);
}
//#####################################################################
// Function HJ_WENO
// phi is (-2,m_3), phix_minus and phix_plus are (1,m)
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION<T_GRID>::
HJ_WENO(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus) const
{
    T epsilon=(T)1e-6; // works only because phi is a distance function
    T one_over_dx=1/dx;
    ARRAY<T,VECTOR<int,1> > D1(-2,m+2);for(int i=-2;i<=m+2;i++) D1(i)=(phi(i+1)-phi(i))*one_over_dx; // 1st divided difference
    for(int i=1;i<=m;i++){
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::WENO(D1(i-3),D1(i-2),D1(i-1),D1(i),D1(i+1),epsilon);
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::WENO(D1(i+2),D1(i+1),D1(i),D1(i-1),D1(i-2),epsilon);}
}
//#####################################################################
// Function HJ_ENO
// order = 1, 2 or 3, phi is (-2,m_3), phix_minus and phix_plus are (1,m)
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION<T_GRID>::
HJ_ENO(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus) const
{
    T one_over_dx=1/dx,one_over_two_dx=(T).5*one_over_dx,one_over_three_dx=(T)one_third*one_over_dx;
    ARRAY<T,VECTOR<int,1> > D1(-2,m+2),D2(-2,m+1),D3(-2,m); // divided differences
    for(int i=-2;i<=m+2;i++) D1(i)=(phi(i+1)-phi(i))*one_over_dx;
    if(order >= 2) for(int i=-2;i<=m+1;i++) D2(i)=(D1(i+1)-D1(i))*one_over_two_dx;
    if(order == 3) for(int i=-2;i<=m;i++) D3(i)=(D2(i+1)-D2(i))*one_over_three_dx;

    if(order == 1) for(int i=1;i<=m;i++){phix_minus(i)=D1(i-1);phix_plus(i)=D1(i);}
    else if(order == 2) for(int i=1;i<=m;i++){
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i-1),D2(i-2),D2(i-1));
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i),-D2(i),-D2(i-1));}
    else if(order == 3) for(int i=1;i<=m;i++){
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i-1),D2(i-2),D2(i-1),D3(i-3),D3(i-2),D3(i-1));
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i),-D2(i),-D2(i-1),D3(i),D3(i-1),D3(i-2));}
}
//#####################################################################
template class LEVELSET_ADVECTION<GRID<VECTOR<float,1> > >;
template class LEVELSET_ADVECTION<GRID<VECTOR<float,2> > >;
template class LEVELSET_ADVECTION<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class LEVELSET_ADVECTION<OCTREE_GRID<float> >;
template class LEVELSET_ADVECTION<QUADTREE_GRID<float> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class LEVELSET_ADVECTION<RLE_GRID_2D<float> >;
template class LEVELSET_ADVECTION<RLE_GRID_3D<float> >;
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION<GRID<VECTOR<double,1> > >;
template class LEVELSET_ADVECTION<GRID<VECTOR<double,2> > >;
template class LEVELSET_ADVECTION<GRID<VECTOR<double,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class LEVELSET_ADVECTION<OCTREE_GRID<double> >;
template class LEVELSET_ADVECTION<QUADTREE_GRID<double> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class LEVELSET_ADVECTION<RLE_GRID_2D<double> >;
template class LEVELSET_ADVECTION<RLE_GRID_3D<double> >;
#endif
#endif
