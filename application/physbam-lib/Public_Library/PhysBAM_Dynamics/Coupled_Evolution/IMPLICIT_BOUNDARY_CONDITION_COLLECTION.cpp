//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>::
IMPLICIT_BOUNDARY_CONDITION_COLLECTION(BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input,bool set_all_neumann_cells_to_dirichlet_input,
    bool zero_all_dirichlet_face_velocities_input,bool use_psi_R_input,bool use_boundary_condition_info_input,TV_BOOL periodic_boundary_input)
    :callback(callback_input),set_all_neumann_cells_to_dirichlet(set_all_neumann_cells_to_dirichlet_input),zero_all_dirichlet_face_velocities(zero_all_dirichlet_face_velocities_input),
    use_psi_R(use_psi_R_input),use_boundary_condition_info(use_boundary_condition_info_input),periodic_boundary(periodic_boundary_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>::
~IMPLICIT_BOUNDARY_CONDITION_COLLECTION()
{
    boundary_conditions.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>::
Compute(const GRID<TV>& grid,ARRAY<T,TV_INT>& p,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time)
{
    psi_D.Resize(grid.Domain_Indices(1));
    psi_N.Resize(grid,1);
    if(use_psi_R) psi_R.Resize(grid,1);
    psi_D.Fill(false);
    psi_N.Fill(false);
    if(use_psi_R) psi_R.Fill((T)0);
    for(int i=1;i<=boundary_conditions.m;i++) boundary_conditions(i)->Update_Boundary_Conditions(grid,psi_D,psi_N,p,face_velocities,time);

    if(set_all_neumann_cells_to_dirichlet){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,1);iterator.Valid();iterator.Next()) if(All_Cell_Faces_Neumann(iterator.Cell_Index())){
            psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=0;}}

    if(zero_all_dirichlet_face_velocities){
        for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(grid);iterator.Valid();iterator.Next()){
            if(psi_D(iterator.First_Cell_Index()) && psi_D(iterator.Second_Cell_Index())){
                face_velocities(iterator.Axis(),iterator.Face_Index())=0;}}}

    if(use_boundary_condition_info) Compute_Boundary_Condition_Info(grid,p,face_velocities);
}
//#####################################################################
// Function Compute_Boundary_Condition_Info
//#####################################################################
template<class TV> void IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>::
Compute_Boundary_Condition_Info(const GRID<TV>& grid,const ARRAY<T,TV_INT>& p,const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities)
{
    boundary_condition_info.Remove_All();
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(grid);iterator.Valid();iterator.Next()){
        Compute_Boundary_Condition_Info(p,face_velocities,iterator.Full_Index(),1);
        Compute_Boundary_Condition_Info(p,face_velocities,iterator.Full_Index(),2);}
}
//#####################################################################
// Function Compute_Boundary_Condition_Info
//#####################################################################
template<class TV> void IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>::
Compute_Boundary_Condition_Info(const ARRAY<T,TV_INT>& p,const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const FACE_INDEX<TV::m>& f,int in_side)
{
    TV_INT a=f.Cell_Index(in_side),b=f.Cell_Index(3-in_side);
    if(psi_D(a) || !psi_D(b)) return;
    BOUNDARY_CONDITION_INFO<TV> bci;
    bci.in_side=in_side;
    bci.f=f;
    bci.value=0;
    bci.type=callback->Get_Boundary_Along_Ray(a,b,bci.theta,bci.value);
    if(bci.type==BOUNDARY_CONDITIONS_CALLBACKS<TV>::unused){
        if(psi_N(f)){bci.type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::noslip;bci.theta=(T).5;bci.value=face_velocities(f);}
        else{bci.type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::free;bci.theta=1;bci.value=p(b);}}
    boundary_condition_info.Append(bci);
}
//#####################################################################
// Function All_Cell_Faces_Neumann
//#####################################################################
template<class TV> bool IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>::
All_Cell_Faces_Neumann(const TV_INT& cell_index) const
{
    for(int axis=1;axis<=TV::dimension;axis++)
        if(!psi_N.Component(axis)(cell_index) || !psi_N.Component(axis)(cell_index+TV_INT::Axis_Vector(axis)))
            return false;
    return true;
}
template class IMPLICIT_BOUNDARY_CONDITION_COLLECTION<VECTOR<float,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLECTION<VECTOR<float,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLECTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMPLICIT_BOUNDARY_CONDITION_COLLECTION<VECTOR<double,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLECTION<VECTOR<double,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLECTION<VECTOR<double,3> >;
#endif
