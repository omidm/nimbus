//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_FACE_INDEX.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_FACE_POISSON_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_INDEX_MAP_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_FACE_POISSON_UNIFORM<TV>::
LEVELSET_FACE_POISSON_UNIFORM(LEVELSET_INDEX_MAP_UNIFORM<TV>& index_map_input)
    :index_map(index_map_input),theta_threshold((T)1e-5)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_FACE_POISSON_UNIFORM<TV>::
~LEVELSET_FACE_POISSON_UNIFORM()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void LEVELSET_FACE_POISSON_UNIFORM<TV>::
Compute(int axis)
{
    TV one_over_dX2=index_map.grid.one_over_dX*index_map.grid.one_over_dX;
    P.Reset(index_map.index_to_face.m);
    b.Resize(index_map.index_to_face.m);
    b.Set_Zero();
    typedef BOUNDARY_CONDITIONS_CALLBACKS<TV> CB;
    for(int i=1;i<=index_map.index_to_face.m;i++){FACE_INDEX<d> face=index_map.index_to_face(i);
        T middle_value=-2*one_over_dX2.Sum();
        for(int a=1;a<=d;a++){
            for(int s=-1;s<=1;s+=2){
                FACE_INDEX<d> neighbor=face;
                neighbor.index(a)+=s;
                if(index_map.inside(neighbor)){P.Append_Entry_To_Current_Row(index_map.face_to_index(neighbor),one_over_dX2(a));continue;}

                T value=0,theta=0;
                switch(index_map.callback->Get_Boundary_Along_Ray(face,neighbor,theta,value)){
                    case CB::noslip:if(theta<theta_threshold) theta=theta_threshold;middle_value+=one_over_dX2(a)*(1-1/theta);b(i)+=one_over_dX2(a)*value/theta;break;
                    case CB::slip:case CB::free:middle_value+=one_over_dX2(a)*(1+value*index_map.grid.dX(a));break;
                    default: PHYSBAM_FATAL_ERROR("Expected boundary condition");break;}}}
        P.Append_Entry_To_Current_Row(i,middle_value);
        P.Finish_Row();}
    P.Sort_Entries();
}
template class LEVELSET_FACE_POISSON_UNIFORM<VECTOR<float,1> >;
template class LEVELSET_FACE_POISSON_UNIFORM<VECTOR<float,2> >;
template class LEVELSET_FACE_POISSON_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_FACE_POISSON_UNIFORM<VECTOR<double,1> >;
template class LEVELSET_FACE_POISSON_UNIFORM<VECTOR<double,2> >;
template class LEVELSET_FACE_POISSON_UNIFORM<VECTOR<double,3> >;
#endif
