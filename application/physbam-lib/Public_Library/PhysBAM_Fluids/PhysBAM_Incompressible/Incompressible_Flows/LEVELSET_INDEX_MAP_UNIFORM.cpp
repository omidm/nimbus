//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_INDEX_MAP_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
using namespace PhysBAM;

template<class TV> LEVELSET_INDEX_MAP_UNIFORM<TV>::
//#####################################################################
// Constructor
//#####################################################################
LEVELSET_INDEX_MAP_UNIFORM(const GRID<TV>& grid_input,BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input)
    :grid(grid_input),callback(callback_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_INDEX_MAP_UNIFORM<TV>::
~LEVELSET_INDEX_MAP_UNIFORM()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void LEVELSET_INDEX_MAP_UNIFORM<TV>::
Compute(int axis,VECTOR<bool,d> periodic_boundary_input)
{
    inside.Resize(grid,1);
    face_to_index.Resize(grid);
    face_to_index.Fill(0);
    index_to_face.Remove_All();
    inside.Fill(false);
    callback->Mark_Outside(inside);
    periodic_boundary=periodic_boundary_input;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::WHOLE_REGION,0,axis);it.Valid();it.Next()){FACE_INDEX<d> index=it.Full_Index();
        bool& b=inside(index);
        b=!b;
        if(b && !(periodic_boundary[index.axis] && index.index[index.axis]==1))
            face_to_index(index)=index_to_face.Append(index);}

    for(int a=1;a<=d;a++) if(periodic_boundary[a])
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,2*a,a);it.Valid();it.Next()){
            FACE_INDEX<d> index=it.Full_Index();
            FACE_INDEX<d> matching_index(index);matching_index.index(a)=1;
            assert(face_to_index(index));
            face_to_index(matching_index)=face_to_index(index);}
}
//#####################################################################
// Function Gather
//#####################################################################
template<class TV> void LEVELSET_INDEX_MAP_UNIFORM<TV>::
Gather(const ARRAY<T,FACE_INDEX<d> >& u,VECTOR_ND<T>& v) const
{
    for(int i=1;i<=index_to_face.m;i++) v(i)=u(index_to_face(i));
}
//#####################################################################
// Function Scatter
//#####################################################################
template<class TV> void LEVELSET_INDEX_MAP_UNIFORM<TV>::
Scatter(const VECTOR_ND<T>& u,ARRAY<T,FACE_INDEX<d> >& v) const
{
    for(int i=1;i<=index_to_face.m;i++) v(index_to_face(i))=u(i);
    for(int a=1;a<=d;a++) if(periodic_boundary[a])
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,2*a,a);it.Valid();it.Next()){
            FACE_INDEX<d> index=it.Full_Index();
            FACE_INDEX<d> matching_index(index);matching_index.index(a)=1;
            v(matching_index)=v(index);}
}
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<float,1> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<float,2> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<double,1> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<double,2> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<double,3> >;
#endif
