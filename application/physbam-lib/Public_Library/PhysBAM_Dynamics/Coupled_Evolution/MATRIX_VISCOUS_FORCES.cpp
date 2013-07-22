//#####################################################################
// Copyright 2010, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_VISCOUS_FORCES
//##################################################################### 
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_FACE_INDEX.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_VISCOUS_FORCES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_VISCOUS_FORCES<TV>::
MATRIX_VISCOUS_FORCES(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :grid(index_map_input.grid),index_map(index_map_input)
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Compute(const T dt,const ARRAY<bool,FACE_INDEX<d> >& psi_N,T mu)
{
    // Setting up the correct gradient
    // TODO: Dirchlet cells & Neumann faces
    entries.Remove_All();

    T coefficient=sqrt(dt*mu*grid.One_Over_Cell_Size());
    TV face_areas=grid.Face_Sizes();
    last_id=VISCOUS_FORCE_ID();
    if(!mu) return;

    for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();FACE_INDEX<d> this_face(iterator.Full_Index());
        for(int other_axis=1;other_axis<=TV::dimension;other_axis++){
            last_id++;
            FACE_INDEX<d> other_face(axis,face_index+TV_INT::Axis_Vector(other_axis));T weight=coefficient*face_areas(axis);
            int this_index(0),other_index(0);
            if(index_map.face_indices.Valid_Index(this_face)) this_index=index_map.face_indices(this_face);
            else if(index_map.constraint_indices.Contains(SIDED_FACE_INDEX<TV::dimension>(2,this_face))) this_index=index_map.indexed_faces.m+index_map.constraint_indices.Get(SIDED_FACE_INDEX<TV::dimension>(2,this_face));
            if(index_map.face_indices.Valid_Index(other_face)) other_index=index_map.face_indices(other_face);
            else if(index_map.constraint_indices.Contains(SIDED_FACE_INDEX<TV::dimension>(1,other_face))) other_index=index_map.indexed_faces.m+index_map.constraint_indices.Get(SIDED_FACE_INDEX<TV::dimension>(1,other_face));
            if((this_index>index_map.indexed_faces.m && other_index>index_map.indexed_faces.m) || !this_index || !other_index){last_id--;continue;}
            if(this_index) entries.Append(ENTRY(-weight,this_index,last_id));
            if(other_index) entries.Append(ENTRY(weight,other_index,last_id));}}
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Times_Add(const VECTOR_ND<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const
{
    for(int i=1;i<=entries.m;i++){
        const ENTRY& entry=entries(i);
        viscous_force_coefficients(entry.viscous_id)+=entry.weight*velocities(entry.face_index);}
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Times(const VECTOR_ND<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const
{
    ARRAYS_COMPUTATIONS::Fill(viscous_force_coefficients,T());
    Times_Add(velocities,viscous_force_coefficients);
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Transpose_Times_Add(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,VECTOR_ND<T>& velocities) const
{
    for(int i=1;i<=entries.m;i++){
        const ENTRY& entry=entries(i);
        velocities(entry.face_index)+=entry.weight*viscous_force_coefficients(entry.viscous_id);}
}
//#####################################################################
// Function Transpose_Times
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Transpose_Times(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,VECTOR_ND<T>& velocities) const
{
    //ARRAYS_COMPUTATIONS::Fill(velocities,T());
    velocities.Fill(T());
    Transpose_Times_Add(viscous_force_coefficients,velocities);
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
// can depend on position too
template<class TV> VISCOUS_FORCE_ID MATRIX_VISCOUS_FORCES<TV>::
Viscous_Forces_Size() const
{
    return last_id;
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Test_Matrix() const
{
    RANDOM_NUMBERS<T> random;
    ARRAY<T,VISCOUS_FORCE_ID> K(last_id),K2(K);
    VECTOR_ND<T> V(index_map.Number_Faces()),V2(V);

    random.Fill_Uniform(K,-1,1);
    random.Fill_Uniform(V,-1,1);
    Times(V,K2);
    Transpose_Times(K,V2);

    T inner_V=V.Dot_Product(V,V2);
    T inner_K=ARRAYS_COMPUTATIONS::Dot_Product(K,K2);
    std::stringstream ss;
    ss<<"MATRIX_VISCOUS_FORCES Symmetry Test: "<<inner_V<<"  vs  "<<inner_K<<"  relative  "<<abs(inner_V-inner_K)/maxabs((T)1e-30,inner_V,inner_K)<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Print_Each_Matrix(int n) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("N-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("N",Value(last_id),index_map.Number_Faces());
    for(int i=1;i<=entries.m;i++){
        const ENTRY& entry=entries(i);
        oo.Add_Sparse_Entry(Value(entry.viscous_id),entry.face_index,entry.weight);}
    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    for(int i=1;i<=entries.m;i++){
        const ENTRY& entry=entries(i);
        data.Append(TRIPLE<int,int,T>(Value(entry.viscous_id),entry.face_index,entry.weight));}
}
//#####################################################################
template class MATRIX_VISCOUS_FORCES<VECTOR<float,1> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<float,2> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_VISCOUS_FORCES<VECTOR<double,1> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<double,2> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<double,3> >;
#endif
