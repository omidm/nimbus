//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_POISSON
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/SIDED_FACE_INDEX.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_POISSON.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_FLUID_POISSON<TV>::
MATRIX_FLUID_POISSON(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input,
    const T_ARRAYS_SCALAR& one_over_rho_c_squared_input)
    :index_map(index_map_input),one_over_rho_c_squared(one_over_rho_c_squared_input),dt(0)
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_FLUID_POISSON<TV>::
Compute(const SPARSE_MATRIX_FLAT_MXN<T>& gradient,const VECTOR_ND<T>& one_over_fluid_mass,const T dt_in,const bool use_preconditioner)
{
    if(!use_preconditioner) return;
    if(!dt_in){
        SPARSE_MATRIX_FLAT_MXN<T> negative_divergence;
        gradient.Transpose(negative_divergence);
        poisson=negative_divergence.Times_Diagonal_Times(one_over_fluid_mass,gradient).Create_NXN_Matrix();}

    // will preconditioning work correctly with this matrix?
    Compute_Preconditioner(dt_in);
}
//#####################################################################
// Function Compute_Preconditioner
//#####################################################################
template<class TV> void MATRIX_FLUID_POISSON<TV>::
Compute_Preconditioner(const T dt_in)
{
    dt=dt_in;
    if(dt){const T dt_squared_by_V=(dt*dt)/index_map.grid.Cell_Size();
        V_over_rho_c_squared_dt_squared_inverse_flat.Resize(index_map.Number_Cells());
        for(int i=1;i<=index_map.indexed_cells.m;i++){TV_INT cell_index=index_map.indexed_cells(i);
            V_over_rho_c_squared_dt_squared_inverse_flat(i)= // TODO(jontg): save and pass in rho_c_squared
                dt_squared_by_V/one_over_rho_c_squared(cell_index);}
    return;}
#if 1
    ARRAY<int> rmap(poisson.n);
    for(int i=1;i<=poisson.n;i++)
        if(poisson.offsets(i)!=poisson.offsets(i+1))
            rmap(i)=map.Append(i);

    delete poisson.C;
    poisson.C=new SPARSE_MATRIX_FLAT_NXN<T>;
    poisson.C->Reset();
    int index=poisson.offsets(1);
    for(int i=1;i<=poisson.n;i++){
        int end=poisson.offsets(i+1);
        if(index==end) continue;
        for(;index<end;index++)
            poisson.C->Append_Entry_To_Current_Row(rmap(poisson.A(index).j),poisson.A(index).a);
        poisson.C->Finish_Row();}

    poisson.C->In_Place_Incomplete_Cholesky_Factorization();
#else
    ARRAY<int> row_lengths(index_map.real_cell_indices.m);
    for(int row=1;row<=index_map.real_cell_indices.m;++row){int row_index=index_map.real_cell_indices(row);
        for(int index=poisson.offsets(row_index);index<poisson.offsets(row_index+1);++index){
            int col_index=poisson.A(index).j,col=index_map.real_cell_indices_reverse_map(col_index);
            if(col >= 0){row_lengths(col)++;}}}

    delete poisson.C;poisson.C=new SPARSE_MATRIX_FLAT_NXN<T>;
    poisson.C->Set_Row_Lengths(row_lengths);

    for(int row=1;row<=index_map.real_cell_indices.m;++row){int row_index=index_map.real_cell_indices(row);
        for(int index=poisson.offsets(row_index);index<poisson.offsets(row_index+1);++index){
            int col_index=poisson.A(index).j,col=index_map.real_cell_indices_reverse_map(col_index);
            if(col >= 0) poisson.C->Set_Element(row,col,poisson.A(index).a);}}

    poisson.C->In_Place_Incomplete_Cholesky_Factorization();
#endif
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MATRIX_FLUID_POISSON<TV>::
Apply_Preconditioner(VECTOR_ND<T>& pressure) const
{
    if(dt){for(int i=1;i<=index_map.indexed_cells.m;i++) pressure(i)*=V_over_rho_c_squared_dt_squared_inverse_flat(i);
        return;}
#if 1
    VECTOR_ND<T> sub_vector(map.m),temp_vector(map.m);
    for(int i=1;i<=map.m;++i) sub_vector(i)=pressure(map(i));

    poisson.C->Solve_Forward_Substitution(sub_vector,temp_vector,true);
    poisson.C->Solve_Backward_Substitution(temp_vector,sub_vector,false,true);

    for(int i=1;i<=map.m;++i) pressure(map(i))=sub_vector(i);
#else
    VECTOR_ND<T> sub_vector(index_map.real_cell_indices.m),temp_vector(index_map.real_cell_indices.m);
    for(int i=1;i<=index_map.real_cell_indices.m;++i) sub_vector(i)=pressure(index_map.real_cell_indices(i));

    poisson.C->Solve_Forward_Substitution(sub_vector,temp_vector,true);
    poisson.C->Solve_Backward_Substitution(temp_vector,sub_vector,false,true);

    for(int i=1;i<=index_map.real_cell_indices.m;++i) pressure(index_map.real_cell_indices(i))=sub_vector(i);
#endif
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_FLUID_POISSON<TV>::
Times_Add(const VECTOR_ND<T>& pressure_in,VECTOR_ND<T>& pressure_out) const
{
    VECTOR_ND<T> result(pressure_in.n);
    poisson.Times(pressure_in,result);
    pressure_out+=result;
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_FLUID_POISSON<TV>::
Times(const VECTOR_ND<T>& pressure_in,VECTOR_ND<T>& pressure_out) const
{
    poisson.Times(pressure_in,pressure_out);
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_POISSON<TV>::
Test_Matrix(const bool print_matrix) const
{
    RANDOM_NUMBERS<T> random;
    VECTOR_ND<T> a(poisson.n),a2(poisson.n),b(poisson.n),b2(poisson.n);
    random.Fill_Uniform(a,-1,1);
    random.Fill_Uniform(b,-1,1);

    poisson.Times(b,a2);
    poisson.Times(a,b2);

    T inner1=VECTOR_ND<T>::Dot_Product(a,a2);
    T inner2=VECTOR_ND<T>::Dot_Product(b,b2);
    std::stringstream ss;
    ss<<"MATRIX_FLUID_POISSON Test: "<<inner1<<"  vs  "<<inner2<<"  relative  "<<abs(inner1-inner2)/maxabs((T)1e-30,inner1,inner2)<<std::endl;

    if(print_matrix) ss<<"poisson:\n"<<poisson<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_POISSON<TV>::
Print_Each_Matrix(int n) const
{
    OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("Poiss-%i.txt",n).c_str()).Write("Poiss",poisson);
}
//#####################################################################
template class MATRIX_FLUID_POISSON<VECTOR<float,1> >;
template class MATRIX_FLUID_POISSON<VECTOR<float,2> >;
template class MATRIX_FLUID_POISSON<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_FLUID_POISSON<VECTOR<double,1> >;
template class MATRIX_FLUID_POISSON<VECTOR<double,2> >;
template class MATRIX_FLUID_POISSON<VECTOR<double,3> >;
#endif
