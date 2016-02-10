//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_GRADIENT
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/SIDED_FACE_INDEX.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/BOUNDARY_CONDITION_INFO.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_FLUID_GRADIENT<TV>::
MATRIX_FLUID_GRADIENT(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :BASE(index_map_input),periodic_boundary(index_map_input.boundary_condition_collection.periodic_boundary)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_FLUID_GRADIENT<TV>::
~MATRIX_FLUID_GRADIENT()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT<TV>::
Compute(const ARRAY<bool,FACE_INDEX<d> >& psi_N)
{
    for(int axis=1;axis<=TV::dimension;++axis)
        if(periodic_boundary(axis) && (index_map.boundary_condition_collection.mpi_boundary(axis)(1) || index_map.boundary_condition_collection.mpi_boundary(axis)(2)))
            PHYSBAM_FATAL_ERROR("Periodic boundaries in MPI not supported in SPD solver");

    gradient.Reset(index_map.Number_Cells());
    VECTOR<TV,2> face_areas(-index_map.grid.Face_Sizes(),index_map.grid.Face_Sizes());
    const GRID<TV>& grid=index_map.grid;
    TV_INT grid_counts=grid.counts;
    
    for(int i=1;i<=index_map.indexed_faces.m;i++){
        const FACE_INDEX<d>& face_index=index_map.indexed_faces(i);
        if(!psi_N(face_index)){
            Add_Cell(i,face_index.axis,face_index.First_Cell_Index(),face_areas(1)(face_index.axis));
            Add_Cell(i,face_index.axis,face_index.Second_Cell_Index(),face_areas(2)(face_index.axis));
            Add_Interface(i,face_index,face_areas(1)(face_index.axis));}

        gradient.Finish_Row();}
    for(int i=1;i<=index_map.indexed_constraints.m;++i){
        const SIDED_FACE_INDEX<d>& face_index=index_map.indexed_constraints(i);
        Add_Cell(index_map.indexed_faces.m+i,face_index.axis,face_index.Real_Cell_Index(),face_areas(face_index.side)(face_index.axis));
        Add_Interface(index_map.indexed_faces.m+i,index_map.indexed_constraints(i),face_areas(1)(face_index.axis));

        gradient.Finish_Row();}
    gradient.Sort_Entries();
}
//#####################################################################
// Function Add_Cell
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT<TV>::
Add_Cell(int face_lookup,int axis,const TV_INT& cell_index,T weight)
{
    int cell_lookup;
    if(periodic_boundary(axis) && index_map.grid.Domain_Indices().Lazy_Outside(cell_index)){
        TV_INT periodic_offset_cell=cell_index;
        int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,index_map.grid.counts[axis]);
        periodic_offset_cell[axis]=axis_periodic_cell;
        cell_lookup=index_map.cell_indices(periodic_offset_cell);}
    else
        cell_lookup=index_map.cell_indices(cell_index);

    if(cell_lookup) gradient.Append_Entry_To_Current_Row(cell_lookup,weight);
    else{
        bool is_mpi_boundary=false;
        for(int axis=1;axis<=TV::dimension;++axis){
            if((cell_index(axis)<1 && index_map.boundary_condition_collection.mpi_boundary(axis)(1)) || (cell_index(axis)>index_map.grid.Counts()(axis) && index_map.boundary_condition_collection.mpi_boundary(axis)(2))) is_mpi_boundary=true; 
            if((cell_index(axis)<1 && !index_map.boundary_condition_collection.mpi_boundary(axis)(1)) || (cell_index(axis)>index_map.grid.Counts()(axis) && !index_map.boundary_condition_collection.mpi_boundary(axis)(2))){is_mpi_boundary=false;break;}}
        if(!is_mpi_boundary){
            typename BASE::GHOST_GRADIENT_ENTRY ge;
            ge.index=cell_index;
            ge.face=face_lookup;
            ge.weight=weight;
            ghost_gradient.Append(ge);}}
}
//#####################################################################
// Function Add_Interface
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT<TV>::
Add_Interface(int face_lookup,const FACE_INDEX<d>& face_index,T weight)
{
    TV_INT cell1=face_index.First_Cell_Index(),cell2=face_index.Second_Cell_Index();
    int ci1=index_map.cell_indices(cell1),ci2=index_map.cell_indices(cell2);
    const ARRAY<bool,TV_INT>& psi_D=index_map.boundary_condition_collection.psi_D;
    if(index_map.two_phase){
        if(!ci1 || !ci2) return;
        bool out1=psi_D(cell1),out2=psi_D(cell2);
        if(out1==out2) return;
        if(out2) weight=-weight;}
    else{
        if(!ci1==!ci2) return;
        if(!ci2) weight=-weight;}

    typename BASE::INTERFACE_ENTRY ie={face_lookup,weight};
    interface_gradient.Append(ie);
}
//#####################################################################
// Function Add_Interface
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT<TV>::
Add_Interface(int face_lookup,const SIDED_FACE_INDEX<d>& face_index,T weight)
{
    TV_INT cell=face_index.Real_Cell_Index();
    int ci=index_map.cell_indices(cell);
    const ARRAY<bool,TV_INT>& psi_D=index_map.boundary_condition_collection.psi_D;
    if(index_map.two_phase){
        if(!ci || psi_D(cell)) return;
        if(face_index.side==1) weight=-weight;}
    else{
        if(ci) return;
        if(face_index.side==1) weight=-weight;}

    typename BASE::INTERFACE_ENTRY ie={face_lookup,weight};
    interface_gradient.Append(ie);
}
template class MATRIX_FLUID_GRADIENT<VECTOR<float,1> >;
template class MATRIX_FLUID_GRADIENT<VECTOR<float,2> >;
template class MATRIX_FLUID_GRADIENT<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_FLUID_GRADIENT<VECTOR<double,1> >;
template class MATRIX_FLUID_GRADIENT<VECTOR<double,2> >;
template class MATRIX_FLUID_GRADIENT<VECTOR<double,3> >;
#endif
