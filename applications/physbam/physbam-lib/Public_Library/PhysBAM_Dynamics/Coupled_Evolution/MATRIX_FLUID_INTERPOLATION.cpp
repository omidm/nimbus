//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_INTERPOLATION
//##################################################################### 
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_ELEMENT_ID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_FACE_INDEX.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/GENERALIZED_FLUID_MASS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_FLUID_INTERPOLATION<TV>::
MATRIX_FLUID_INTERPOLATION(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :MATRIX_FLUID_INTERPOLATION_BASE<TV>(index_map_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_FLUID_INTERPOLATION<TV>::
~MATRIX_FLUID_INTERPOLATION()
{
}
//#####################################################################
// Function Number_Of_Constraints
//#####################################################################
template<class TV> COUPLING_CONSTRAINT_ID MATRIX_FLUID_INTERPOLATION<TV>::
Number_Of_Constraints() const
{
    return rows.m;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION<TV>::
Compute(int ghost_cells)
{
    VECTOR<SIDED_FACE_INDEX<d>,2*d-2> neighbor_faces;
    rows.Remove_All();
    for(int i=1;i<=index_map.indexed_constraints.m;++i) rows.Append(index_map.indexed_faces.m+i);
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION<TV>::
Times_Add(const VECTOR_ND<T>& faces,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++)
        constraints(i)+=faces(rows(i));
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION<TV>::
Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,VECTOR_ND<T>& faces) const
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++)
        faces(rows(i))+=constraints(i);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION<TV>::
Print() const
{
    std::stringstream ss;
    for(int i=1;i<=index_map.indexed_faces.m;i++)
        ss<<i<<": "<<index_map.indexed_faces(i)<<std::endl;

    ARRAY<T,COUPLING_CONSTRAINT_ID> constraints(rows.m);
    VECTOR_ND<T> faces(index_map.Number_Faces());
    for(int i=1;i<=faces.n;i++){
        faces(i)=1;
        Times(faces,constraints);
        ss<<constraints<<std::endl;
        faces(i)=0;}
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION<TV>::
Print_Each_Matrix(int n) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("W-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("W",Value(rows.m),index_map.Number_Faces());

    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++)
        oo.Add_Sparse_Entry(Value(i),rows(i),1);

    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Add_Diagonal
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION<TV>::
Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_FLUID_MASS<TV>& fluid_mass) const
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++)
        diagonal(i)+=fluid_mass.one_over_fluid_mass_at_faces(rows(i));
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++)
        data.Append(TRIPLE<int,int,T>(Value(i),rows(i),1));
}
//#####################################################################
template class MATRIX_FLUID_INTERPOLATION<VECTOR<float,1> >;
template class MATRIX_FLUID_INTERPOLATION<VECTOR<float,2> >;
template class MATRIX_FLUID_INTERPOLATION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_FLUID_INTERPOLATION<VECTOR<double,1> >;
template class MATRIX_FLUID_INTERPOLATION<VECTOR<double,2> >;
template class MATRIX_FLUID_INTERPOLATION<VECTOR<double,3> >;
#endif
