//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED
//##################################################################### 
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info)
    :MATRIX_SOLID_INTERPOLATION_BASE<TV>(info)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
~MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED()
{
}
//#####################################################################
// Function Number_Of_Constraints
//#####################################################################
template<class TV> COUPLING_CONSTRAINT_ID MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Number_Of_Constraints() const
{
    return COUPLING_CONSTRAINT_ID(entries.m*TV::m);
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Compute(const int ghost_cells)
{
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Times_Add(const GENERALIZED_VELOCITY<TV>& solids,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const
{
    for(int i=1;i<=entries.m;i++) for(int axis=1;axis<=TV::m;axis++)
        constraints(COUPLING_CONSTRAINT_ID((i-1)*TV::m+axis))+=solids.V.array(entries(i))(axis);
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,GENERALIZED_VELOCITY<TV>& solids) const
{
    for(int i=1;i<=entries.m;i++) for(int axis=1;axis<=TV::m;axis++)
        solids.V.array(entries(i))(axis)+=constraints(COUPLING_CONSTRAINT_ID((i-1)*TV::m+axis));
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Print_Each_Matrix(int n,GENERALIZED_VELOCITY<TV>& G) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("J-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("J",Value(Number_Of_Constraints()),G.Raw_Size());
    ARRAY<int> reverse_map_deformable(G.V.array.Size());
    reverse_map_deformable.Subset(G.V.indices)=IDENTITY_ARRAY<>(G.V.Size());

    for(int i=1;i<=entries.m;i++)
        for(int axis=1;axis<=TV::m;axis++)
            oo.Add_Sparse_Entry((i-1)*TV::m+axis,(reverse_map_deformable(entries(i))-1)*TV::m+axis,1);

    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    ARRAY<int> reverse_map_deformable(this->V_size);
    reverse_map_deformable.Subset(*this->V_indices)=IDENTITY_ARRAY<>(this->V_size);

    for(int i=1;i<=entries.m;i++)
        for(int axis=1;axis<=TV::m;axis++)
            data.Append(TRIPLE<int,int,T>((i-1)*TV::m+axis,(reverse_map_deformable(entries(i))-1)*TV::m+axis,1));
}
//#####################################################################
// Function Add_Diagonal
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_MASS<TV>& solid_mass) const
{
    for(int i=1;i<=entries.m;i++) for(int axis=1;axis<=TV::m;axis++)
        diagonal(COUPLING_CONSTRAINT_ID((i-1)*TV::m+axis))+=solid_mass.one_over_mass.array(entries(i));
}
//#####################################################################
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,1> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,2> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,1> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,2> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,3> >;
#endif
