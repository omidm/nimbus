//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_INTERPOLATION_BASE
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
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION_BASE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_FLUID_INTERPOLATION_BASE<TV>::
MATRIX_FLUID_INTERPOLATION_BASE(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :index_map(index_map_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_FLUID_INTERPOLATION_BASE<TV>::
~MATRIX_FLUID_INTERPOLATION_BASE()
{
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_BASE<TV>::
Times(const VECTOR_ND<T>& faces,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const
{
    ARRAYS_COMPUTATIONS::Fill(constraints,T());
    Times_Add(faces,constraints);
}
//#####################################################################
// Function Transpose_Times
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_BASE<TV>::
Transpose_Times(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,VECTOR_ND<T>& faces) const
{
    faces.Fill(T());
    Transpose_Times_Add(constraints,faces);
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_BASE<TV>::
Test_Matrix() const
{
    RANDOM_NUMBERS<T> random;
    ARRAY<T,COUPLING_CONSTRAINT_ID> constraints(Number_Of_Constraints()),constraints2(Number_Of_Constraints());
    random.Fill_Uniform(constraints,-1,1);

    VECTOR_ND<T> faces(index_map.Number_Faces()),faces2(index_map.Number_Faces());
    random.Fill_Uniform(faces,-1,1);

    Times(faces,constraints2);
    Transpose_Times(constraints,faces2);

    T inner_faces=faces.Dot_Product(faces,faces2);
    T inner_constraints=ARRAYS_COMPUTATIONS::Dot_Product(constraints,constraints2);

    std::stringstream ss;
    ss<<"MATRIX_FLUID_INTERPOLATION_BASE Test: "<<inner_faces<<"  vs  "<<inner_constraints<<"  relative  "<<
        abs(inner_faces-inner_constraints)/maxabs((T)1e-30,inner_faces,inner_constraints)<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
template class MATRIX_FLUID_INTERPOLATION_BASE<VECTOR<float,1> >;
template class MATRIX_FLUID_INTERPOLATION_BASE<VECTOR<float,2> >;
template class MATRIX_FLUID_INTERPOLATION_BASE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_FLUID_INTERPOLATION_BASE<VECTOR<double,1> >;
template class MATRIX_FLUID_INTERPOLATION_BASE<VECTOR<double,2> >;
template class MATRIX_FLUID_INTERPOLATION_BASE<VECTOR<double,3> >;
#endif
