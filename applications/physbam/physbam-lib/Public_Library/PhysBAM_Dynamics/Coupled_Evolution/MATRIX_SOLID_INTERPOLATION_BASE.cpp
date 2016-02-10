//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_INTERPOLATION_BASE
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
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_BASE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION_BASE<TV>::
MATRIX_SOLID_INTERPOLATION_BASE(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info)
    :iterator_info(info),V_size(0),rigid_V_size(0),V_indices(0),rigid_V_indices(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION_BASE<TV>::
~MATRIX_SOLID_INTERPOLATION_BASE()
{
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_BASE<TV>::
Times(const GENERALIZED_VELOCITY<TV>& solids,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const
{
    ARRAYS_COMPUTATIONS::Fill(constraints,T());
    Times_Add(solids,constraints);
}
//#####################################################################
// Function Transpose_Times
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_BASE<TV>::
Transpose_Times(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,GENERALIZED_VELOCITY<TV>& solids) const
{
    // TODO: Careful to zero out enough of the solids state.
    solids*=(T)0;
    Transpose_Times_Add(constraints,solids);
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_BASE<TV>::
Test_Matrix(int number_particles,int number_rigid_particles) const
{
    RANDOM_NUMBERS<T> random;

    ARRAY<TV> V(number_particles),V2(number_particles);
    random.Fill_Uniform(V,-1,1);

    ARRAY<TWIST<TV> > twist(number_rigid_particles),twist2(number_rigid_particles);
    random.Fill_Uniform(twist,-1,1);

    ARRAY<int> empty;
    GENERALIZED_VELOCITY<TV> solids(V,empty,twist,empty,empty),solids2(V2,empty,twist2,empty,empty);

    ARRAY<T,COUPLING_CONSTRAINT_ID> constraints(Number_Of_Constraints()),constraints2(Number_Of_Constraints());
    random.Fill_Uniform(constraints,-1,1);

    Times(solids,constraints2);
    Transpose_Times(constraints,solids2);

    CONSTANT_ARRAY<RIGID_BODY_MASS<TV,true> > rigid_mass(twist.m,RIGID_BODY_MASS<TV,true>(1,typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR()+1));
    T inner_solids=ARRAYS_COMPUTATIONS::Dot_Product(V,V2)+ARRAYS_COMPUTATIONS::Inner_Product(rigid_mass,twist,twist2);
    T inner_constraints=ARRAYS_COMPUTATIONS::Dot_Product(constraints,constraints2);

    std::stringstream ss;
    ss<<"MATRIX_SOLID_INTERPOLATION_BASE Test: "<<inner_solids<<"  vs  "<<inner_constraints<<"  relative  "<<
        abs(inner_solids-inner_constraints)/maxabs((T)1e-30,inner_solids,inner_constraints)<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Store_Maps
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_BASE<TV>::
Store_Maps(const GENERALIZED_VELOCITY<TV>& G)
{
    V_size=G.V.array.Size();
    rigid_V_size=G.rigid_V.Size();
    rigid_V_indices=&G.rigid_V.indices;
    V_indices=&G.V.indices;
}
//#####################################################################
template class MATRIX_SOLID_INTERPOLATION_BASE<VECTOR<float,1> >;
template class MATRIX_SOLID_INTERPOLATION_BASE<VECTOR<float,2> >;
template class MATRIX_SOLID_INTERPOLATION_BASE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_SOLID_INTERPOLATION_BASE<VECTOR<double,1> >;
template class MATRIX_SOLID_INTERPOLATION_BASE<VECTOR<double,2> >;
template class MATRIX_SOLID_INTERPOLATION_BASE<VECTOR<double,3> >;
#endif
