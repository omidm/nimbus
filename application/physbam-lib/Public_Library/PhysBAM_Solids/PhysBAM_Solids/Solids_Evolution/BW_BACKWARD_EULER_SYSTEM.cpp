//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_BACKWARD_EULER_SYSTEM
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/BW_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/ASYNCHRONOUS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BW_BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BW_BACKWARD_EULER_SYSTEM<TV>::
BW_BACKWARD_EULER_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection,BW_COLLISIONS<TV>& bw_collisions_input,const T dt_input,const T time_input)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),solid_body_collection(solid_body_collection),bw_collisions(bw_collisions_input),dt(dt_input),time(time_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BW_BACKWARD_EULER_SYSTEM<TV>::
~BW_BACKWARD_EULER_SYSTEM()
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
    Project_Helper(V,true);
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Force(const VECTOR_T& V,VECTOR_T& F) const
{
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&> V_subset=F.V.array.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> >,ARRAY<int>&> rigid_V_subset=F.rigid_V.array.Subset(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles);
    ARRAYS_COMPUTATIONS::Fill(V_subset,TV());ARRAYS_COMPUTATIONS::Fill(rigid_V_subset,TWIST<TV>());
    solid_body_collection.Implicit_Velocity_Independent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,dt,time);
    solid_body_collection.Add_Velocity_Dependent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,time);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    Force(V,F);
    for(int i=1;i<=V.V.Size();i++) F.V(i)=V.V(i)-dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(i)*F.V(i);
/*    for(int i=1;i<=V.rigid_V.Size();i++) F.rigid_V(i)=V.rigid_V(i)-dt*(projection_data.mass.world_space_rigid_mass_inverse(i)*F.rigid_V(i));*/
}
//#####################################################################
// Function Project_Helper
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Project_Helper(KRYLOV_VECTOR_BASE<T>& BV,const bool negate) const
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    // User constrained nodes
    VECTOR_T& V=debug_cast<VECTOR_T&>(BV);
    solid_body_collection.example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(V.V.array,time,time);

    // Cloth/body contacts
    for(int i=1;i<=bw_collisions.cloth_body_constraints.m;i++){
        int particle_index=bw_collisions.cloth_body_constraints(i).x;
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(bw_collisions.cloth_body_constraints(i).y);
        TV particle_location=particles.X(particle_index);
        TV normal=rigid_body.Implicit_Geometry_Normal(particle_location);
        TV relative_velocity=rigid_body.Pointwise_Object_Velocity(particle_location)-particles.V(particle_index);
        if(negate) V.V.array(particle_index)=TV::Dot_Product(relative_velocity,normal)*normal;
        else V.V.array(particle_index)=(MATRIX<T,TV::m>::Identity_Matrix()-MATRIX<T,TV::m>::Outer_Product(normal,normal))*V.V.array(particle_index);
    }
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
    Project_Helper(BV,false);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double BW_BACKWARD_EULER_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(BV1),&V2=debug_cast<const VECTOR_T&>(BV2);
    double inner_product=ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(solid_body_collection.deformable_body_collection.particles.mass,V1.V,V2.V);
//+ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(projection_data.mass.world_space_rigid_mass,V1.rigid_V,V2.rigid_V);
    //TODO this should not be mass scaled
    return inner_product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR BW_BACKWARD_EULER_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(BR);
    T convergence_norm_squared=sqr(ARRAYS_COMPUTATIONS::Maximum_Magnitude(R.V));
/*    for(int p=1;p<=R.rigid_V.Size();p++){
        const TWIST<TV>& twist=R.rigid_V(p);
        const RIGID_BODY_MASS<TV,true> &rigid_mass=projection_data.mass.world_space_rigid_mass(p),&rigid_mass_inverse=projection_data.mass.world_space_rigid_mass_inverse(p);
        convergence_norm_squared=max(convergence_norm_squared,
        twist.linear.Magnitude_Squared()+rigid_mass_inverse.mass*Dot_Product(twist.angular,rigid_mass.inertia_tensor*twist.angular));}*/
    T convergence_norm=sqrt(convergence_norm_squared);
    return convergence_norm;
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
}
//#####################################################################
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<float,1> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<float,2> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<double,1> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<double,2> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<double,3> >;
#endif
