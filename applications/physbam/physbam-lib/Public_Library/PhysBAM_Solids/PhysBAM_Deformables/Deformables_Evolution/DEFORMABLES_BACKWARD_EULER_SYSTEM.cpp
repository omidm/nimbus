//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_BACKWARD_EULER_SYSTEM
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLES_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_EVOLUTION.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
DEFORMABLES_BACKWARD_EULER_SYSTEM(DEFORMABLES_EVOLUTION<TV>& deformables_evolution_input,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,const T dt,
    const T current_velocity_time,const T current_position_time,MPI_SOLIDS<TV>* mpi_solids_input,TRIANGLE_REPULSIONS<TV>* repulsions_input,const bool velocity_update_input)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),deformables_evolution(deformables_evolution_input),deformable_body_collection(deformable_body_collection),
    dt(dt),current_velocity_time(current_velocity_time),current_position_time(current_position_time),repulsions(repulsions_input),mpi_solids(mpi_solids_input),
    velocity_update(velocity_update_input),project_nullspace_frequency(INT_MAX),projection_data(deformable_body_collection)
{
    if(repulsions)
        repulsions->Set_Collision_Pairs(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,projection_data.point_face_pairs,
            projection_data.edge_edge_pairs,(T)1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
~DEFORMABLES_BACKWARD_EULER_SYSTEM()
{
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> void DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Force(const VECTOR_T& V,VECTOR_T& F) const
{
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(V.V.array);

    INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&> V_subset=F.V.array.Subset(deformable_body_collection.simulated_particles);
    ARRAYS_COMPUTATIONS::Fill(V_subset,TV());

    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(V.V.array);
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(V.V.array);
    deformable_body_collection.Implicit_Velocity_Independent_Forces(V.V.array,F.V.array,dt,current_velocity_time+dt);
    // we don't inherit velocity dependent forces to the drifted particles (contrary to the explicit velocity independent ones) because that would compromise symmetry
    deformable_body_collection.Add_Velocity_Dependent_Forces(V.V.array,F.V.array,current_velocity_time+dt);
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(F.V.array);
    deformable_body_collection.binding_list.Distribute_Force_To_Parents(F.V.array);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    Force(V,F);
    for(int i=1;i<=V.V.Size();i++) F.V(i)=V.V(i)-dt*projection_data.mass.one_over_mass(i)*F.V(i);    
}
//#####################################################################
// Function Set_Global_Boundary_Conditions
//#####################################################################
template<class TV> void DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Set_Global_Boundary_Conditions(VECTOR_T& V) const
{
    deformables_evolution.Set_External_Velocities(V.V.array,current_velocity_time+dt,current_position_time);
    deformables_evolution.kinematic_evolution.Set_External_Velocities(V.rigid_V.array,current_velocity_time+dt,current_position_time);
    if(velocity_update){
        PHYSBAM_ASSERT(ARRAY_VIEW<TV>::Same_Array(V.V.array,deformable_body_collection.particles.V));
        if(deformables_evolution.deformables_parameters.use_post_cg_constraints && repulsions) repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,false);}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
    VECTOR_T& V=debug_cast<VECTOR_T&>(BV);
    deformables_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.V.array,current_velocity_time+dt,current_position_time);
    for(int i=1;i<=deformables_evolution.deformables_parameters.implicit_solve_parameters.cg_projection_iterations;i++){
        int middle_projection=1;
        if(projection_data.point_face_precomputed.m || projection_data.edge_edge_precomputed.m){
            TRIANGLE_REPULSIONS<TV>::Project_All_Moving_Constraints(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,V.V.array);middle_projection=4;}
        if(middle_projection!=1) deformables_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.V.array,current_velocity_time+dt,current_position_time);}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(BV1),&V2=debug_cast<const VECTOR_T&>(BV2);
    double inner_product=ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(projection_data.mass.mass,V1.V,V2.V);
    if(mpi_solids) inner_product=mpi_solids->Reduce_Add(inner_product);
    return inner_product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(BR);
    T convergence_norm=ARRAYS_COMPUTATIONS::Maximum_Magnitude(R.V);
    if(mpi_solids) convergence_norm=mpi_solids->Reduce_Max(convergence_norm);
    return convergence_norm;
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
    if(++projection_data.project_nullspace_counter>=project_nullspace_frequency){Project(V);projection_data.project_nullspace_counter=0;}
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
    R.Copy((T)1,V);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_MASS<TV>::
DEFORMABLES_MASS(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection)
    :mass(deformable_body_collection.particles.mass,deformable_body_collection.dynamic_particles),
    one_over_mass(deformable_body_collection.particles.one_over_mass,deformable_body_collection.dynamic_particles)
{
}
//#####################################################################
// Function Inverse_Multiply
//#####################################################################
template<class TV> void DEFORMABLES_MASS<TV>::
Inverse_Multiply(const DEFORMABLES_VELOCITY<TV>& V,DEFORMABLES_VELOCITY<TV>& F) const
{
    F.V=one_over_mass*V.V;
}
//#####################################################################
template class DEFORMABLES_BACKWARD_EULER_SYSTEM<VECTOR<float,1> >;
template class DEFORMABLES_BACKWARD_EULER_SYSTEM<VECTOR<float,2> >;
template class DEFORMABLES_BACKWARD_EULER_SYSTEM<VECTOR<float,3> >;
template class DEFORMABLES_MASS<VECTOR<float,1> >;
template class DEFORMABLES_MASS<VECTOR<float,2> >;
template class DEFORMABLES_MASS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLES_BACKWARD_EULER_SYSTEM<VECTOR<double,1> >;
template class DEFORMABLES_BACKWARD_EULER_SYSTEM<VECTOR<double,2> >;
template class DEFORMABLES_BACKWARD_EULER_SYSTEM<VECTOR<double,3> >;
template class DEFORMABLES_MASS<VECTOR<double,1> >;
template class DEFORMABLES_MASS<VECTOR<double,2> >;
template class DEFORMABLES_MASS<VECTOR<double,3> >;
#endif
