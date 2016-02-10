//#####################################################################
// Copyright 2006, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/RIGID_BODY_FRACTURE_OBJECT_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION_S3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_S3D<T>::
Initialize_Bodies()
{
#if 0
    if(!fracture_object) PHYSBAM_FATAL_ERROR(); // this should be called from derived class after fracture objects are created

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    FINITE_VOLUME<TV,2>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
    TRIANGLE_BENDING_ELEMENTS<T>* bending_elements=solids_parameters.solid_body_collection.template Find_Force<TRIANGLE_BENDING_ELEMENTS<T>*>();

    // material properties save
    finite_volume.strain_measure.Initialize_Dm_Inverse_Save();
    finite_volume.Initialize_Be_Scales_Save();
    finite_volume.Save_V();
    if(bending_elements) bending_elements->Initialize_Reference_Quantities();

    // map fracture_bias_direction_to_rotated_space
    if(spatial_fracture_bias_direction.m != 0){
        for(int t=1;t<=fracture_object->fracture_bias_direction.m;t++){
            MATRIX<T,3,2> Q=finite_volume.strain_measure.Ds(deformable_body_collection.particles.X,t)*finite_volume.strain_measure.Dm_inverse(t);
            fracture_object->fracture_bias_direction(t)=Q.Transpose_Times(spatial_fracture_bias_direction(t));}
        spatial_fracture_bias_direction.Resize(0);}

    Reinitialize_Bodies();
    deformable_body_collection.deformable_geometry.template Find_Structure<T_EMBEDDED_MATERIAL_SURFACE&>().Create_Material_Surface();

    // initialize collisons (collisions must always be reinitialized after creating a material surface)
    Initialize_Self_Collision();
#endif
}
//#####################################################################
// Function Initialize_Self_Collision
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_S3D<T>::
Initialize_Self_Collision()
{
#if 0
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    if(solids_parameters.perform_self_collision || push_out){
        ARRAY<bool>& particle_on_surface=deformable_body_collection.collisions.check_collision;
        solids_parameters.solid_body_collection.deformable_body_collection.triangle_collisions.geometry.structures.Remove_All();
        solids_parameters.Initialize_Triangle_Collisions(); // false - do not clamp repulsion thickness
        ARRAY<TV> rest_X(deformable_body_collection.particles.array_collection->Size());
        for(int p=1;p<=rest_X.m;p++) if(particle_on_surface(p)) rest_X(p)=fracture_object->Rest_Position_Of_Material_Surface_Particle(p);
        solids_parameters.solid_body_collection.deformable_body_collection.triangle_repulsions.Clamp_Repulsion_Thickness_With_Meshes(rest_X,solids_parameters.collisions_repulsion_clamp_fraction);}
#endif
}

//#####################################################################
// Function Reinitialize_Bodies
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_S3D<T>::
Reinitialize_Bodies()
{
#if 0
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    FINITE_VOLUME<TV,2>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
    TRIANGLE_BENDING_ELEMENTS<T>* bending_elements=solids_parameters.solid_body_collection.template Find_Force<TRIANGLE_BENDING_ELEMENTS<T>*>();

    // copy back information
    finite_volume.strain_measure.Copy_Dm_Inverse_Save_Into_Dm_Inverse(fracture_object->corresponding_simplex_in_reference);
    finite_volume.Copy_Be_Scales_Save_Into_Be_Scales(fracture_object->corresponding_simplex_in_reference);
    if(bending_elements){
        bending_elements->Set_Quadruples_From_Reference_Triangle_Mesh(fracture_object->embedded_object.simplicial_object.mesh,
            fracture_object->corresponding_simplex_in_reference);
        bending_elements->Copy_Back_Reference_Quantities(fracture_object->corresponding_node_in_reference);}
    // resize collision arrays (Have to collide with the right particles !!!).
    deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,2>&>().Update_Number_Nodes();
    deformable_body_collection.collisions.Initialize_Object_Collisions(solids_parameters.collide_with_interior);

    // update particles and binding_list
    deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,2>&>().Update_Binding_List_From_Embedding(solids_parameters.solid_body_collection.deformable_body_collection);
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solids_parameters.solid_body_collection.deformable_body_collection.soft_bindings);

    solids_parameters.solid_body_collection.Update_Simulated_Particles();

    solids_parameters.solid_body_collection.Update_Position_Based_State(0);    // TODO: pass in correct time
#endif
}
//#####################################################################
// Function Rebuild_Topology
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_S3D<T>::
Rebuild_Topology()
{
#if 0

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    TRIANGLE_BENDING_ELEMENTS<T>* bending_elements=solids_parameters.solid_body_collection.template Find_Force<TRIANGLE_BENDING_ELEMENTS<T>*>();

    {std::stringstream ss;ss<<"number of elements = "<<fracture_object->embedded_object.simplicial_object.mesh.elements.m;
    {std::stringstream ss;ss<<"   ("<<fracture_object->embedded_object.Fraction_Of_Elements_With_Embedded_Subelements()<<")  ";
    {std::stringstream ss;ss<<"building embedded_object: ";

    if(bending_elements && bending_elements->plastic_yield) bending_elements->Initialize_Save_Quantities();
    ARRAY<int> map_to_old_particles,map_to_old_embedded_particles,map_to_old_simplices;
    fracture_object->Rebuild_Embedded_Object(map_to_old_particles,map_to_old_embedded_particles,map_to_old_simplices);
    if(plasticity_model) plasticity_model->Fp_inverse=ARRAY<SYMMETRIC_MATRIX<T,2> >(plasticity_model->Fp_inverse).Subset(map_to_old_simplices);
    if(bending_elements && bending_elements->plastic_yield) bending_elements->Copy_Back_Save_Quantities(map_to_old_particles);

    Reinitialize_Bodies();

    deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,2>&>().Create_Material_Surface();
    deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,2>&>().Perturb_Nodes_For_Collision_Freeness(perturb_amount_for_collision_freeness);

    // rebuild triangle collisions
    Initialize_Self_Collision();
#endif
}
//#####################################################################
// Function Fracture_Where_High_Stress
//#####################################################################
template<class T> int FRACTURE_EVOLUTION_S3D<T>::
Fracture_Where_High_Stress(const T small_number)
{
#if 0
    FINITE_VOLUME<TV,2>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
    {std::stringstream ss;ss << "computing sigma" << std::endl;LOG::filecout(ss.str());}

    // compute sigma and map bias directions
    ARRAY<SYMMETRIC_MATRIX<T,3> > sigma(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    ARRAY<TV> spatial_fracture_bias_direction(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    T max_eigenvalue_seen=0,min_eigenvalue_seen=0,max_sum=0,max_sum_sqr=0,min_sum=0,min_sum_sqr=0;
    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>& isotropic_model=dynamic_cast<ISOTROPIC_CONSTITUTIVE_MODEL<T,2>&>(finite_volume.constitutive_model);
    for(int t=1;t<=sigma.m;t++){
        DIAGONAL_MATRIX<T,2> Fe_hat_clipped=finite_volume.Fe_hat(t).Clamp_Min(small_number);
        T one_over_clipped_J=1/Fe_hat_clipped.Determinant();
        DIAGONAL_MATRIX<T,2> P_hat=isotropic_model.P_From_Strain(finite_volume.Fe_hat(t),(T)1,t); // scale for volume too?
        DIAGONAL_MATRIX<T,2> sigma_hat=one_over_clipped_J*P_hat*Fe_hat_clipped;
        sigma(t)=SYMMETRIC_MATRIX<T,3>::Conjugate(finite_volume.U(t),sigma_hat);
        max_eigenvalue_seen=max(sigma_hat.Max(),max_eigenvalue_seen);max_sum+=sigma_hat.Max();max_sum_sqr+=sqr(sigma_hat.Max());
        min_eigenvalue_seen=min(sigma_hat.Min(),min_eigenvalue_seen);min_sum+=sigma_hat.Min();min_sum_sqr+=sqr(sigma_hat.Min());
        spatial_fracture_bias_direction(t)=Spatial_Fracture_Bias_Direction(t,small_number);}

    // print out statistics
    T max_eigenvalue_mean=max_sum/sigma.m,min_eigenvalue_mean=min_sum/sigma.m;
    T max_eigenvalue_variance=max_sum_sqr/sigma.m-sqr(max_eigenvalue_mean),min_eigenvalue_variance=min_sum_sqr/sigma.m-sqr(min_eigenvalue_mean);
    if(max_eigenvalue_variance<0) max_eigenvalue_variance=0;if(min_eigenvalue_variance<0) min_eigenvalue_variance=0;
    T max_eigenvalue_standard_deviation=sqrt(max_eigenvalue_variance),min_eigenvalue_standard_deviation=sqrt(min_eigenvalue_variance);
    {std::stringstream ss;ss<<"max_eigenvalue_seen = "<<max_eigenvalue_seen<<"   ------   "<<"min_eigenvalue_seen = "<<min_eigenvalue_seen<<std::endl;LOG::filecout(ss.str());}
    {std::stringstream ss;ss<<"max:("<<max_eigenvalue_mean<<","<<max_eigenvalue_variance<<","<<max_eigenvalue_standard_deviation<<")"<<std::endl;LOG::filecout(ss.str());}
    {std::stringstream ss;ss<<"min:("<<min_eigenvalue_mean<<","<<min_eigenvalue_variance<<","<<min_eigenvalue_standard_deviation<<")"<<std::endl;LOG::filecout(ss.str());}

    // make cuts
    return fracture_object->Fracture_Where_High_Stress(sigma,spatial_fracture_bias_direction);
#endif
    return 0;
}
//#####################################################################
// Function Spatial_Fracture_Bias_Direction
//#####################################################################
template<class T> VECTOR<T,3> FRACTURE_EVOLUTION_S3D<T>::
Spatial_Fracture_Bias_Direction(const int t,const T small_number) const
{
#if 0
    FINITE_VOLUME<TV,2>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();

    TV transformed_fracture_bias_direction;
    int ref_t=fracture_object->corresponding_simplex_in_reference(t);
    if(plasticity_model){
        MATRIX<T,3,2> U;MATRIX<T,2> V;DIAGONAL_MATRIX<T,2>F_hat;finite_volume.strain_measure.F(t).Fast_Singular_Value_Decomposition(U,F_hat,V);
        transformed_fracture_bias_direction=U*(F_hat.Clamp_Min(small_number).Inverse()*V.Transpose_Times(fracture_object->fracture_bias_direction(ref_t)));}
    else{
        DIAGONAL_MATRIX<T,2> Fe_hat_inverse=finite_volume.Fe_hat(t).Clamp_Min(small_number).Inverse();
        MATRIX<T,3,2>& U=finite_volume.U(t);MATRIX<T,2>& V=(*finite_volume.V)(t);
        transformed_fracture_bias_direction=U*(Fe_hat_inverse*V.Transpose_Times(fracture_object->fracture_bias_direction(ref_t)));}
    return transformed_fracture_bias_direction;
#endif
    return TV();
}
//#####################################################################
template class FRACTURE_EVOLUTION_S3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_EVOLUTION_S3D<double>;
#endif

