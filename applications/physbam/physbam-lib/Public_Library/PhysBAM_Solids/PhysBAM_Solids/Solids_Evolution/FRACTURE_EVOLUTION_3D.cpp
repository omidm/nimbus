//#####################################################################
// Copyright 2006-2007, Frank Losasso, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/LEVELSET_GRAIN_BOUNDARIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/RIGID_BODY_FRACTURE_OBJECT_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/RIGID_FRACTURE_QUASISTATICS_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_3D<T>::
Initialize_Bodies()
{
#if 0
    if(!fracture_object) PHYSBAM_FATAL_ERROR(); // this should be called from derived class after fracture objects are created

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    FINITE_VOLUME<TV,3>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();

    // material properties save
    finite_volume.strain_measure.Initialize_Dm_Inverse_Save();
    finite_volume.Initialize_Be_Scales_Save();
    finite_volume.Save_V();

    // map fracture_bias_direction_to_rotated_space
    if(true/*spatial_fracture_bias_direction.m!=0*/){
        assert(fracture_object->fracture_bias_direction.m==fracture_object->reference_simplicial_object.mesh.elements.m);
        for(int t=1;t<=fracture_object->fracture_bias_direction.m;t++){
            MATRIX<T,3> Q=finite_volume.strain_measure.Ds(deformable_body_collection.particles.X,t)*finite_volume.strain_measure.Dm_inverse(t);Q.Transpose();
            fracture_object->fracture_bias_direction(t)=Q*fracture_object->fracture_bias_direction(t);}
    }//spatial_fracture_bias_direction.Resize(0);}

    Reinitialize_Bodies();
    deformable_body_collection.deformable_geometry.template Find_Structure<T_EMBEDDED_MATERIAL_SURFACE&>().Create_Material_Surface();

    // initialize collisons (collisions must always be reinitialized after creating a material surface)
    Initialize_Self_Collision();
#endif
}
//#####################################################################
// Function Initialize_Self_Collision
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_3D<T>::
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
template<class T> void FRACTURE_EVOLUTION_3D<T>::
Reinitialize_Bodies()
{
#if 0
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    FINITE_VOLUME<TV,3>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();

    // copy back information
    finite_volume.strain_measure.Copy_Dm_Inverse_Save_Into_Dm_Inverse(fracture_object->corresponding_simplex_in_reference);
    finite_volume.Copy_Be_Scales_Save_Into_Be_Scales(fracture_object->corresponding_simplex_in_reference);

    T_EMBEDDED_MATERIAL_SURFACE& embedding=deformable_body_collection.deformable_geometry.template Find_Structure<T_EMBEDDED_MATERIAL_SURFACE&>();
    // update particles and binding_list
    embedding.Update_Binding_List_From_Embedding(solids_parameters.solid_body_collection.deformable_body_collection);  // TAKEN OUT TO GET RIGID BODY AND DEFORMABLE TO ALIGN
    //solid_body_collection.deformable_body_collection.soft_bindings.Clean_Memory();  // TODO: reuse old soft bound particles
    // TODO: fix this to be same as solids standard tests
    //Substitute_Soft_Bindings_For_Embedded_Nodes(embedding,solid_body_collection.binding_list,solid_body_collection.deformable_body_collection.soft_bindings);
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solids_parameters.solid_body_collection.deformable_body_collection.soft_bindings);
    embedding.Update_Number_Nodes();embedding.material_surface.Update_Number_Nodes();
    deformable_body_collection.collisions.Initialize_Object_Collisions(solids_parameters.collide_with_interior);

    solids_parameters.solid_body_collection.Update_Simulated_Particles();

    solids_parameters.solid_body_collection.Update_Position_Based_State(0);    // TODO: pass in correct time
#endif
}
//#####################################################################
// Function Rebuild_Topology
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_3D<T>::
Rebuild_Topology()
{
#if 0
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;

    {std::stringstream ss;ss<<"number of elements = "<<fracture_object->embedded_object.simplicial_object.mesh.elements.m;
    {std::stringstream ss;ss<<"   ("<<fracture_object->embedded_object.Fraction_Of_Elements_With_Embedded_Subelements()<<")  ";
    {std::stringstream ss;ss<<"building embedded_object: ";
    int index;if(fracture_object->embedded_object.simplicial_object.Minimum_Signed_Volume(&index)<0)
        {std::stringstream ss;ss<<"MINIMUM SIGNED VOLUME IS NEGATIVE JUST BEFORE BULID CURRENT ETV\n\n"<<std::endl;LOG::filecout(ss.str());}

    ARRAY<int> map_to_old_particles,map_to_old_embedded_particles,map_to_old_simplices;
    PHYSBAM_DEBUG_PRINT("rebuild_topology",map_to_old_particles.m,map_to_old_embedded_particles.m,deformable_body_collection.particles.array_collection->Size(),(deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,3>&>()).particles.array_collection->Size());
    fracture_object->Rebuild_Embedded_Object(map_to_old_particles,map_to_old_embedded_particles,map_to_old_simplices,true);
    PHYSBAM_DEBUG_PRINT("rebuild_topology",map_to_old_particles.m,map_to_old_embedded_particles.m,deformable_body_collection.particles.array_collection->Size(),(deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,3>&>()).particles.array_collection->Size());
    fracture_object->embedded_object.simplicial_object.mesh.Initialize_Segment_Mesh();
    fracture_object->embedded_object.simplicial_object.mesh.Initialize_Element_Edges();
    if(plasticity_model) plasticity_model->Fp_inverse=ARRAY<SYMMETRIC_MATRIX<T,3> >(plasticity_model->Fp_inverse).Subset(map_to_old_simplices);

    {std::stringstream ss;ss<<"building boundary surface:";
    EMBEDDED_MATERIAL_SURFACE<TV,3>& embedding=deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,3>&>(); // TODO: can only fracture one object.
    embedding.Create_Material_Surface();embedding.Perturb_Nodes_For_Collision_Freeness(perturb_amount_for_collision_freeness);

    Reinitialize_Bodies();
    if(rigid_fracture && fractured_after_rebuild_topology && create_new_rigid_bodies){ // TODO: add this parameter for one time fracture.
        Create_New_Rigid_Bodies_From_Fracture(map_to_old_particles);
        solids_parameters.solid_body_collection.rigid_body_collection.rigid_geometry_collection.Destroy_Unreferenced_Geometry();
        //for(int i=1;i<=structure_list.Number_Of_Active_Elements();i++) structure_list.Element(structure_list.Active_Element_Id(i))->Update_Number_Nodes();
    }

    if(fracture_object->embedded_object.simplicial_object.Minimum_Signed_Volume(&index)<0)
        {std::stringstream ss;ss<<"MINIMUM SIGNED VOLUME IS NEGATIVE JUST AFTER BULID CURRENT ETV\n\n"<<std::endl;LOG::filecout(ss.str());}

    solids_parameters.solid_body_collection.Update_Simulated_Particles();
    embedding.Update_Number_Nodes();

    // rebuild triangle collisions
    Initialize_Self_Collision();
#endif
}
//#####################################################################
// Function Fracture_Where_High_Stress
//#####################################################################
template<class T> int FRACTURE_EVOLUTION_3D<T>::
Fracture_Where_High_Stress(const T small_number)
{
#if 0
    FINITE_VOLUME<TV,3>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
    {std::stringstream ss;ss<<"computing sigma"<<std::endl;LOG::filecout(ss.str());}

    // compute sigma and map bias directions
    ARRAY<SYMMETRIC_MATRIX<T,3> > sigma(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    ARRAY<TV> spatial_fracture_bias_direction(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    T max_eigenvalue_seen=0,min_eigenvalue_seen=0,max_sum=0,max_sum_sqr=0,min_sum=0,min_sum_sqr=0;

    ISOTROPIC_CONSTITUTIVE_MODEL<T,3>& isotropic_model=dynamic_cast<ISOTROPIC_CONSTITUTIVE_MODEL<T,3>&>(finite_volume.constitutive_model);
    for(int t=1;t<=sigma.m;t++){
        DIAGONAL_MATRIX<T,3> Fe_hat_clipped=finite_volume.Fe_hat(t).Clamp_Min(small_number);
        T one_over_clipped_J=1/Fe_hat_clipped.Determinant();
        DIAGONAL_MATRIX<T,3> P_hat=isotropic_model.P_From_Strain(finite_volume.Fe_hat(t),(T)1,t); // scale for volume too?
        DIAGONAL_MATRIX<T,3> sigma_hat=one_over_clipped_J*P_hat*Fe_hat_clipped;
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

    return fracture_object->Fracture_Where_High_Stress(sigma,spatial_fracture_bias_direction);
#endif
    return 0;
}
//#####################################################################
// Function Rigid_Fracture_Where_High_Stress
//#####################################################################
template<class T> int FRACTURE_EVOLUTION_3D<T>::
Rigid_Fracture_Where_High_Stress(const T small_number)
{
#if 0
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solids_parameters.solid_body_collection.rigid_body_collection;

    {std::stringstream ss;ss<<"computing sigma"<<std::endl;LOG::filecout(ss.str());}

    // compute sigma and map bias directions
    ARRAY<SYMMETRIC_MATRIX<T,3> > sigma(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    ARRAY<TV> spatial_fracture_bias_direction(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    T max_eigenvalue_seen=0,min_eigenvalue_seen=0,max_sum=0,max_sum_sqr=0,min_sum=0,min_sum_sqr=0;

    LINEAR_FINITE_VOLUME<TV,3> linear_finite_volume(fracture_object->embedded_object.simplicial_object,(T)1e5,(T).3,(T).02);
    for(int i=1;i<=rigid_bodies_with_impulse.m;i++){
        RIGID_BODY_FRACTURE_OBJECT_3D<T>& rigid_body_fracture_object=dynamic_cast<RIGID_BODY_FRACTURE_OBJECT_3D<T>&>(rigid_body_collection.Rigid_Body(rigid_bodies_with_impulse(i)));
        for(int t=1;t<=sigma.m;t++){
            // Compute Sigma From Quasistatics Computation
            sigma(t)=linear_finite_volume.Stress_Differential(rigid_body_fracture_object.average_dX,t);
            DIAGONAL_MATRIX<T,3> eigenvalues;MATRIX<T,3> eigenvectors;
            sigma(t).Solve_Eigenproblem(eigenvalues,eigenvectors);
            max_eigenvalue_seen=PhysBAM::max((T)eigenvalues.Max(),(T)max_eigenvalue_seen);max_sum+=eigenvalues.Max();max_sum_sqr+=sqr(eigenvalues.Max());
            min_eigenvalue_seen=PhysBAM::min((T)eigenvalues.Min(),(T)min_eigenvalue_seen);min_sum+=eigenvalues.Min();min_sum_sqr+=sqr(eigenvalues.Min());
            spatial_fracture_bias_direction(t)=Spatial_Fracture_Bias_Direction(t,small_number);}}

    // print out statistics
    T max_eigenvalue_mean=max_sum/sigma.m,min_eigenvalue_mean=min_sum/sigma.m;
    T max_eigenvalue_variance=max_sum_sqr/sigma.m-sqr(max_eigenvalue_mean),min_eigenvalue_variance=min_sum_sqr/sigma.m-sqr(min_eigenvalue_mean);
    if(max_eigenvalue_variance<0) max_eigenvalue_variance=0;if(min_eigenvalue_variance<0) min_eigenvalue_variance=0;
    T max_eigenvalue_standard_deviation=sqrt(max_eigenvalue_variance),min_eigenvalue_standard_deviation=sqrt(min_eigenvalue_variance);
    {std::stringstream ss;ss<<"max_eigenvalue_seen = "<<max_eigenvalue_seen<<"   ------   "<<"min_eigenvalue_seen = "<<min_eigenvalue_seen<<std::endl;LOG::filecout(ss.str());}
    {std::stringstream ss;ss<<"max:("<<max_eigenvalue_mean<<","<<max_eigenvalue_variance<<","<<max_eigenvalue_standard_deviation<<")"<<std::endl;LOG::filecout(ss.str());}
    {std::stringstream ss;ss<<"min:("<<min_eigenvalue_mean<<","<<min_eigenvalue_variance<<","<<min_eigenvalue_standard_deviation<<")"<<std::endl;LOG::filecout(ss.str());}

    int initial_number_of_embedded_subelements=fracture_object->embedded_object.embedded_object.mesh.elements.m;

    /* This is stress based fracture... theoretically correct, but region based fracture below seems to work better.
    for(int i=1;i<=rigid_bodies_with_impulse.m;i++){
        if(!fracture_object->embedded_object.simplicial_object.mesh.incident_elements) fracture_object->embedded_object.simplicial_object.mesh.Initialize_Incident_Elements();
        RIGID_BODY_FRACTURE_OBJECT_3D<T>* rigid_body=dynamic_cast<RIGID_BODY_FRACTURE_OBJECT_3D<T>*>(&rigid_body_collection.Rigid_Body(rigid_bodies_with_impulse(i)));
        SIMPLEX_MESH<3>& mesh=rigid_body->template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;

        if(rigid_body->levelset_grain_boundaries){
            for(int t=1;t<=mesh.elements.m;t++){
                int deformable_tet=rigid_body->rigid_to_deformable_tets(t);
                int number_of_cuts=fracture_object->embedded_object.Number_Of_Embedded_Cuts(deformable_tet);
                if(number_of_cuts<3){
                    VECTOR<int,4> regions;
                    int number_of_regions=rigid_body->levelset_grain_boundaries->Regions_Intersecting_Element(t,regions);
                    if(number_of_regions==1)continue;
                    T_DIAGONAL_MATRIX eigenvalues=sigma(deformable_tet).Fast_Eigenvalues();
                    T threshold=fracture_object->fracture_threshold[number_of_cuts+1]*rigid_body->levelset_grain_boundaries->Element_Weakness_Multiplier(number_of_regions,regions);
                    T amt_over=eigenvalues.First()-threshold,amt_under=fracture_object->compressive_threshold[number_of_cuts+1]-eigenvalues.Last();
                    if((amt_over>0 && amt_over>amt_under) || amt_under>0)for(int i=1;i<=number_of_regions-1;i++){
                        VECTOR<T,4> element_phi;rigid_body->levelset_grain_boundaries->Phi_For_Region_In_Element(t,regions[i],element_phi);
                        fracture_object->Add_Cut_Based_On_Phi(rigid_body->rigid_to_deformable_tets(t),element_phi);}}}}}
    fracture_object->Fracture_Where_High_Stress(sigma,spatial_fracture_bias_direction);*/

    for(int i=1;i<=rigid_bodies_with_impulse.m;i++){ // NOTE: this tries to make full cuts to break off full pieces
        if(!fracture_object->embedded_object.simplicial_object.mesh.incident_elements) fracture_object->embedded_object.simplicial_object.mesh.Initialize_Incident_Elements();
        RIGID_BODY_FRACTURE_OBJECT_3D<T>& rigid_body=dynamic_cast<RIGID_BODY_FRACTURE_OBJECT_3D<T>&>(rigid_body_collection.Rigid_Body(rigid_bodies_with_impulse(i)));
        SIMPLEX_MESH<3>& mesh=rigid_body.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
        if(rigid_body.grain_boundaries.m>0){
            for(int grain_boundary=1;grain_boundary<=rigid_body.grain_boundaries.m;grain_boundary++){
                ARRAY<TV> tet_centers(mesh.elements.m);
                for(int t=1;t<=mesh.elements.m;t++){
                    const VECTOR<int,4>& nodes=mesh.elements(t);
                    tet_centers(t)=rigid_body.Frame()*TETRAHEDRON<T>::Center(rigid_body.particles.X(nodes[1]),rigid_body.particles.X(nodes[2]),rigid_body.particles.X(nodes[3]),rigid_body.particles.X(nodes[4]));}
                rigid_body.grain_boundaries(grain_boundary)->Initialize_Element_Weakness_Multiplier(tet_centers);
                ARRAY<T> accumulator(rigid_body.grain_boundaries(grain_boundary)->Number_Of_Regions());
                ARRAY<int> number_of_nodes(rigid_body.grain_boundaries(grain_boundary)->Number_Of_Regions());
                for(int t=1;t<=mesh.elements.m;t++){
                    int deformable_tet=rigid_body.rigid_to_deformable_tets(t);
                    int number_of_cuts=fracture_object->embedded_object.Number_Of_Embedded_Cuts(deformable_tet);
                    VECTOR<int,4> regions;
                    int number_of_regions=rigid_body.grain_boundaries(grain_boundary)->Regions_Intersecting_Element(t,regions);
                    T_DIAGONAL_MATRIX eigenvalues=sigma(deformable_tet).Fast_Eigenvalues();
                    T threshold=fracture_object->fracture_threshold[min(number_of_cuts+1,3)]*rigid_body.grain_boundaries(grain_boundary)->Element_Weakness_Multiplier(t,number_of_regions,regions);
                    T amt_over=eigenvalues.First()-threshold;//amt_under=fracture_object->compressive_threshold[number_of_cuts+1]-eigenvalues.Last_Element();
                    for(int region=1;region<=number_of_regions;region++){
                        accumulator(regions(region))+=amt_over;number_of_nodes(regions(region))++;}}
                //TODO: why iis this commented out?
                //for(int region=1;region<=accumulator.m;region++)accumulator(region)/=number_of_nodes(region);
                rigid_body.grain_boundaries(grain_boundary)->Initialize_Breakability(tet_centers);
                for(int t=1;t<=mesh.elements.m;t++){
                    int deformable_tet=rigid_body.rigid_to_deformable_tets(t);
                    int number_of_cuts=fracture_object->embedded_object.Number_Of_Embedded_Cuts(deformable_tet);
                    if(number_of_cuts<3 && rigid_body.grain_boundaries(grain_boundary)->is_breakable(t)){
                        VECTOR<int,4> regions;
                        int number_of_regions=rigid_body.grain_boundaries(grain_boundary)->Regions_Intersecting_Element(t,regions);
                        if(number_of_regions==1) continue;
                        bool fracture=false;
                        // can't do region check for 
                        for(int region=1;region<=number_of_regions;region++)if(accumulator(regions(region))>0){fracture=true;break;}
                        if(fracture){
                            for(int i=1;i<=number_of_regions;i++){
                                VECTOR<T,4> element_phi;rigid_body.grain_boundaries(grain_boundary)->Phi_For_Region_In_Element(t,regions[i],element_phi);
                                fracture_object->Add_Cut_Based_On_Phi(rigid_body.rigid_to_deformable_tets(t),element_phi);
                                if(number_of_regions==2) break;}}}}}}}
    return fracture_object->embedded_object.embedded_object.mesh.elements.m-initial_number_of_embedded_subelements;
#endif
    return 0;
}
//#####################################################################
// Function Spatial_Fracture_Bias_Direction
//#####################################################################
template<class T> VECTOR<T,3> FRACTURE_EVOLUTION_3D<T>::
Spatial_Fracture_Bias_Direction(const int t,const T small_number) const
{
#if 0
    FINITE_VOLUME<TV,3>& finite_volume=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();

    TV transformed_fracture_bias_direction;
    int ref_t=fracture_object->corresponding_simplex_in_reference(t);
    if(plasticity_model){
        MATRIX<T,3>U,V;DIAGONAL_MATRIX<T,3>F_hat;finite_volume.strain_measure.F(t).Fast_Singular_Value_Decomposition(U,F_hat,V);
        transformed_fracture_bias_direction=U*(F_hat.Clamp_Min(small_number).Inverse()*V.Transpose_Times(fracture_object->fracture_bias_direction(ref_t)));}
    else{
        DIAGONAL_MATRIX<T,3> Fe_hat_inverse=finite_volume.Fe_hat(t).Clamp_Min(small_number).Inverse();
        MATRIX<T,3> &U=finite_volume.U(t),&V=(*finite_volume.V)(t);
        transformed_fracture_bias_direction=U*(Fe_hat_inverse*V.Transpose_Times(fracture_object->fracture_bias_direction(ref_t)));}
    return transformed_fracture_bias_direction;
#endif
    return TV();
}
//#####################################################################
// Function Adjust_Nodes_For_Segment_Triangle_Intersections (TODO: move to Deformable_Object_Collisions ?)
//#####################################################################
template<class T> int FRACTURE_EVOLUTION_3D<T>::
Adjust_Nodes_For_Segment_Triangle_Intersections(T threshhold)
{
#if 0
    int processed_nodes=0;

    TRIANGULATED_SURFACE<T>& material_surface=solids_parameters.solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<T_EMBEDDED_MATERIAL_SURFACE&>().material_surface;
    PARTICLES<TV >& particles=solids_parameters.solid_body_collection.deformable_body_collection.particles;
    ARRAY<VECTOR<int,2> > intersecting_segment_triangle_pairs;
    material_surface.Segment_Triangle_Intersection(*material_surface.mesh.segment_mesh,particles.X,solids_parameters.solid_body_collection.deformable_body_collection.triangle_collisions.geometry.small_number,true,
        &intersecting_segment_triangle_pairs);
    {std::stringstream ss;ss<<"!!!! found "<<intersecting_segment_triangle_pairs.m<<" intersecting edge triangle pairs"<<std::endl;LOG::filecout(ss.str());}

    // loop through intersecting edge triangle pairs
    ARRAY<TRIPLE<int,int,T> > edge_triangle_perturbations(particles.array_collection->Size());
    for(int pair=1;pair<=intersecting_segment_triangle_pairs.m;pair++){
        int e,t;intersecting_segment_triangle_pairs(pair).Get(e,t);
        VECTOR<int,2> segment=material_surface.mesh.segment_mesh->elements(e);
        int i,j,k;material_surface.mesh.elements(t).Get(i,j,k);
        TRIANGLE_3D<T> triangle(particles.X(i),particles.X(j),particles.X(k));
        T fraction;TV weights;
        if(INTERSECTION::Intersects(SEGMENT_3D<T>(particles.X(segment[1]),particles.X(segment[2])),triangle,fraction,weights,solids_parameters.solid_body_collection.deformable_body_collection.triangle_collisions.geometry.small_number)){
            for(int a=1;a<=2;a++){
                T perturb_amount=triangle.Signed_Distance(particles.X(segment[a]))+(T)1.5*solids_parameters.solid_body_collection.deformable_body_collection.triangle_collisions.collision_thickness;
                if(fraction<.5 && perturb_amount>0 && (!edge_triangle_perturbations(segment[a]).x || perturb_amount<edge_triangle_perturbations(segment[a]).z)){
                    edge_triangle_perturbations(segment[a])=TRIPLE<int,int,T>(e,t,perturb_amount);}}}}

    // sort by perturbation amount - always perturb with least amounts ones first - NOTE: testing with a VERY VERY primitive sort first.
    ARRAY<int> permutation(particles.array_collection->Size());for(int i=1;i<=permutation.m;i++) permutation(i)=i;
    Sort(permutation,Indirect_Comparison(edge_triangle_perturbations,Field_Comparison(&TRIPLE<int,int,T>::z))); // sort permutation by edge_triangle_perturbation.z

    // loop through edge_triangle_perturbation
    for(int i=1;i<=particles.array_collection->Size();i++){
        int p=permutation(i);
        TRIPLE<int,int,T>& edge_triangle_perturbation=edge_triangle_perturbations(p);
        if(edge_triangle_perturbation.x){
            int e=edge_triangle_perturbation.x,t=edge_triangle_perturbation.y;
            {std::stringstream ss;ss<<"p "<<p<<" with edge "<<e<<" triangle "<<t<<" and distance "<<edge_triangle_perturbation.z<<std::endl;LOG::filecout(ss.str());}
            if(edge_triangle_perturbation.z>1e-5){
                {std::stringstream ss;ss<<"skipping big perturbation here"<<std::endl;LOG::filecout(ss.str());}continue;}
            int a,b;material_surface.mesh.segment_mesh->elements(e).Get(a,b);
            int ti,tj,tk;material_surface.mesh.elements(t).Get(ti,tj,tk);
            TRIANGLE_3D<T> triangle(particles.X(ti),particles.X(tj),particles.X(tk));
            if(INTERSECTION::Intersects(SEGMENT_3D<T>(particles.X(a),particles.X(b)),triangle,solids_parameters.solid_body_collection.deformable_body_collection.triangle_collisions.geometry.small_number)){
                TV triangle_normal=triangle.Normal();
                particles.X(p)+=edge_triangle_perturbation.z*triangle_normal;
                T old_normal_velocity=TV::Dot_Product(particles.V(p),triangle_normal);
                TV weights=triangle.Barycentric_Coordinates(particles.X(p));
                TV triangle_velocity=TRIANGLE_3D<T>::Point_From_Barycentric_Coordinates(weights,particles.V(ti),particles.V(tj),particles.V(tk));
                T new_normal_velocity=TV::Dot_Product(triangle_velocity,triangle_normal);
                particles.V(p)+=(new_normal_velocity-old_normal_velocity)*triangle_normal;
                processed_nodes++;}
            else {std::stringstream ss;ss<<"no longer an intersection with node "<<p<<std::endl;LOG::filecout(ss.str());}}}

    return processed_nodes;
#endif
    return 0;
}
//#####################################################################
// Create_Rigid_Body_Fracture_Objects
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_3D<T>::
Create_Rigid_Body_Fracture_Object(const TV velocity,const TV angular_velocity,const ARRAY<TV>* seed_positions,const ARRAY<T>* seed_weakness_multipliers,const FRACTURE_CALLBACKS<TV>* fracture_callbacks)
{
#if 0
    T density=100;T uniform_levelset_cell_size=(T).05;int triangle_subdivision_loops=0;int levels_of_octree=0; // these can be passed in as parameters later; ditto for velocities
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,3>&>().embedded_object.simplicial_object;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solids_parameters.solid_body_collection.rigid_body_collection;
    RIGID_BODY_FRACTURE_OBJECT_3D<T>* rigid_body=new RIGID_BODY_FRACTURE_OBJECT_3D<T>(rigid_body_deformable_body_particles,solids_parameters.solid_body_collection,
        particle_to_rigid_body_id,deformable_to_rigid_particles,rigid_body_collection);
    rigid_body->Initialize_Rigid_Body_From_Fragment(density,uniform_levelset_cell_size,triangle_subdivision_loops,true,levels_of_octree);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);
    rigid_body->Update_Particle_To_Rigid_Body_Id_Mapping();
    rigid_body->Set_Coefficient_Of_Friction(1);
    rigid_body->Set_Coefficient_Of_Restitution(.5);
    rigid_body->Twist().linear=velocity;
    rigid_body->Twist().angular=angular_velocity;rigid_body->Update_Angular_Momentum();
    // TODO: need to find a way to do impulse body force
    //rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>* impulse_accumulator=new RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>(*rigid_body);
    TETRAHEDRALIZED_VOLUME<T>& rigid_body_tetrahedralized_volume=rigid_body->template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    rigid_body_tetrahedralized_volume.Initialize_Hierarchy(true);
    impulse_accumulator->Initialize_Simplicial_Object(&rigid_body_tetrahedralized_volume,new ARRAY<TV>(tetrahedralized_volume.particles.array_collection->Size()));
    rigid_body_collection.collision_body_list->Get_Collision_Geometry(rigid_body->particle_index)->impulse_accumulator=impulse_accumulator;
    if(fracture_callbacks && fracture_callbacks->Has_Crack_Map()) rigid_body->Initialize_Grain_Boundaries(*seed_positions,*seed_weakness_multipliers,fracture_callbacks,false);
    else if(seed_positions) rigid_body->Initialize_Grain_Boundaries(*seed_positions,*seed_weakness_multipliers,fracture_callbacks);
#endif
}
//#####################################################################
// Create_New_Rigid_Bodies_From_Fracture
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_3D<T>::
Create_New_Rigid_Bodies_From_Fracture(ARRAY<int>& map_to_old_particles)
{
#if 0
    {std::stringstream ss;ss<<"CREATE NEW BODIES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<std::endl;LOG::filecout(ss.str());}
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solids_parameters.solid_body_collection.deformable_body_collection;
    EMBEDDED_MATERIAL_SURFACE<TV,3>& embedding=deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,3>&>();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedding.embedded_object.simplicial_object;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solids_parameters.solid_body_collection.rigid_body_collection;
    particle_to_rigid_body_id.Resize(deformable_body_collection.particles.array_collection->Size());
    deformable_to_rigid_particles.Resize(deformable_body_collection.particles.array_collection->Size());
    for(int i=1;i<=rigid_bodies_with_impulse.m;i++){
        RIGID_BODY_FRACTURE_OBJECT_3D<T>& rigid_body=dynamic_cast<RIGID_BODY_FRACTURE_OBJECT_3D<T>&>(rigid_body_collection.Rigid_Body(rigid_bodies_with_impulse(i)));

        PHYSBAM_FATAL_ERROR("fix saved states, saved states are now driven through the fluid collision geometry");
        //T time_plus_dt=save_state?rigid_body.saved_states(COLLISION_GEOMETRY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_NEW_STATE)->time:0;
        T density=1;T uniform_levelset_cell_size=(T).05;int triangle_subdivision_loops=0;int levels_of_octree=0; // these should really come from the old rigid body
        // THIS SHOULD NOT MOVE THE CURRENT RIGID BODY. WE MAY NEED TO ALIGN FIRST OR SOMETHING SIMILAR. AFTER CREATION, BODY SHOULD BE AT OLD POSITION
        // IT WILL GET EVOLVED RIGHT AFTER THIS
        if(create_new_rigid_bodies){
            RIGID_BODY_FRACTURE_OBJECT_3D<T>* new_rigid_body=new RIGID_BODY_FRACTURE_OBJECT_3D<T>(rigid_body_deformable_body_particles,solids_parameters.solid_body_collection,
                particle_to_rigid_body_id,deformable_to_rigid_particles,rigid_body_collection);
            new_rigid_body->Initialize_Rigid_Body_From_Fragment(density,uniform_levelset_cell_size,triangle_subdivision_loops,true,levels_of_octree); // may need to remove extra stuff
            if(rigid_body.grain_boundaries.m>0){
                for(int grain_boundary=1;grain_boundary<=rigid_body.grain_boundaries.m;grain_boundary++){
                    ARRAY<TV> world_space_seed_positions(rigid_body.grain_boundaries(grain_boundary)->seed_positions.m);
                    for(int seed=1;seed<=world_space_seed_positions.m;seed++)world_space_seed_positions(seed)=rigid_body.Frame()*rigid_body.grain_boundaries(grain_boundary)->seed_positions(seed);
                    new_rigid_body->Initialize_Grain_Boundaries(world_space_seed_positions,rigid_body.grain_boundaries(grain_boundary)->seed_weakness_multipliers,
                        rigid_body.grain_boundaries(grain_boundary)->fracture_callbacks,rigid_body.grain_boundaries(grain_boundary)->is_levelset_grain_boundary);}}
            rigid_body_collection.Add_Rigid_Body_And_Geometry(new_rigid_body);
            new_rigid_body->Update_Particle_To_Rigid_Body_Id_Mapping();
            // TODO: fix adding of forces
            //new_rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
            PHYSBAM_FATAL_ERROR("TODO: fix adding of forces");
            new_rigid_body->parent_rigid_body_id=rigid_body.particle_index;
            new_rigid_body->Set_Coefficient_Of_Friction(rigid_body.coefficient_of_friction);
            new_rigid_body->Twist().linear=rigid_body.Pointwise_Object_Velocity(new_rigid_body->X());
            new_rigid_body->Twist().angular=rigid_body.Twist().angular;new_rigid_body->Update_Angular_Momentum();
            {std::stringstream ss;ss<<"+++++++++++++++ NEW RIGID BODY VELOCITY: "<<new_rigid_body->Twist().linear<<", angular: "<<new_rigid_body->Twist().angular<<std::endl;LOG::filecout(ss.str());}
            RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>* impulse_accumulator=new RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>(*new_rigid_body);
            TETRAHEDRALIZED_VOLUME<T>& new_rigid_body_tetrahedralized_volume=new_rigid_body->template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            new_rigid_body_tetrahedralized_volume.Initialize_Hierarchy(true);
            impulse_accumulator->Initialize_Simplicial_Object(&new_rigid_body_tetrahedralized_volume,new ARRAY<TV>(tetrahedralized_volume.particles.array_collection->Size()));

            rigid_body_collection.collision_body_list->Get_Collision_Geometry(new_rigid_body->particle_index)->impulse_accumulator=impulse_accumulator;
            PHYSBAM_FATAL_ERROR("fix saved states, saved states are now driven through the fluid collision geometry");
            //if(save_state) new_rigid_body->Save_State(COLLISION_GEOMETRY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_NEW_STATE,time_plus_dt);
            assert(new_rigid_body_tetrahedralized_volume.mesh.number_nodes==new_rigid_body_tetrahedralized_volume.particles.array_collection->Size());}
        rigid_body_collection.rigid_body_particle.Remove_Body(rigid_body.particle_index);}// may need it to be without destroying
    // geometry--depends on usage
#endif
}
//#####################################################################
// Function Substitute_Soft_Bindings_For_Embedded_Nodes
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_3D<T>::
Substitute_Soft_Bindings_For_Embedded_Nodes(T_EMBEDDED_MATERIAL_SURFACE& object,const BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,
    HASHTABLE<int,int>* persistent_soft_bindings)
{
#if 0
    ARRAY<int> nodes;object.material_surface.mesh.elements.Flattened().Get_Unique(nodes);
    ARRAY<int> map_to_new_particles(IDENTITY_ARRAY<>(object.material_surface.particles.array_collection->Size()));
    for(int i=1;i<=nodes.m;i++) if(binding_list.Binding_Index_From_Particle_Index(nodes(i))){
        int embedded_node=nodes(i),bound_node;
        if(persistent_soft_bindings && persistent_soft_bindings->Get(embedded_node,bound_node))
            persistent_soft_bindings->Delete(embedded_node);
        else{
            bound_node=object.material_surface.particles.array_collection->Add_Element_From_Deletion_List();
            object.material_surface.particles.array_collection->Copy_Element(*object.material_surface.particles.array_collection,embedded_node,bound_node);}
        map_to_new_particles(embedded_node)=bound_node;
        soft_bindings.Add_Binding(VECTOR<int,2>(bound_node,embedded_node),true);}
    for(int i=1;i<=object.material_surface.mesh.elements.m;i++) object.material_surface.mesh.elements(i)=VECTOR<int,3>::Map(map_to_new_particles,object.material_surface.mesh.elements(i));
#endif
}
//#####################################################################
// Function Process_Rigid_Fracture
//#####################################################################
template<class T> void FRACTURE_EVOLUTION_3D<T>::
Process_Rigid_Fracture(const T dt,const T time,SOLIDS_EVOLUTION<TV>* rigid_deformable_evolution_old,SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks)
{
#if 0
    rigid_bodies_with_impulse.Resize(0);

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solids_parameters.solid_body_collection.rigid_body_collection;
    for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i))
        if(RIGID_BODY_FRACTURE_OBJECT_3D<T>* rigid_body_fracture_object=dynamic_cast<RIGID_BODY_FRACTURE_OBJECT_3D<T>*>(&rigid_body_collection.Rigid_Body(i))){
            T threshold_squared=rigid_body_fracture_object->fracture_threshold*rigid_body_fracture_object->fracture_threshold;
            COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>* collision_body_impulse_accumulator=
                rigid_body_collection.collision_body_list->Get_Collision_Geometry(rigid_body_fracture_object->particle_index)->impulse_accumulator;
            RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>& impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>&>(*collision_body_impulse_accumulator);
            for(int rp=1;rp<=impulse_accumulator.accumulated_node_impulses->m;rp++)if((*impulse_accumulator.accumulated_node_impulses)(rp).Magnitude_Squared()>threshold_squared){
                rigid_bodies_with_impulse.Append(i);break;}}

    if(rigid_bodies_with_impulse.m){
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Before restore state",1,0);
        PHYSBAM_FATAL_ERROR("fix saved states, saved states are now driven through the fluid collision geometry");
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Before allign bodies",1,0);
        for(int i=1;i<=rigid_bodies_with_impulse.m;i++) dynamic_cast<RIGID_BODY_FRACTURE_OBJECT_3D<T>&>(rigid_body_collection.Rigid_Body(rigid_bodies_with_impulse(i))).Align_Deformable_Object_With_Rigid_Body();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("After allign bodies",1,0);


        bool need_to_simulate_rigid_bodies=Run_Quasistatics_And_Fracture(time,INT_MAX);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("After fracture and rebuild",1,0);

        if(need_to_simulate_rigid_bodies){
            LOG::Time("advancing rigid bodies in Postprocess_Fracture");
            PHYSBAM_FATAL_ERROR("fix saved states, saved states are now driven through the fluid collision geometry");
        }
    }
#endif
}
//#####################################################################
// Function Run_Quasistatics
//#####################################################################
template<class T> bool FRACTURE_EVOLUTION_3D<T>::
Run_Quasistatics_And_Fracture(const T time,const int max_number_of_fracture_iterations)
{
#if 0
    // run quasistatics
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solids_parameters.solid_body_collection.rigid_body_collection;
    for(int i=1;i<=rigid_bodies_with_impulse.m;i++){
        RIGID_BODY_FRACTURE_OBJECT_3D<T>& rigid_body_fracture_object=dynamic_cast<RIGID_BODY_FRACTURE_OBJECT_3D<T>&>(rigid_body_collection.Rigid_Body(rigid_bodies_with_impulse(i)));
        COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>* collision_body_impulse_accumulator=
            rigid_body_collection.collision_body_list->Get_Collision_Geometry(rigid_body_fracture_object.particle_index)->impulse_accumulator;
        RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>& impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>&>(*collision_body_impulse_accumulator);
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=rigid_body_fracture_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        SIMPLEX_MESH<3>& mesh=tetrahedralized_volume.mesh;mesh.Initialize_Neighbor_Nodes();
        assert(impulse_accumulator.accumulated_node_impulses->m==mesh.number_nodes && tetrahedralized_volume.particles.array_collection->Size()==mesh.number_nodes);
        ARRAY<int> distance_from_impulse(mesh.number_nodes);
        QUEUE<int> queue(mesh.number_nodes);
        for(int i=1;i<=mesh.number_nodes;i++)if((*impulse_accumulator.accumulated_node_impulses)(i).Magnitude_Squared()>0){
            distance_from_impulse(i)=1;queue.Enqueue(i);}
        while(!queue.Empty()){
            int node=queue.Dequeue();
            for(int neighbor=1;neighbor<=(*mesh.neighbor_nodes)(node).m;neighbor++){
                int neighbor_node=(*mesh.neighbor_nodes)(node)(neighbor);
                if(distance_from_impulse(neighbor_node)==0){distance_from_impulse(neighbor_node)=distance_from_impulse(node)+1;queue.Enqueue(neighbor_node);}}}

        // find three nodes to constrain (the ones furthest away from impulses)
        int max_node1,max_node2,max_node3;
        if(distance_from_impulse(1)>distance_from_impulse(2) && distance_from_impulse(1)>distance_from_impulse(3)){max_node1=1;
            if(distance_from_impulse(2)>distance_from_impulse(3)){max_node2=2;max_node3=3;}else{max_node2=3;max_node3=2;}}
        else if(distance_from_impulse(2)>distance_from_impulse(3)){max_node1=2;
            if(distance_from_impulse(1)>distance_from_impulse(3)){max_node2=1;max_node3=3;}else{max_node2=3;max_node3=1;}}
        else{max_node1=3;
            if(distance_from_impulse(1)>distance_from_impulse(2)){max_node2=1;max_node3=2;}else{max_node2=2;max_node3=1;}}

        for(int i=1;i<=mesh.number_nodes;i++){
            if(distance_from_impulse(i)>distance_from_impulse(max_node3)){
                if(distance_from_impulse(i)>distance_from_impulse(max_node2)){
                    if(distance_from_impulse(i)>distance_from_impulse(max_node1)){max_node3=max_node2;max_node2=max_node1;max_node1=i;}
                    else{max_node3=max_node2;max_node2=i;}}
                else{max_node3=i;}}}

        VECTOR<int,3> largest_nodes(max_node1,max_node2,max_node3);
        dynamic_cast<RIGID_FRACTURE_QUASISTATICS_FORCES<T>&>(*solids_parameters.solid_body_collection.example_forces_and_velocities).Initialize(rigid_body_fracture_object,largest_nodes);
        PARTICLES<TV>& particles=solids_parameters.solid_body_collection.deformable_body_collection.particles;
        ARRAY<TV>& average_dX=rigid_body_fracture_object.average_dX;ARRAYS_COMPUTATIONS::Fill(average_dX,TV());average_dX.Resize(particles.array_collection->Size());
        // TODO: should transform impulses to world space?
        QUASISTATIC_EVOLUTION<TV> quasistatic_evolution(solids_parameters);quasistatic_evolution.balance_external_forces_only=true;
        quasistatic_evolution.One_Newton_Step_Toward_Steady_State(time,average_dX);
        particles.X.Subset(rigid_body_fracture_object.rigid_to_deformable_particles)+=average_dX.Subset(rigid_body_fracture_object.rigid_to_deformable_particles);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("After showing averageDX",1,0);
        particles.X.Subset(rigid_body_fracture_object.rigid_to_deformable_particles)-=average_dX.Subset(rigid_body_fracture_object.rigid_to_deformable_particles);}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("After quasistatics",1,0);

    int number_of_fracture_iterations=0;
    while(Rigid_Fracture_Where_High_Stress()){
        fractured_after_rebuild_topology=true;
        number_of_fracture_iterations++;
        if(number_of_fracture_iterations>=max_number_of_fracture_iterations) break;}
    {std::stringstream ss;ss<<"fracture_iterations "<<number_of_fracture_iterations<<std::endl;LOG::filecout(ss.str());}
    bool need_to_simulate_rigid_bodies=fractured_after_rebuild_topology;
    if(fractured_after_rebuild_topology){
        Rebuild_Topology();fractured_after_rebuild_topology=false;}
    return need_to_simulate_rigid_bodies;
#endif
    return false;
}
//#####################################################################
template class FRACTURE_EVOLUTION_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_EVOLUTION_3D<double>;
#endif

