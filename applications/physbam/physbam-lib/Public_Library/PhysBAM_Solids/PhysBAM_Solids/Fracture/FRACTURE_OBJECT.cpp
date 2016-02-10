//#####################################################################
// Copyright 2006, Geoffrey Irving, Frank Losasso, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_OBJECT
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/LEVELSET_GRAIN_BOUNDARIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/VIRTUAL_NODE_ALGORITHM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> FRACTURE_OBJECT<TV,d>::
FRACTURE_OBJECT(T_EMBEDDED_OBJECT& embedded_object_input)
    :eigenvector_coefficient((T)1),fracture_bias_direction_coefficient(0),fracture_bias_propagation_coefficient(0),number_of_fracture_initiations(0),
    embedded_object(embedded_object_input),reference_simplicial_object(reference_mesh,reference_particles),max_number_of_cuts(d),
    number_of_smoothing_passes(4),number_of_fracture_initiation_points(1),bias_stress(true),
    fracture_quality_threshold(0),extra_levelset_propagation_direction_scaling((T)100),force_edge_connected_fracture(true),
    initiation_point_positions(0),initiation_point_reference_seed_positions(0),initiation_point_radii(0)
{
    random_numbers.Set_Seed(84396);
    reference_mesh.Initialize_Mesh(embedded_object.simplicial_object.mesh);
    reference_particles.array_collection->Initialize(*embedded_object.particles.array_collection);
    Set_Fracture_Bias_Propagation();Set_Fracture_Bias_Stress_Scaling();Set_Fracture_Bias_Magnitude();
    fracture_bias_direction.Resize(embedded_object.simplicial_object.mesh.elements.m);
    fracture_phi_index.Resize(embedded_object.simplicial_object.mesh.elements.m);
    embedded_object.Initialize_Node_In_Simplex_Is_Material();
    corresponding_node_in_reference=IDENTITY_ARRAY<>(embedded_object.particles.array_collection->Size());
    corresponding_simplex_in_reference=IDENTITY_ARRAY<>(embedded_object.simplicial_object.mesh.elements.m);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> FRACTURE_OBJECT<TV,d>::
~FRACTURE_OBJECT()
{}
//#####################################################################
// Function Rebuild_Embedded_Object
//#####################################################################
template<class TV,int d> void FRACTURE_OBJECT<TV,d>::
Rebuild_Embedded_Object(ARRAY<int>& map_to_old_particles,ARRAY<int>& map_to_old_embedded_particles,ARRAY<int>& map_to_old_simplices,const bool verbose)
{
    VIRTUAL_NODE_ALGORITHM::Rebuild_Embedded_Object(embedded_object,map_to_old_particles,map_to_old_embedded_particles,map_to_old_simplices,verbose);
    corresponding_node_in_reference.Resize(map_to_old_particles.m);
    corresponding_node_in_reference=ARRAY<int>(corresponding_node_in_reference).Subset(map_to_old_particles);
    corresponding_simplex_in_reference=ARRAY<int>(corresponding_simplex_in_reference).Subset(map_to_old_simplices);
}
//#####################################################################
// Function Rest_Position_Of_Material_Surface_Particle
//#####################################################################
template<class TV,int d> TV FRACTURE_OBJECT<TV,d>::
Rest_Position_Of_Material_Surface_Particle(const int material_surface_particle)
{
    if(int embedded_material_surface_particle=embedded_object.embedded_particles.subset_index_from_point_cloud_index(material_surface_particle)){
        T lambda=embedded_object.interpolation_fraction(embedded_material_surface_particle);
        VECTOR<int,2> ref_parents(corresponding_node_in_reference.Subset(embedded_object.parent_particles(embedded_material_surface_particle)));
        return LINEAR_INTERPOLATION<T,TV>::Linear(reference_particles.X(ref_parents[1]),reference_particles.X(ref_parents[2]),lambda);}
    else{
        int reference_particle=corresponding_node_in_reference(material_surface_particle);
        return reference_particles.X(reference_particle);}
}
//#####################################################################
// Function Embedded_Subelement_Normal
//#####################################################################
template<class TV> static TV
Embedded_Segment_Vector(const EMBEDDED_OBJECT<TV,2>& embedded_object,const int emb_seg,const int current_triangle)
{
    assert(embedded_object.Element_Containing_Subelement(emb_seg)==current_triangle);
    int a,b;VECTOR<int,2>::Map(embedded_object.embedded_particles.subset_index_from_point_cloud_index,embedded_object.embedded_mesh.elements(emb_seg)).Get(a,b);
    return embedded_object.Position_Of_Embedded_Particle(b)-embedded_object.Position_Of_Embedded_Particle(a);
}
// returns a vector normal to emb_seg and the triangle_normal
template<class T> static VECTOR<T,2>
Embedded_Subelement_Normal(const EMBEDDED_OBJECT<VECTOR<T,2>,2>& embedded_object,const int emb_seg,const int current_triangle)
{
    return Embedded_Segment_Vector(embedded_object,emb_seg,current_triangle).Rotate_Clockwise_90().Normalized();
}
template<class T> static VECTOR<T,3>
Embedded_Subelement_Normal(const EMBEDDED_OBJECT<VECTOR<T,3>,2>& embedded_object,const int emb_seg,const int current_triangle)
{
    int i,j,k;embedded_object.simplicial_object.mesh.elements(current_triangle).Get(i,j,k);
    VECTOR<T,3> triangle_normal=TRIANGLE_3D<T>::Normal_Direction(embedded_object.particles.X(i),embedded_object.particles.X(j),embedded_object.particles.X(k));
    return VECTOR<T,3>::Cross_Product(triangle_normal,Embedded_Segment_Vector(embedded_object,emb_seg,current_triangle)).Normalized();
}
template<class T> static VECTOR<T,3>
Embedded_Subelement_Normal(const EMBEDDED_OBJECT<VECTOR<T,3>,3>& embedded_object,const int emb_tri,const int tetrahedron)
{
    assert(embedded_object.Element_Containing_Subelement(emb_tri)==tetrahedron);
    int a,b,c;VECTOR<int,3>::Map(embedded_object.embedded_particles.subset_index_from_point_cloud_index,embedded_object.embedded_object.mesh.elements(emb_tri)).Get(a,b,c);
    return TRIANGLE_3D<T>::Normal(embedded_object.Position_Of_Embedded_Particle(a),embedded_object.Position_Of_Embedded_Particle(b),embedded_object.Position_Of_Embedded_Particle(c));
}
//#####################################################################
// Function Project_Onto_Simplex_Span
//#####################################################################
template<class T> static void Project_Onto_Simplex_Span(const TRIANGULATED_SURFACE<T>& triangulated_surface,const int triangle,VECTOR<T,3>& fracture_normal)
{
    int i,j,k;triangulated_surface.mesh.elements(triangle).Get(i,j,k);
    VECTOR<T,3> triangle_normal=TRIANGLE_3D<T>::Normal(triangulated_surface.particles.X(i),triangulated_surface.particles.X(j),triangulated_surface.particles.X(k));
    fracture_normal-=VECTOR<T,3>::Dot_Product(fracture_normal,triangle_normal)*triangle_normal;
}
template<class T> static void Project_Onto_Simplex_Span(const TRIANGULATED_AREA<T>&,const int,VECTOR<T,2>&){}
template<class T> static void Project_Onto_Simplex_Span(const TETRAHEDRALIZED_VOLUME<T>&,const int,VECTOR<T,3>&){}
//#####################################################################
// Function Fracture_Where_High_Stress
//#####################################################################
template<class TV,int d> int FRACTURE_OBJECT<TV,d>::
Fracture_Where_High_Stress(ARRAY<T_SYMMETRIC_MATRIX>& sigma,ARRAY<TV>& spatial_fracture_bias_direction)
{
    int initial_number_of_embedded_subelements=embedded_object.embedded_object.mesh.elements.m;
    if(number_of_smoothing_passes || (bias_stress && fracture_bias_propagation)) embedded_object.simplicial_object.mesh.Initialize_Adjacent_Elements();
    if(!embedded_object.simplicial_object.mesh.incident_elements) embedded_object.simplicial_object.mesh.Initialize_Incident_Elements();
    if(!embedded_object.embedded_mesh.incident_elements) embedded_object.embedded_mesh.Initialize_Incident_Elements();
    if(!embedded_object.embedded_subelements_in_parent_element) embedded_object.Initialize_Embedded_Subelements_In_Parent_Element();
    embedded_object.Initialize_Parents_To_Embedded_Particles_Hash_Table();

    // smooth sigma
    {std::stringstream ss;ss<<"smoothing sigma. Smoothing passes: "<<number_of_smoothing_passes<<std::endl;LOG::filecout(ss.str());}
    if(number_of_smoothing_passes){
        ARRAY<T> size(embedded_object.simplicial_object.mesh.elements.m);
        for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++) size(t)=embedded_object.simplicial_object.Element_Size(t);
        for(int pass=1;pass<=number_of_smoothing_passes;pass++) for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++){
            T total_size=size(t);
            sigma(t)*=size(t);
            for(int a=1;a<=(*embedded_object.simplicial_object.mesh.adjacent_elements)(t).m;a++){
                int adj_t=(*embedded_object.simplicial_object.mesh.adjacent_elements)(t)(a);
                sigma(t)+=size(adj_t)*sigma(adj_t);total_size+=size(adj_t);}
            sigma(t)/=total_size;}}

    // bias stress
    {std::stringstream ss;ss<<"biasing stress: "<<bias_stress<<std::endl;LOG::filecout(ss.str());}
    ARRAY<TV> propagation_direction(embedded_object.simplicial_object.mesh.elements.m);
    if(bias_stress) for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++){
        int ref_t=corresponding_simplex_in_reference(t);
        if(embedded_object.Number_Of_Embedded_Cuts(t)!=0)continue;
        if(fracture_bias_stress_scaling(ref_t)!=1) sigma(t)*=fracture_bias_stress_scaling(ref_t);
        if(fracture_bias_magnitude(ref_t) && fracture_bias_direction(ref_t)!=VECTOR<T,d>())
            sigma(t)+=fracture_bias_magnitude(ref_t)*T_SYMMETRIC_MATRIX::Outer_Product(spatial_fracture_bias_direction(t));
        if(fracture_bias_propagation){ // bias stress based on cuts in adjacent elements
            for(int a=1;a<=(*embedded_object.simplicial_object.mesh.adjacent_elements)(t).m;a++){
                int adj_t=(*embedded_object.simplicial_object.mesh.adjacent_elements)(t)(a);
                VECTOR<int,2*d-2> emb_elements=embedded_object.Embedded_Subelements_In_Element(adj_t);
                for(int i=1;i<=emb_elements.m && emb_elements[i];i++){
                    int isolated_node=embedded_object.Node_Separated_By_Embedded_Subelement(emb_elements[i]);
                    if(!isolated_node || embedded_object.simplicial_object.mesh.Node_In_Simplex(isolated_node,t)){
                        TV emb_face_normal=Embedded_Subelement_Normal(embedded_object,emb_elements[i],adj_t);
                        propagation_direction(t)+=sign(TV::Dot_Product(emb_face_normal,propagation_direction(t)))*emb_face_normal;}}}
            T magnitude=propagation_direction(t).Magnitude();
            if(magnitude){
                if(Fracture_Phi_Index(t)) propagation_direction(t)*=extra_levelset_propagation_direction_scaling;
                propagation_direction(t)/=magnitude;sigma(t)+=fracture_bias_propagation*T_SYMMETRIC_MATRIX::Outer_Product(propagation_direction(t));}}}

    // make cuts
    assert(max_number_of_cuts<=d);int count=0;
    for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++){
        T_DIAGONAL_MATRIX eigenvalues=sigma(t).Fast_Eigenvalues();
        int number_of_cuts=embedded_object.Number_Of_Embedded_Cuts(t);
        if(number_of_cuts<max_number_of_cuts){ // room for more cuts
            int cut=number_of_cuts+1;
            // check Rankine condition
            TV fracture_normal;
            T threshold=fracture_threshold[cut];
            T amt_over=eigenvalues.First()-threshold,amt_under=compressive_threshold[cut]-eigenvalues.Last();
            if(amt_over>0 && amt_over>amt_under) fracture_normal=sigma(t).First_Eigenvector_From_Ordered_Eigenvalues(eigenvalues);
            else if(amt_under>0) fracture_normal=sigma(t).Last_Eigenvector_From_Ordered_Eigenvalues(eigenvalues);
            else continue;
            // add cuts
            count++;
            if(cut==1){
                int ref_t=corresponding_simplex_in_reference(t);
                fracture_normal=eigenvector_coefficient*fracture_normal+fracture_bias_propagation_coefficient*propagation_direction(t);
                if(fracture_bias_direction_coefficient && fracture_bias_direction(ref_t)!=VECTOR<T,d>()) fracture_normal+=fracture_bias_direction_coefficient*spatial_fracture_bias_direction(t);}
            fracture_normal.Normalize();
            Project_Onto_Simplex_Span(embedded_object.simplicial_object,t,fracture_normal);
            Add_Cut(t,fracture_normal);}}

    {std::stringstream ss;ss<<"fraction above threshold="<<(T)count/(T)embedded_object.simplicial_object.mesh.elements.m<<std::endl;LOG::filecout(ss.str());}
    return embedded_object.embedded_object.mesh.elements.m-initial_number_of_embedded_subelements;
}
//#####################################################################
// Function Set_Random_Fracture_Bias_Stress_Scaling_Constant
//#####################################################################
template<class TV,int d> void FRACTURE_OBJECT<TV,d>::
Set_Random_Fracture_Bias_Stress_Scaling_Constant(const T fracture_threshold,const int averaging_iterations) // TODO: this function is apparently deprecated
{
    bool adjacent_elements_defined=(embedded_object.simplicial_object.mesh.adjacent_elements!=0);
    if(!adjacent_elements_defined) embedded_object.simplicial_object.mesh.Initialize_Adjacent_Elements();

    T fracture_threshold_over_two=(T).5*fracture_threshold;
    for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++)
        fracture_bias_stress_scaling(t)=random_numbers.Get_Uniform_Number((T)-fracture_threshold_over_two,(T)fracture_threshold_over_two);
    for(int iteration=1;iteration<=averaging_iterations;iteration++) for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++){
        int number_of_adjacent_elements=(*embedded_object.simplicial_object.mesh.adjacent_elements)(t).m;
        T average=fracture_bias_stress_scaling(t);
        for(int i=1;i<=number_of_adjacent_elements;i++){
            int adjacent_t=(*embedded_object.simplicial_object.mesh.adjacent_elements)(t)(i);
            average+=fracture_bias_stress_scaling(adjacent_t);}
        average/=(T)(number_of_adjacent_elements+1);
        fracture_bias_stress_scaling(t)=average;}

    if(!adjacent_elements_defined){
        delete embedded_object.simplicial_object.mesh.adjacent_elements;
        embedded_object.simplicial_object.mesh.adjacent_elements=0;}
}
//#####################################################################
// Function Initiation_Point
//#####################################################################
template<class TV,int d> bool FRACTURE_OBJECT<TV,d>::
Initiation_Point(const int element)
{
    if(number_of_fracture_initiations>=number_of_fracture_initiation_points) return false;
    if(embedded_object.Number_Of_Edges_With_Embedded_Particles(element)) return false;
    if(initiation_point_positions){
        assert(initiation_point_positions->m==number_of_fracture_initiations);
        TV reference_centroid=reference_simplicial_object.Centroid(corresponding_simplex_in_reference(element));
        for(int i=1;i<=number_of_fracture_initiations;i++) if(((*initiation_point_positions)(i)-reference_centroid).Magnitude()<(*initiation_point_radii)(i)) return false;
        if(!initiation_point_reference_seed_positions) initiation_point_positions->Append(reference_centroid);
        else if(number_of_fracture_initiations){int i=number_of_fracture_initiations+1;
            if(((*initiation_point_reference_seed_positions)(i)-reference_centroid).Magnitude()<(*initiation_point_radii)(i)) return false;
            initiation_point_positions->Append((*initiation_point_reference_seed_positions)(i));}}
    number_of_fracture_initiations++;return true;
}
//#####################################################################
template class FRACTURE_OBJECT<VECTOR<float,2>,2>;
template class FRACTURE_OBJECT<VECTOR<float,3>,2>;
template class FRACTURE_OBJECT<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_OBJECT<VECTOR<double,2>,2>;
template class FRACTURE_OBJECT<VECTOR<double,3>,2>;
template class FRACTURE_OBJECT<VECTOR<double,3>,3>;
#endif
