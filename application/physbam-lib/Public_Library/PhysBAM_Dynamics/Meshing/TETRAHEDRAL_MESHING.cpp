//#####################################################################
// Copyright 2003-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Igor Neverov, Mike Rodgers, Tamar Shinar, Eftychios Sifakis, Joey Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRAL_MESHING
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/INTERRUPTS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TETRAHEDRA.h>
#include <PhysBAM_Dynamics/Meshing/TETRAHEDRAL_MESHING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>

using namespace PhysBAM;
template<class T> TETRAHEDRAL_MESHING<T>::
TETRAHEDRAL_MESHING(const STREAM_TYPE stream_type)
    :solids_parameters(*new SOLIDS_PARAMETERS<TV>),solid_body_collection(*new SOLID_BODY_COLLECTION<TV>(new EXAMPLE_FORCES_AND_VELOCITIES<TV>(),0)),
    solids_evolution(new NEWMARK_EVOLUTION<TV>(solids_parameters,solid_body_collection)),
    implicit_surface(0),level_set_forces_and_velocities(0),
    stream_type(stream_type),output_directory("meshing_data"),frame(0),extra_refinement_criteria(0),dependent_nodes(0),boundary_mesh(0)
{
    Use_Masses_And_Springs();Set_Curvature_Subdivision_Threshold();Set_Interpolation_Error_Subdivision_Threshold();Set_Maximum_Boundary_Edge_Length();
    Set_Density();Increase_Mass_On_Boundary();Use_Dynamic_Ether_Viscosity();Use_Global_Quality_Criteria_For_Early_Exit();
    Replace_Green_Refinement_With_Embedded_T_Junctions();
    use_constant_mass=true;solids_parameters.cfl=(T).5;solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=false;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
    symmetric_initial_grid=false;
}
template<class T> TETRAHEDRAL_MESHING<T>::
~TETRAHEDRAL_MESHING()
{
    layers.Delete_Pointers_And_Clean_Memory();
    delete extra_refinement_criteria;
    delete solids_evolution;
    delete &solid_body_collection;
    delete &solids_parameters;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Initialize(IMPLICIT_OBJECT<TV>* implicit_surface_input)
{
    implicit_surface=implicit_surface_input;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
    deformable_body_collection.deformable_geometry.Add_Structure(tetrahedralized_volume);
    level_set_forces_and_velocities=new LEVEL_SET_FORCES_AND_VELOCITIES<TV>(*tetrahedralized_volume,*implicit_surface);
    solid_body_collection.example_forces_and_velocities=level_set_forces_and_velocities;
    solid_body_collection.rigid_body_collection.rigids_example_forces_and_velocities=level_set_forces_and_velocities;
    solid_body_collection.Set_CFL_Number(solids_parameters.cfl);
}
//#####################################################################
// Function Snap_Nodes_To_Level_Set_Boundary
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Snap_Nodes_To_Level_Set_Boundary(const int iterations)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    TETRAHEDRON_MESH& mesh=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    if(!mesh.boundary_nodes) mesh.Initialize_Boundary_Nodes();
    for(int t=1;t<=mesh.boundary_nodes->m;t++) for(int k=1;k<=iterations;k++){
        int node=(*mesh.boundary_nodes)(t);TV X=deformable_body_collection.particles.X(node);
        deformable_body_collection.particles.X(node)-=implicit_surface->Extended_Phi(X)*implicit_surface->Extended_Normal(X);}
    FILE_UTILITIES::Create_Directory(output_directory);
    Write_Output_Files(++frame);
}
//#####################################################################
// Function Initialize_Optimization
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Initialize_Optimization(const bool verbose)
{
    TETRAHEDRON_MESH& mesh=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    mesh.Initialize_Incident_Elements();if(!boundary_mesh) mesh.Initialize_Boundary_Mesh();else mesh.boundary_mesh=boundary_mesh;
    mesh.boundary_mesh->Initialize_Neighbor_Nodes();mesh.boundary_mesh->Initialize_Incident_Elements();
    mesh.Initialize_Boundary_Nodes(); // assumes that Initialize_Boundary_Nodes will use boundary_mesh
    map_from_nodes_to_boundary_list.Resize(mesh.number_nodes);
    for(int i=1;i<=mesh.boundary_nodes->m;i++) map_from_nodes_to_boundary_list((*mesh.boundary_nodes)(i))=i;
    for(int i=1;i<=layers.m;i++) delete layers(i);layers.Resize(1);layers(1)=mesh.boundary_nodes;
    mesh.boundary_nodes=0; // we don't need it hanging off the mesh object any more
    if(verbose) {std::stringstream ss;ss<<"boundary layer has "<<layers(1)->m<<" nodes"<<std::endl;LOG::filecout(ss.str());}
    ARRAY<bool,VECTOR<int,1> > marked(1,mesh.number_nodes);for(int i=1;i<=layers(1)->m;i++) marked((*layers(1))(i))=true;
    for(int l=2;;l++){
        layers.Append(new ARRAY<int>);
        for(int i=1;i<=layers(l-1)->m;i++){
            int j=(*layers(l-1))(i);
            for(int k=1;k<=(*mesh.incident_elements)(j).m;k++) for(int a=1;a<=4;a++){
                int b=mesh.elements((*mesh.incident_elements)(j)(k))(a);
                if(!marked(b)){layers(l)->Append(b);marked(b)=true;}}}
        if(layers(l)->m==0){delete layers(l);layers.Remove_End();break;}
        if(verbose) {std::stringstream ss;ss<<"layer "<<l<<" has "<<layers(l)->m<<" nodes"<<std::endl;LOG::filecout(ss.str());}}
    boundary_mesh_normals.Resize(layers(1)->m);
    if(replace_green_refinement_with_embedded_t_junctions)
        for(int i=1;i<=layers.m;i++) for(int j=layers(i)->m;j>=1;j--) if(!(*mesh.incident_elements)((*layers(i))(j)).m) layers(i)->Remove_Index_Lazy(j);
    Compute_Boundary_Mesh_Normals();
}
//#####################################################################
// Function Create_Final_Mesh_With_Optimization
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Create_Final_Mesh_With_Optimization(const int number_of_initial_steps,const int number_of_final_steps,const bool verbose)
{
    Write_Output_Files(frame);
    for(int i=1;i<=number_of_initial_steps;i++){
        if(verbose) {std::stringstream ss;ss<<"Working on initial iteration "<<i<<" of "<<number_of_initial_steps<<"+"<<number_of_final_steps<<std::endl;LOG::filecout(ss.str());}
        Optimization_Sweep(i/(T)(i+2),verbose);
        Write_Output_Files(++frame);}
    for(int i=1;i<=number_of_final_steps;i++){
        if(verbose) {std::stringstream ss;ss<<"Working on iteration "<<i<<" of "<<number_of_final_steps<<" (full step towards boundary)"<<std::endl;LOG::filecout(ss.str());}
        Optimization_Sweep(1,verbose);
        Write_Output_Files(++frame);}
}
//#####################################################################
// Function Optimization_Sweep
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Optimization_Sweep(const T compression_fraction,const bool verbose)
{
    worst_boundary_quality=worst_interior_quality=FLT_MAX;
    Optimize_Boundary_Layer(compression_fraction);
    std::stringstream ss;
    if(verbose) ss<<'.';
    for(int j=2;j<=layers.m;j++){Optimize_Interior_Layer(j);if(verbose) ss<<'.';}
    for(int j=layers.m;j>=2;j--){Optimize_Interior_Layer(j,true);if(verbose) ss<<'.';}
    Optimize_Boundary_Layer(compression_fraction,true);
    if(verbose) ss<<'.'<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Optimize_Boundary_Layer
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Optimize_Boundary_Layer(const T compression_fraction,const bool reverse)
{
    Check_For_Interrupts();
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    ARRAY<TV> directions(5);ARRAY<int>& nodes=*layers(1);
    for(int i=1;i<=nodes.m;i++){
        particles.X(nodes(i))-=compression_fraction*(*implicit_surface)(particles.X(nodes(i)))*boundary_mesh_normals(map_from_nodes_to_boundary_list(nodes(i)));
        Update_Dependent_Nodes(nodes(i));}
    Compute_Boundary_Mesh_Normals();
    for(int i=1;i<=nodes.m;i++){
        int p=nodes(reverse?nodes.m+1-i:i);
        TV normal=boundary_mesh_normals(map_from_nodes_to_boundary_list(p));
        if(abs(normal.x)>abs(normal.z) || abs(normal.y)>abs(normal.z)) directions(1)=TV(normal.y,-normal.x,0);
        else directions(1)=TV(normal.z,0,-normal.x);
        directions(1).Normalize();
        TV b=TV::Cross_Product(normal,directions(1));
        directions(2)=(T).30901699437494742410229341718282*directions(1)+(T).95105651629515357211643933337938*b;
        directions(3)=(T)-.80901699437494742410229341718282*directions(1)+(T).58778525229247312916870595463907*b;
        directions(4)=(T)-.80901699437494742410229341718282*directions(1)-(T).58778525229247312916870595463907*b;
        directions(5)=(T).30901699437494742410229341718282*directions(1)-(T).95105651629515357211643933337938*b;
        Search_For_Best_Position(p,directions,true);}
}
//#####################################################################
// Function Optimize_Interior_Layer
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Optimize_Interior_Layer(const int layer,const bool reverse)
{
    Check_For_Interrupts();
    ARRAY<TV> directions(7);
    directions(1)=TV(1,0,0);
    directions(2)=TV((T).21884275609895,(T).97576013861144,0);
    directions(3)=TV((T).21884282342443,-(T).59716303919821,-(T).77168913640869);
    directions(4)=TV(-(T).05406424975064,(T).23640431413810,(T).97014950247670);
    directions(5)=TV(-(T).90174918437566,-(T).34565929710323,(T).25955357597988);
    directions(6)=TV((T).21863854196520,-(T).86642677947325,(T).44888954518784);
    directions(7)=TV(-(T).58820534751966,(T).35620151622547,-(T).72604059734147);
    for(int i=1;i<=layers(layer)->m;i++){int j;if(reverse) j=(*layers(layer))(layers(layer)->m+1-i);else j=(*layers(layer))(i);Search_For_Best_Position(j,directions);}
}
//#####################################################################
// Function Search_For_Best_Position
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Search_For_Best_Position(const int node,const ARRAY<TV>& directions,bool include_boundary_terms)
{
    TETRAHEDRON_MESH& mesh=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    T best_quality=Quality_Of_Worst_Dependent_Tetrahedron(node);if(include_boundary_terms) best_quality+=Quality_Of_Worst_Dependent_Boundary_Triangle(node);
    if(use_global_quality_criteria_for_early_exit){
        if(include_boundary_terms){
            if(best_quality<worst_boundary_quality) worst_boundary_quality=best_quality;
            else if(best_quality>worst_boundary_quality+.15) return;} // early exit if good enough relative to rest of mesh
        else{
            if(best_quality<worst_interior_quality) worst_interior_quality=best_quality;
            else if(best_quality>worst_interior_quality+.1) return;}} // early exit if good enough relative to rest of mesh
    TV best_x(particles.X(node)),xi,xj,xk,xl;T alpha=FLT_MAX;PLANE<T> p;
    for(int s=1;s<=(*mesh.incident_elements)(node).m;s++){
        int t=(*mesh.incident_elements)(node)(s);
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);xi=particles.X(i);xj=particles.X(j);xk=particles.X(k);xl=particles.X(l);
        if(i==node){p.Specify_Three_Points(xj,xl,xk);alpha=min(alpha,TV::Dot_Product(xi-xj,p.normal));}
        else if(j==node){p.Specify_Three_Points(xi,xk,xl);alpha=min(alpha,TV::Dot_Product(xj-xi,p.normal));}
        else if(k==node){p.Specify_Three_Points(xi,xl,xj);alpha=min(alpha,TV::Dot_Product(xk-xi,p.normal));}
        else{p.Specify_Three_Points(xi,xj,xk);alpha=min(alpha,TV::Dot_Product(xl-xi,p.normal));}}
    alpha*=(T).05;
    int strikes=0,last_direction=1;
    while(strikes<=3){
        T localbest_quality=best_quality;TV localbest_x=best_x;
        for(int d=1;d<=directions.m;d++){
            int this_direction;if(d%2) this_direction=last_direction+d/2;else this_direction=last_direction-d/2;
            this_direction=(this_direction+directions.m-1)%directions.m+1;
            particles.X(node)=best_x+alpha*directions(this_direction);Update_Dependent_Nodes(node);
            T q=Quality_Of_Worst_Dependent_Tetrahedron(node);if(include_boundary_terms) q+=Quality_Of_Worst_Dependent_Boundary_Triangle(node);
            if(q>localbest_quality){localbest_quality=q;localbest_x=particles.X(node);last_direction=this_direction;break;}}
        if(localbest_quality>best_quality){best_quality=localbest_quality;best_x=localbest_x;}
        else{strikes++;alpha*=(T).45;}}
    particles.X(node)=best_x;Update_Dependent_Nodes(node);
}
//#####################################################################
// Function Quality_Of_Worst_Incident_Tetrahedron
//#####################################################################
template<class T> T TETRAHEDRAL_MESHING<T>::
Quality_Of_Worst_Incident_Tetrahedron(const int node)
{
    TETRAHEDRON_MESH& mesh=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    PLANE<T> p;TV xi,xj,xk,xl,n1,n2,n3,n4;T worst_quality=1;
    for(int s=1;s<=(*mesh.incident_elements)(node).m;s++){
        int t=(*mesh.incident_elements)(node)(s);
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);xi=particles.X(i);xj=particles.X(j);xk=particles.X(k);xl=particles.X(l);
        T max_edge_length=max((xi-xj).Magnitude(),(xj-xk).Magnitude(),(xk-xi).Magnitude(),(xi-xl).Magnitude(),(xj-xl).Magnitude(),(xk-xl).Magnitude());
        p.Specify_Three_Points(xj,xl,xk);n1=p.normal;T min_altitude=TV::Dot_Product(xi-xj,p.normal);
        p.Specify_Three_Points(xi,xk,xl);n2=p.normal;min_altitude=min(min_altitude,TV::Dot_Product(xj-xi,p.normal));
        p.Specify_Three_Points(xi,xl,xj);n3=p.normal;min_altitude=min(min_altitude,TV::Dot_Product(xk-xi,p.normal));
        p.Specify_Three_Points(xi,xj,xk);n4=p.normal;min_altitude=min(min_altitude,TV::Dot_Product(xl-xi,p.normal));
        T min_dihedral=min(TV::Dot_Product(n1,n2),TV::Dot_Product(n1,n3),TV::Dot_Product(n1,n4),TV::Dot_Product(n2,n3),
            TV::Dot_Product(n2,n4),TV::Dot_Product(n3,n4));
        worst_quality=min(worst_quality,min_altitude/max_edge_length+(T).1*min_dihedral);}
    return worst_quality;
}
//#####################################################################
// Function Quality_Of_Worst_Incident_Boundary_Triangle
//#####################################################################
template<class T> T TETRAHEDRAL_MESHING<T>::
Quality_Of_Worst_Incident_Boundary_Triangle(const int node)
{
    TETRAHEDRON_MESH& mesh=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    TRIANGLE_3D<T> triangle;T worst_quality=1;
    for(int s=1;s<=(*mesh.boundary_mesh->incident_elements)(node).m;s++){
        int t=(*mesh.boundary_mesh->incident_elements)(node)(s);
        int i,j,k;mesh.boundary_mesh->elements(t).Get(i,j,k);
        triangle.Specify_Three_Points(particles.X(i),particles.X(j),particles.X(k));
        worst_quality=min(worst_quality,1/triangle.Aspect_Ratio()+1/triangle.Maximum_Angle());}
    return worst_quality;
}
//#####################################################################
// Function Quality_Of_Worst_Dependent_Tetrahedron
//#####################################################################
template<class T> T TETRAHEDRAL_MESHING<T>::
Quality_Of_Worst_Dependent_Tetrahedron(const int node)
{
    T worst_quality=Quality_Of_Worst_Incident_Tetrahedron(node);
    if(replace_green_refinement_with_embedded_t_junctions) for(int i=1;i<=(*dependent_nodes)(node).m;i++){
        int dependent_node=(*dependent_nodes)(node)(i);
        worst_quality=min(worst_quality,Quality_Of_Worst_Incident_Tetrahedron(dependent_node));}
    return worst_quality;
}
//#####################################################################
// Function Quality_Of_Worst_Dependent_Boundary_Triangle
//#####################################################################
template<class T> T TETRAHEDRAL_MESHING<T>::
Quality_Of_Worst_Dependent_Boundary_Triangle(const int node)
{
    T worst_quality=Quality_Of_Worst_Incident_Boundary_Triangle(node);
    if(replace_green_refinement_with_embedded_t_junctions) for(int i=1;i<=(*dependent_nodes)(node).m;i++){
        int dependent_node=(*dependent_nodes)(node)(i);
        worst_quality=min(worst_quality,Quality_Of_Worst_Incident_Boundary_Triangle(dependent_node));}
    return worst_quality;
}
//#####################################################################
// Function Compute_Boundary_Mesh_Normals
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Compute_Boundary_Mesh_Normals()
{
    TETRAHEDRON_MESH& mesh=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    ARRAYS_COMPUTATIONS::Fill(boundary_mesh_normals,TV());PLANE<T> p;
    for(int t=1;t<=mesh.boundary_mesh->elements.m;t++){
        int i,j,k;mesh.boundary_mesh->elements(t).Get(i,j,k);
        p.Specify_Three_Points(particles.X(i),particles.X(j),particles.X(k));
        boundary_mesh_normals(map_from_nodes_to_boundary_list(i))+=p.normal;
        boundary_mesh_normals(map_from_nodes_to_boundary_list(j))+=p.normal;
        boundary_mesh_normals(map_from_nodes_to_boundary_list(k))+=p.normal;}
    for(int i=1;i<=boundary_mesh_normals.m;i++) boundary_mesh_normals(i).Normalize();
}
//#####################################################################
// Function Update_Dependent_Nodes
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Update_Dependent_Nodes(const int node)
{
    if(!replace_green_refinement_with_embedded_t_junctions) return;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    for(int i=1;i<=(*dependent_nodes)(node).m;i++)
        particles.X((*dependent_nodes)(node)(i))=binding_list.Embedded_Position((*dependent_nodes)(node)(i));
}
//#####################################################################
// Function Initialize_Dynamics
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Initialize_Dynamics()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    TETRAHEDRON_MESH& mesh=tetrahedralized_volume.mesh;
    mesh.Initialize_Adjacent_Elements();if(!boundary_mesh) mesh.Initialize_Boundary_Mesh();else mesh.boundary_mesh=boundary_mesh;
    mesh.Initialize_Boundary_Nodes();
    tetrahedralized_volume.Initialize_Triangulated_Surface();
    tetrahedralized_volume.triangulated_surface->mesh.Initialize_Incident_Elements();

    // set up dynamics
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(tetrahedralized_volume,density,use_constant_mass);
    if(boundary_mass_multiplicative_factor!=1){
        bool boundary_nodes_defined=mesh.boundary_nodes!=0;if(!boundary_nodes_defined) mesh.Initialize_Boundary_Nodes();
        for(int i=1;i<=mesh.boundary_nodes->m;i++) tetrahedralized_volume.particles.X(i)*=boundary_mass_multiplicative_factor;
        if(!boundary_nodes_defined){delete mesh.boundary_nodes;mesh.boundary_nodes=0;}}
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(deformable_body_collection.particles.mass);
    solid_body_collection.deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    if(dynamic_ether_viscosity!=0)
        solid_body_collection.Add_Force(new ETHER_DRAG<GRID<TV> >(dynamic_cast<PARTICLES<TV>&>(tetrahedralized_volume.particles),
                solid_body_collection.rigid_body_collection,true,true,dynamic_ether_viscosity));
    if(use_finite_volume) solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,
        new ROTATED_LINEAR<T,3>(youngs_modulus,poissons_ratio,Rayleigh_coefficient)));
    else if(use_masses_and_springs){
        solid_body_collection.Add_Force(Create_Edge_Springs(tetrahedralized_volume,edge_spring_stiffness,edge_spring_overdamping_fraction));
        solid_body_collection.Add_Force(Create_Altitude_Springs(tetrahedralized_volume,altitude_spring_stiffness,
            altitude_spring_overdamping_fraction,true,(T).1,true,(T).1,true,(T)0,false));}
    solid_body_collection.Update_Simulated_Particles();
    solids_evolution->Initialize_Rigid_Bodies((T)24,false);    
}
//#####################################################################
// Function Create_Final_Mesh_With_Dynamics
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Create_Final_Mesh_With_Dynamics(const T time_step,const int number_of_force_steps,const int number_of_velocity_steps,const bool verbose)
{
    Write_Output_Files(frame);

    // forces
    for(int k=1;k<=number_of_force_steps;k++){
        Advance_Dynamics((k-1)*time_step,k*time_step,verbose);
        if(verbose) {std::stringstream ss;ss<<"TIME STEP = "<<k<<", TIME = "<<" "<<k*time_step<<std::endl;LOG::filecout(ss.str());}
        Write_Output_Files(++frame);}

    // enslaved velocities
    if(verbose) {std::stringstream ss;ss<<"\n\n\n SWITCHING TO SETTING EXTERNAL VELOCITIES RATHER THAN FORCES!!!\n\n"<<std::endl;LOG::filecout(ss.str());}
    level_set_forces_and_velocities->Use_External_Velocities();
    for(int k=number_of_force_steps+1;k<=number_of_force_steps+number_of_velocity_steps;k++){
        Advance_Dynamics((k-1)*time_step,k*time_step,verbose);
        if(verbose) {std::stringstream ss;ss<<"TIME STEP = "<<k<<", TIME = "<<" "<<k*time_step<<std::endl;LOG::filecout(ss.str());}
        Write_Output_Files(++frame);}
}
//#####################################################################
// Function Advance_Dynamics
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Advance_Dynamics(const T time,const T stopping_time,const bool verbose)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    // prepare for force computation
    solid_body_collection.Update_Position_Based_State(time,true);

    T new_time=time;
    int substep=0;bool done=false;
    while(!done){substep++;
        T dt=solid_body_collection.CFL();
        if(new_time+dt>=stopping_time){dt=stopping_time-new_time;done=true;}else if(new_time+2*dt>=stopping_time) dt=(T).51*(stopping_time-new_time);
        if(verbose) {std::stringstream ss;ss<<"dt="<<dt<<"   substep="<<substep<<std::endl;LOG::filecout(ss.str());}
        solids_evolution->Advance_One_Time_Step_Position(dt,new_time,true);
        solids_evolution->Advance_One_Time_Step_Velocity(dt,new_time,true);new_time+=dt;}

    if(verbose){int index=0;
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        std::stringstream ss;
        ss<<" maxPhi="<<tetrahedralized_volume.Maximum_Magnitude_Phi_On_Boundary(*implicit_surface,&index)<<"("<<index<<")";
        index=ARRAYS_COMPUTATIONS::Arg_Maximum_Magnitude(deformable_body_collection.particles.V);
        ss<<"  maxV="<<deformable_body_collection.particles.V(index).Magnitude()<<"("<<index<<")";
        ss<<"  maxAR="<<tetrahedralized_volume.Maximum_Aspect_Ratio(&index)<<"("<<index<<")\n";
        LOG::filecout(ss.str());}
}
//#####################################################################
// Function Discard_Valence_Zero_Particles_And_Renumber
//#####################################################################
template<class TV,class T_MESH1,class T_MESH2>
void Discard_Valence_Zero_Particles_And_Renumber(PARTICLES<TV>& particles,T_MESH1& mesh1,T_MESH2& mesh2,ARRAY<int>& condensation_mapping)
{
    assert(mesh1.number_nodes==mesh2.number_nodes);

    // mark which nodes are used
    ARRAY<bool> node_is_used(mesh1.number_nodes);
    for(int t=1;t<=mesh1.elements.m;t++){
        INDIRECT_ARRAY<ARRAY<bool>,typename T_MESH1::ELEMENT_TYPE&> node_is_used_subset=node_is_used.Subset(mesh1.elements(t));ARRAYS_COMPUTATIONS::Fill(node_is_used_subset,true);}
    for(int t=1;t<=mesh2.elements.m;t++){
        INDIRECT_ARRAY<ARRAY<bool>,typename T_MESH2::ELEMENT_TYPE&> node_is_used_subset=node_is_used.Subset(mesh2.elements(t));ARRAYS_COMPUTATIONS::Fill(node_is_used_subset,true);}

    // make condensation mapping
    condensation_mapping.Resize(mesh1.number_nodes,false,false);ARRAYS_COMPUTATIONS::Fill(condensation_mapping,0);
    for(int t=1,counter=0;t<=mesh1.number_nodes;t++) if(node_is_used(t)) condensation_mapping(t)=++counter;

    // make new triangle mesh
    mesh1.number_nodes=0;
    for(int t=1;t<=mesh1.elements.m;t++){
        mesh1.elements(t)=condensation_mapping.Subset(mesh1.elements(t));
        mesh1.number_nodes=max(mesh1.number_nodes,mesh1.elements(t).Max());}
    for(int t=1;t<=mesh2.elements.m;t++){
        mesh2.elements(t)=condensation_mapping.Subset(mesh2.elements(t));
        mesh1.number_nodes=max(mesh1.number_nodes,mesh2.elements(t).Max());}
    mesh2.number_nodes=mesh1.number_nodes;

    // do particles same way
    for(int p=1;p<=condensation_mapping.m;p++) if(!condensation_mapping(p)) particles.array_collection->Add_To_Deletion_List(p);
    for(int p=condensation_mapping.m+1;p<=particles.array_collection->Size();p++) particles.array_collection->Add_To_Deletion_List(p);
    particles.array_collection->Delete_Elements_On_Deletion_List(true);particles.array_collection->Compact();

    mesh1.Refresh_Auxiliary_Structures();mesh2.Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Create_Initial_Mesh
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Create_Initial_Mesh(const T bcc_lattice_cell_size,const bool use_adaptive_refinement,const int max_subdivision_levels,const bool discard_to_get_nice_topology,const bool verbose,
    const bool use_aggressive_tet_pruning_globally,const ARRAY<RANGE<TV> >* bounding_boxes_for_aggressive_pruning)
{
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    TETRAHEDRON_MESH& mesh=tetrahedralized_volume.mesh;
    PARTICLES<TV>& particles=dynamic_cast<PARTICLES<TV>&>(tetrahedralized_volume.particles);

    // initial bcc mesh
    RANGE<TV>& box=implicit_surface->box;
    TV size=box.Edge_Lengths();T cell_size=bcc_lattice_cell_size;
    if(!bcc_lattice_cell_size){
        if(use_adaptive_refinement) cell_size=(T).1*min(size.x,size.y,size.z); // default is about 10 grid cells
        else cell_size=implicit_surface->Minimum_Cell_Size();} // use the cell size of the implicit surface
    int m=(int)ceil(size.x/cell_size),n=(int)ceil(size.y/cell_size),mn=(int)ceil(size.z/cell_size);
    GRID<TV> bcc_grid;
    if(!symmetric_initial_grid)
        bcc_grid=GRID<TV>(m+1,n+1,mn+1,box.min_corner.x,box.min_corner.x+cell_size*m,box.min_corner.y,box.min_corner.y+cell_size*n,box.min_corner.z,box.min_corner.z+cell_size*mn);
    else{
        TV center=box.Center();
        TV shift=cell_size/2*TV((T)m,(T)n,(T)mn);
        bcc_grid=GRID<TV>(m+1,n+1,mn+1,RANGE<TV>(center-shift,center+shift));}
    tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(bcc_grid);
    if(use_aggressive_tet_pruning_globally) tetrahedralized_volume.Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(*implicit_surface);
    else tetrahedralized_volume.Discard_Tetrahedrons_Outside_Implicit_Surface(*implicit_surface);
    if(bounding_boxes_for_aggressive_pruning) tetrahedralized_volume.Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(*implicit_surface,*bounding_boxes_for_aggressive_pruning);
    Check_For_Interrupts();

    // refine further if adaptive
    RED_GREEN_TETRAHEDRA<T> redgreen(tetrahedralized_volume);
    if(use_adaptive_refinement){
        ARRAY<int> tets_to_refine;tets_to_refine.Preallocate(5000);
        for(int iterations=1;iterations<=max_subdivision_levels;iterations++){
            tets_to_refine.Remove_All();
            std::stringstream ss;
            ss<<"Checking for refinement "<<std::flush;
            for(int t=1;t<=mesh.elements.m;t++){
                if(t%10000==0){
                    ss<<"."<<std::flush;
                    Check_For_Interrupts();}
                if(Tetrahedron_Refinement_Criteria(t)) tets_to_refine.Append(t);}
            ss<<std::endl;if(tets_to_refine.m==0) break;
            if(verbose) ss<<"Refining "<<tets_to_refine.m<<" out of "<<mesh.elements.m<<" tets."<<std::endl;
            redgreen.Refine_Simplex_List(tets_to_refine);
            if(verbose) ss<<"(done with iteration="<<iterations<<")"<<std::endl;LOG::filecout(ss.str());}}

    // cut off tetrahedra to get the initial mesh, optionally substitute T-junctions for green refinements and clean memory

    ARRAY<bool> keep_tet_flag;if(discard_to_get_nice_topology) Discard_To_Get_Nice_Topology(redgreen,keep_tet_flag);

    if(replace_green_refinement_with_embedded_t_junctions){
        TETRAHEDRON_MESH final_mesh;ARRAY<int> t_junctions;ARRAY<VECTOR<int,2> > t_junction_parents;
        if(discard_to_get_nice_topology){
            TETRAHEDRON_MESH minimal_mesh;
            redgreen.Coarsen_Complete_Refinements_Of_Subset(minimal_mesh,keep_tet_flag,t_junctions,t_junction_parents,allow_coarsening_to_non_graded_mesh);
            TRIANGLE_MESH minimal_boundary_mesh;
            minimal_mesh.Initialize_Boundary_Mesh_With_T_Junctions(minimal_boundary_mesh,t_junctions,t_junction_parents);
            ARRAY<bool> node_is_uncoarsenable(mesh.number_nodes);
            for(int t=1;t<=minimal_boundary_mesh.elements.m;t++){
                INDIRECT_ARRAY<ARRAY<bool>,TRIANGLE_MESH::ELEMENT_TYPE&> node_is_uncoarsenable_subset=node_is_uncoarsenable.Subset(minimal_boundary_mesh.elements(t));
                ARRAYS_COMPUTATIONS::Fill(node_is_uncoarsenable_subset,true);}
            redgreen.Coarsen_Complete_Refinements_Of_Subset(final_mesh,keep_tet_flag,t_junctions,t_junction_parents,allow_coarsening_to_non_graded_mesh,&node_is_uncoarsenable);}
        else{assert(!allow_coarsening_to_non_graded_mesh); // In the absence of discarding coarsening would just produce the unrefined BCC lattice
            redgreen.Coarsen_Green_Refinements(final_mesh,t_junctions,t_junction_parents);}
        mesh.Initialize_Mesh(final_mesh);
        boundary_mesh=new TRIANGLE_MESH;mesh.Initialize_Boundary_Mesh_With_T_Junctions(*boundary_mesh,t_junctions,t_junction_parents);
        mesh.Initialize_Incident_Elements();boundary_mesh->Initialize_Incident_Elements();
        BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
        ARRAY<int> particle_to_t_junction(particles.array_collection->Size());
        for(int i=1;i<=t_junctions.m;i++) if((*mesh.incident_elements)(t_junctions(i)).m || (*boundary_mesh->incident_elements)(t_junctions(i)).m) particle_to_t_junction(t_junctions(i))=i;
        for(int p=1;p<=particles.array_collection->Size();p++) if(particle_to_t_junction(p)){
            ARRAY<int> parents;ARRAY<T> weights;
            int t_junction=particle_to_t_junction(p);
            parents.Append(t_junction_parents(t_junction)(1));weights.Append((T).5);
            parents.Append(t_junction_parents(t_junction)(2));weights.Append((T).5);
            for(int i=1;i<=parents.m;){
                if(!particle_to_t_junction(parents(i))){i++;continue;}
                T old_weight=weights(i);t_junction=particle_to_t_junction(parents(i));
                parents.Remove_Index_Lazy(i);weights.Remove_Index_Lazy(i);
                for(int j=1;j<=2;j++){
                    int new_parent=t_junction_parents(t_junction)(j);
                    int index=parents.Find(new_parent);if(!index){index=parents.Append(new_parent);weights.Append((T)0);}
                    weights(index)+=(T).5*old_weight;}}
            switch(parents.m){
              case 2: binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,p,VECTOR<int,2>(parents(1),parents(2)),VECTOR<T,2>(weights(1),weights(2))));break;
              case 3: binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,p,VECTOR<int,3>(parents(1),parents(2),parents(3)),TV(weights(1),weights(2),weights(3))));break;
              case 4: binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,p,VECTOR<int,4>(parents(1),parents(2),parents(3),parents(4)),
                VECTOR<T,4>(weights(1),weights(2),weights(3),weights(4))));break;
              default: PHYSBAM_FATAL_ERROR();}}
        ARRAY<int> condensation_mapping;
        Discard_Valence_Zero_Particles_And_Renumber(particles,mesh,*boundary_mesh,condensation_mapping);
        for(int b=1;b<=binding_list.bindings.m;b++)
            if(LINEAR_BINDING<TV,2>* binding=dynamic_cast<LINEAR_BINDING<TV,2>*>(binding_list.bindings(b))){
                binding->particle_index=condensation_mapping(binding->particle_index);binding->parents=condensation_mapping.Subset(binding->parents);}
            else if(LINEAR_BINDING<TV,3>* binding=dynamic_cast<LINEAR_BINDING<TV,3>*>(binding_list.bindings(b))){
                binding->particle_index=condensation_mapping(binding->particle_index);binding->parents=condensation_mapping.Subset(binding->parents);}
            else if(LINEAR_BINDING<TV,4>* binding=dynamic_cast<LINEAR_BINDING<TV,4>*>(binding_list.bindings(b))){
                binding->particle_index=condensation_mapping(binding->particle_index);binding->parents=condensation_mapping.Subset(binding->parents);}
            else PHYSBAM_NOT_IMPLEMENTED();
        dependent_nodes=new ARRAY<ARRAY<int> >(mesh.number_nodes);
        for(int b=1;b<=binding_list.bindings.m;b++){
            ARRAY<int> parents=binding_list.bindings(b)->Parents();
            for(int p=1;p<=parents.m;p++) (*dependent_nodes)(parents(p)).Append(binding_list.bindings(b)->particle_index);}
        binding_list.Update_Binding_Index_From_Particle_Index();}
    else{
        for(int t=mesh.elements.m;t>=1;t--) if(!keep_tet_flag(t)) mesh.elements.Remove_Index_Lazy(t);mesh.elements.Compact();
        mesh.Delete_Auxiliary_Structures();tetrahedralized_volume.Discard_Valence_Zero_Particles_And_Renumber();}

}
//#####################################################################
// Function Tetrahedron_Refinement_Criteria
//#####################################################################
template<class T> bool TETRAHEDRAL_MESHING<T>::
Tetrahedron_Refinement_Criteria(const int index) const
{
    if(extra_refinement_criteria){
        int extra_criteria=(*extra_refinement_criteria)(index);
        if(extra_criteria) return extra_criteria>0;}

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    GEOMETRY_PARTICLES<TV>& particles=tetrahedralized_volume.particles;

    int i,j,k,l;tetrahedralized_volume.mesh.elements(index).Get(i,j,k,l);
    TV xi=particles.X(i),xj=particles.X(j),xk=particles.X(k),xl=particles.X(l);
    T max_length=sqrt(max((xi-xj).Magnitude_Squared(),(xj-xk).Magnitude_Squared(),(xk-xi).Magnitude_Squared(),(xi-xl).Magnitude_Squared(),
                                         (xj-xl).Magnitude_Squared(),(xk-xl).Magnitude_Squared()));
    T minimum_cell_size_in_tetrahedron=implicit_surface->Minimum_Cell_Size_Within_Box(TETRAHEDRON<T>(xi,xj,xk,xl).Bounding_Box());
    if(max_length<minimum_cell_size_in_tetrahedron) return false; // early exit if this cell is maximally refined
    T phi_i=implicit_surface->Extended_Phi(xi),phi_j=implicit_surface->Extended_Phi(xj),phi_k=implicit_surface->Extended_Phi(xk),phi_l=implicit_surface->Extended_Phi(xl);
    if(min(abs(phi_i),abs(phi_j),abs(phi_k),abs(phi_l)) > max_length) return false; // early exit if the surface cannot pass through the tet

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    // check sample points near interface for high curvature or interpolation error
    if(typeid(*implicit_surface)==typeid(DYADIC_IMPLICIT_OBJECT<TV>)){
        OCTREE_GRID<T>& octree_grid=dynamic_cast<DYADIC_IMPLICIT_OBJECT<TV>*>(implicit_surface)->levelset.grid;
        TETRAHEDRON<T> tetrahedron(xi,xj,xk,xl);
        ARRAY<OCTREE_CELL<T>*> intersecting_cells;
        octree_grid.Get_Cells_Intersecting_Box(tetrahedron.Bounding_Box(),intersecting_cells);
        bool seen_positive=false,seen_negative=false,seen_big_error=false;
        for(int i=1;i<=intersecting_cells.m;i++){
            TV weights,x=tetrahedron.Closest_Point(intersecting_cells(i)->Center(),weights); // sample point
            T phi=implicit_surface->Extended_Phi(x);
            if(abs(phi)<minimum_cell_size_in_tetrahedron){ // close to the interface
                VECTOR<T,2> curvatures=implicit_surface->Principal_Curvatures(x);
                if(max_length*(abs(curvatures[1])+abs(curvatures[2]))>curvature_subdivision_threshold) return true;}
            if(phi>=0) seen_positive=true;if(phi<=0) seen_negative=true;
            if(!seen_big_error){ // figure out linear interpolation of phi through the corners of the tet
                MATRIX<T,4> A(1,1,1,1,xi.x,xj.x,xk.x,xl.x,xi.y,xj.y,xk.y,xl.y,xi.z,xj.z,xk.z,xl.z);A.Invert();
                T phi0=A(1,1)*phi_i+A(1,2)*phi_j+A(1,3)*phi_k+A(1,4)*phi_l;
                TV average_normal(A(2,1)*phi_i+A(2,2)*phi_j+A(2,3)*phi_k+A(2,4)*phi_l,A(3,1)*phi_i+A(3,2)*phi_j+A(3,3)*phi_k+A(3,4)*phi_l,
                                            A(4,1)*phi_i+A(4,2)*phi_j+A(4,3)*phi_k+A(4,4)*phi_l);
                if(abs(phi-(phi0+TV::Dot_Product(average_normal,x)))>interpolation_error_subdivision_threshold*max_length) seen_big_error=true;}
            if((seen_big_error || max_length>maximum_boundary_edge_length) && seen_positive && seen_negative) return true;}}
    else
#endif
    {
        int n=(int)ceil(max_length/minimum_cell_size_in_tetrahedron);T one_over_n=(T)1/n;bool seen_positive=false,seen_negative=false,seen_big_error=false;
        for(int p=0;p<=n;p++){T a=p*one_over_n;for(int q=0;q<=n-p;q++){T b=q*one_over_n;for(int r=0;r<=n-p-q;r++){T c=r*one_over_n;
            TV x=a*xi+b*xj+c*xk+(1-a-b-c)*xl;T phi=implicit_surface->Extended_Phi(x); // sample point
            if(abs(phi)<minimum_cell_size_in_tetrahedron){ // close to the interface
                VECTOR<T,2> curvatures=implicit_surface->Principal_Curvatures(x);
                if(max_length*(abs(curvatures[1])+abs(curvatures[2]))>curvature_subdivision_threshold) return true;}
            if(phi>=0) seen_positive=true;if(phi<=0) seen_negative=true;
            if(!seen_big_error){ // figure out linear interpolation of phi through the corners of the tet
                MATRIX<T,4> A(1,1,1,1,xi.x,xj.x,xk.x,xl.x,xi.y,xj.y,xk.y,xl.y,xi.z,xj.z,xk.z,xl.z);A.Invert();
                T phi0=A(1,1)*phi_i+A(1,2)*phi_j+A(1,3)*phi_k+A(1,4)*phi_l;
                TV average_normal(A(2,1)*phi_i+A(2,2)*phi_j+A(2,3)*phi_k+A(2,4)*phi_l,A(3,1)*phi_i+A(3,2)*phi_j+A(3,3)*phi_k+A(3,4)*phi_l,
                                            A(4,1)*phi_i+A(4,2)*phi_j+A(4,3)*phi_k+A(4,4)*phi_l);
                if(abs(phi-(phi0+TV::Dot_Product(average_normal,x)))>interpolation_error_subdivision_threshold*max_length) seen_big_error=true;}
            if((seen_big_error || max_length>maximum_boundary_edge_length) && seen_positive && seen_negative) return true;}}}}

    return false;
}
//#####################################################################
// Function Discard_To_Get_Nice_Topology
//#####################################################################
// discard to guarantee no overconstrained tets, and discourage bad interior edges and non-manifold boundary nodes
template<class T> void TETRAHEDRAL_MESHING<T>::
Discard_To_Get_Nice_Topology(RED_GREEN_TETRAHEDRA<T>& redgreen,ARRAY<bool>& keep_tet_flag,const bool verbose)
{
    TETRAHEDRON_MESH& mesh=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    keep_tet_flag.Resize(mesh.elements.m,false,false);ARRAYS_COMPUTATIONS::Fill(keep_tet_flag,false);
    Envelope_Interior_Nodes(keep_tet_flag);

    TRIANGLE_MESH boundary_mesh;mesh.Initialize_Boundary_Mesh_Of_Subset(boundary_mesh,keep_tet_flag);
    boundary_mesh.Initialize_Segment_Mesh();boundary_mesh.Initialize_Neighbor_Nodes();

    ARRAY<VECTOR<int,2> > edges_to_refine;
    for(int i=1;i<=particles.array_collection->Size();i++) if((*boundary_mesh.neighbor_nodes)(i).m==3) for(int j=1;j<=3;j++) // refine degree 3 nodes
        edges_to_refine.Append(VECTOR<int,2>(i,(*boundary_mesh.neighbor_nodes)(i)(j)));
    if(verbose) {std::stringstream ss;ss<<"Subdividing "<<edges_to_refine.m<<" undesirable surface edges."<<std::endl;LOG::filecout(ss.str());}
    redgreen.Subdivide_Segment_List(edges_to_refine);edges_to_refine.Clean_Memory();

    Envelope_Interior_Nodes(keep_tet_flag);
    mesh.Initialize_Boundary_Mesh_Of_Subset(boundary_mesh,keep_tet_flag);boundary_mesh.Initialize_Segment_Mesh();
    ARRAY<bool> node_on_boundary(mesh.number_nodes);
    for(int t=1;t<=boundary_mesh.elements.m;t++){
        int i,j,k;boundary_mesh.elements(t).Get(i,j,k);node_on_boundary(i)=true;node_on_boundary(j)=true;node_on_boundary(k)=true;}
    ARRAY<int> boundary_nodes;boundary_nodes.Preallocate(boundary_mesh.elements.m);
    for(int i=1;i<=node_on_boundary.m;i++) if(node_on_boundary(i)) boundary_nodes.Append(i);
    boundary_mesh.segment_mesh->Initialize_Incident_Elements(); // for fast Segment() calls
    SEGMENT_MESH subset_segment_mesh;mesh.Initialize_Segment_Mesh_Of_Subset(subset_segment_mesh,keep_tet_flag);
    subset_segment_mesh.Initialize_Neighbor_Nodes(); // so we can look at just neighbors of boundary nodes
    ARRAY<VECTOR<int,2> > bad_segment_list;
    for(int i=1;i<=boundary_nodes.m;i++){
        int node1=boundary_nodes(i);
        for(int j=1;j<=(*subset_segment_mesh.neighbor_nodes)(node1).m;j++){
            int node2=(*subset_segment_mesh.neighbor_nodes)(node1)(j);
            if(node1<node2 && node_on_boundary(node2) && !boundary_mesh.segment_mesh->Segment(node1,node2)) bad_segment_list.Append(VECTOR<int,2>(node1,node2));}}
    int number_bad_elements=bad_segment_list.m;
    if(number_bad_elements){
        if(verbose) {std::stringstream ss;ss<<"Subdividing "<<bad_segment_list.m<<" bad interior edges."<<std::endl;LOG::filecout(ss.str());}
        redgreen.Subdivide_Segment_List(bad_segment_list);
        mesh.Initialize_Incident_Elements();Envelope_Interior_Nodes(keep_tet_flag);}

    if(number_bad_elements) mesh.Initialize_Boundary_Mesh_Of_Subset(boundary_mesh,keep_tet_flag);
    ARRAY<int> non_manifold_nodes;boundary_mesh.Non_Manifold_Nodes(non_manifold_nodes);
    if(verbose) {std::stringstream ss;ss<<"Subdividing around "<<non_manifold_nodes.m<<" bad nodes."<<std::endl;LOG::filecout(ss.str());}
    number_bad_elements+=non_manifold_nodes.m;
    if(non_manifold_nodes.m){
        ARRAY<int> tets_to_refine;tets_to_refine.Preallocate(25*non_manifold_nodes.m);
        for(int i=1;i<=non_manifold_nodes.m;i++) for(int j=1;j<=(*mesh.incident_elements)(non_manifold_nodes(i)).m;j++)
            tets_to_refine.Append((*mesh.incident_elements)(non_manifold_nodes(i))(j));
        redgreen.Refine_Simplex_List(tets_to_refine);Envelope_Interior_Nodes(keep_tet_flag);}

    while(number_bad_elements){ // continue to envelope while there may be bad things
        mesh.Initialize_Boundary_Mesh_Of_Subset(boundary_mesh,keep_tet_flag);
        ARRAYS_COMPUTATIONS::Fill(node_on_boundary,false);node_on_boundary.Resize(mesh.number_nodes);
        for(int t=1;t<=boundary_mesh.elements.m;t++){
            int i,j,k;boundary_mesh.elements(t).Get(i,j,k);node_on_boundary(i)=true;node_on_boundary(j)=true;node_on_boundary(k)=true;}
        boundary_nodes.Resize(0);int i;for(i=1;i<=node_on_boundary.m;i++) if(node_on_boundary(i)) boundary_nodes.Append(i);
        boundary_mesh.Initialize_Segment_Mesh();boundary_mesh.segment_mesh->Initialize_Incident_Elements(); // for fast Segment() calls
        mesh.Initialize_Segment_Mesh_Of_Subset(subset_segment_mesh,keep_tet_flag);
        subset_segment_mesh.Initialize_Neighbor_Nodes(); // so we can look at just neighbors of boundary nodes
        bad_segment_list.Resize(0);
        for(int i=1;i<=boundary_nodes.m;i++){
            int node1=boundary_nodes(i);
            for(int j=1;j<=(*subset_segment_mesh.neighbor_nodes)(node1).m;j++){
                int node2=(*subset_segment_mesh.neighbor_nodes)(node1)(j);
                if(node1<node2 && node_on_boundary(node2) && !boundary_mesh.segment_mesh->Segment(node1,node2)) bad_segment_list.Append(VECTOR<int,2>(node1,node2));}}
        if(verbose) {std::stringstream ss;ss<<"Enveloping "<<bad_segment_list.m<<" bad interior edges."<<std::endl;LOG::filecout(ss.str());}
        number_bad_elements=bad_segment_list.m;
        for(int i=1;i<=bad_segment_list.m;i++){
            int node1,node2;bad_segment_list(i).Get(node1,node2);T maxphi1=-(T)FLT_MAX,maxphi2=-(T)FLT_MAX;
            for(int j=1;j<=(*mesh.incident_elements)(node1).m;j++) if(!keep_tet_flag((*mesh.incident_elements)(node1)(j))){
                int a,b,c,d;mesh.elements((*mesh.incident_elements)(node1)(j)).Get(a,b,c,d);
                maxphi1=max(maxphi1,implicit_surface->Extended_Phi((T).25*(particles.X(a)+particles.X(b)+particles.X(c)+particles.X(d))));}
            for(int j=1;j<=(*mesh.incident_elements)(node2).m;j++) if(!keep_tet_flag((*mesh.incident_elements)(node2)(j))){
                int a,b,c,d;mesh.elements((*mesh.incident_elements)(node2)(j)).Get(a,b,c,d);
                maxphi2=max(maxphi2,implicit_surface->Extended_Phi((T).25*(particles.X(a)+particles.X(b)+particles.X(c)+particles.X(d))));}
            if(maxphi1<maxphi2) for(int j=1;j<=(*mesh.incident_elements)(node1).m;j++) keep_tet_flag((*mesh.incident_elements)(node1)(j))=true;
            else for(int j=1;j<=(*mesh.incident_elements)(node2).m;j++) keep_tet_flag((*mesh.incident_elements)(node2)(j))=true;}
        if(number_bad_elements) mesh.Initialize_Boundary_Mesh_Of_Subset(boundary_mesh,keep_tet_flag);
        boundary_mesh.Non_Manifold_Nodes(non_manifold_nodes);
        if(verbose) {std::stringstream ss;ss<<"Fixing "<<non_manifold_nodes.m<<" bad nodes."<<std::endl;LOG::filecout(ss.str());}
        number_bad_elements+=non_manifold_nodes.m;
        for(int i=1;i<=non_manifold_nodes.m;i++) for(int j=1;j<=(*mesh.incident_elements)(non_manifold_nodes(i)).m;j++)
            keep_tet_flag((*mesh.incident_elements)(non_manifold_nodes(i))(j))=true;}
}
//##############################################################################
// Function Envelope_Interior_Nodes
//##############################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Envelope_Interior_Nodes(ARRAY<bool>& keep_tet_flag)
{
    TETRAHEDRON_MESH& mesh=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    keep_tet_flag.Resize(mesh.elements.m);ARRAYS_COMPUTATIONS::Fill(keep_tet_flag,false);
    mesh.Initialize_Incident_Elements();mesh.Initialize_Neighbor_Nodes();
    for(int i=1;i<=mesh.number_nodes;i++){
        T phi=implicit_surface->Extended_Phi(particles.X(i));
        if(phi<0){bool envelope_node=true;
            for(int j=1;j<=(*mesh.neighbor_nodes)(i).m;j++) if(implicit_surface->Extended_Phi(particles.X((*mesh.neighbor_nodes)(i)(j)))>-3*phi){
                envelope_node=false;break;}
            if(envelope_node) for(int j=1;j<=(*mesh.incident_elements)(i).m;j++) keep_tet_flag((*mesh.incident_elements)(i)(j))=true;}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void TETRAHEDRAL_MESHING<T>::
Write_Output_Files(const int frame)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    // write state
    solid_body_collection.Write(stream_type,output_directory,frame,1,solids_parameters.write_static_variables_every_frame,solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,
        solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,false);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/tetrahedralized_volume_"+f+".tet",tetrahedralized_volume);
    // write boundary mesh and bindings
    if(replace_green_refinement_with_embedded_t_junctions && frame==0){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/boundary_mesh",*boundary_mesh);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/bindings",solid_body_collection.deformable_body_collection.binding_list);}
    // write diagnostics
    {std::ostream* output(FILE_UTILITIES::Safe_Open_Output(output_directory+"/diagnostics."+f,false));
    Read_Write<TETRAHEDRALIZED_VOLUME<T>,T>::Print_Statistics(*output,tetrahedralized_volume);
    int index;
    *output<<"max_phi = "<<tetrahedralized_volume.Maximum_Magnitude_Phi_On_Boundary(*implicit_surface,&index);*output<<" ("<<index<<")"<<std::endl;
    LINEAR_SPRINGS<TV>* linear_springs=solid_body_collection.template Find_Force<LINEAR_SPRINGS<TV>*>();
    if(linear_springs){*output<<"max_edge_compression = "<<linear_springs->Maximum_Compression_Or_Expansion_Fraction(&index);*output<<" ("<<index<<")"<<std::endl;}
    delete output;}
    // write last frame
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
template class TETRAHEDRAL_MESHING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TETRAHEDRAL_MESHING<double>;
#endif
