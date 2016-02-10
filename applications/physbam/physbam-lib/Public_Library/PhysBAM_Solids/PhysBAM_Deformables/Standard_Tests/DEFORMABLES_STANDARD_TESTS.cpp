//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_STANDARD_TESTS
//#####################################################################
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/SEGMENTED_CURVE_2D_SIGNED_DISTANCE.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/TORUS_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Standard_Tests/DEFORMABLES_STANDARD_TESTS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_STANDARD_TESTS<TV>::
DEFORMABLES_STANDARD_TESTS(EXAMPLE<TV>& example_input,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input)
    :example(example_input),deformable_body_collection(deformable_body_collection_input)
{
}
//#####################################################################
// Function Add_Gravity
//#####################################################################
template<class TV> void DEFORMABLES_STANDARD_TESTS<TV>::
Add_Gravity()
{
    // add gravity on all deformable particles
    deformable_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,true,true));
}
//#####################################################################
// Function Copy_And_Add_Structure
//#####################################################################
template<class TV> template<class T_STRUCTURE> T_STRUCTURE& DEFORMABLES_STANDARD_TESTS<TV>::
Copy_And_Add_Structure(T_STRUCTURE& structure,ARRAY<int>* particle_indices)
{
    deformable_body_collection.deformable_geometry.Add_Structure(structure.Append_Particles_And_Create_Copy(deformable_body_collection.particles,particle_indices));
    deformable_body_collection.deformable_geometry.structures.Last()->Update_Number_Nodes();
    delete &structure;
    return dynamic_cast<T_STRUCTURE&>(*deformable_body_collection.deformable_geometry.structures.Last());
}
//#####################################################################
// Function Set_Initial_Particle_Configuration
//#####################################################################
template<class TV> void DEFORMABLES_STANDARD_TESTS<TV>::
Set_Initial_Particle_Configuration(GEOMETRY_PARTICLES<TV>& particles,const RIGID_GEOMETRY_STATE<TV>& state,const bool relative_to_box_center)
{
    std::stringstream ss;ss<<"Deformable body - Total Particles : "<<particles.array_collection->Size()<<std::endl;LOG::filecout(ss.str());
    particles.Store_Velocity();
    if(relative_to_box_center){
        particles.X-=RANGE<TV>::Bounding_Box(particles.X).Center();}
    for(int p=1;p<=particles.array_collection->Size();p++){
        particles.X(p)=state.frame*(particles.X(p));
        particles.V(p)=state.Pointwise_Object_Velocity(particles.X(p));}
}
//#####################################################################
// Function Create_Segmented_Curve
//#####################################################################
template<class TV> typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Segmented_Curve(const GRID<VECTOR<T,1> >& square_grid,const RIGID_GEOMETRY_STATE<TV>& initial_state,const T density)
{
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    T_SEGMENTED_CURVE& segmented_curve=*T_SEGMENTED_CURVE::Create(particles);
    segmented_curve.Initialize_Straight_Mesh_And_Particles(square_grid);
    std::stringstream ss;ss<<"Adding Segmented Curve - Segments = "<<segmented_curve.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    Set_Mass_Of_Particles(segmented_curve,density,false);
    Set_Initial_Particle_Configuration(particles,initial_state,true);
    typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE& copy=Copy_And_Add_Structure(segmented_curve);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Segmented_Curve
//#####################################################################
template<class TV> typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Segmented_Curve(const int m,const RIGID_GEOMETRY_STATE<TV>& initial_state,const T initial_radius,const T density)
{
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    T_SEGMENTED_CURVE& segmented_curve=*T_SEGMENTED_CURVE::Create(particles);
    segmented_curve.Initialize_Circle_Mesh_And_Particles(m,initial_radius);
    std::stringstream ss;ss<<"Adding Segmented Curve - Segments = "<<segmented_curve.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    Set_Mass_Of_Particles(segmented_curve,density,false);
    Set_Initial_Particle_Configuration(particles,initial_state,true);
    typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE& copy=Copy_And_Add_Structure(segmented_curve);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Tetrahedralized_Volume
//#####################################################################
template<class TV> TETRAHEDRALIZED_VOLUME<typename TV::SCALAR>& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Tetrahedralized_Volume(const std::string& filename,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center,const bool use_constant_mass,const T density,const T scale)
{
    PHYSBAM_ASSERT(scale>0);
    PHYSBAM_ASSERT(density>0);
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    FILE_UTILITIES::Read_From_File(example.stream_type,filename,tetrahedralized_volume);
    tetrahedralized_volume.Rescale(scale);
    std::stringstream ss;ss<<"Adding Tetrahedralized Volume - Tetrahedra = "<<tetrahedralized_volume.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    Set_Mass_Of_Particles(tetrahedralized_volume,density,use_constant_mass);
    Set_Initial_Particle_Configuration(particles,initial_state,relative_to_box_center);
    TETRAHEDRALIZED_VOLUME<T>& copy=Copy_And_Add_Structure(tetrahedralized_volume);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Triangulated_Object
//#####################################################################
template<class TV> typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Triangulated_Object(const GRID<TV>& square_grid,const RIGID_GEOMETRY_STATE<TV>& initial_state,const T density)
{
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    T_TRIANGULATED_OBJECT& triangulated_object=*T_TRIANGULATED_OBJECT::Create(particles);
    triangulated_object.Initialize_Square_Mesh_And_Particles(square_grid);
    triangulated_object.Check_Signed_Areas_And_Make_Consistent(true);
    std::stringstream ss;ss<<"Adding Triangulated Object - Triangles = "<<triangulated_object.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    Set_Mass_Of_Particles(triangulated_object,density,false);
    Set_Initial_Particle_Configuration(particles,initial_state,true);
    typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT& copy=Copy_And_Add_Structure(triangulated_object);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Triangulated_Object
//#####################################################################
template<class TV> typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Triangulated_Object(const std::string& filename,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center,const bool use_constant_mass,const T scale)
{
    PHYSBAM_ASSERT(scale>0);
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    T_TRIANGULATED_OBJECT& triangulated_object=*T_TRIANGULATED_OBJECT::Create(particles);
    FILE_UTILITIES::Read_From_File(example.stream_type,filename,triangulated_object);
    triangulated_object.Rescale(scale);
    std::stringstream ss;ss<<"Adding Triangulated Object - Triangles = "<<triangulated_object.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    Set_Mass_Of_Particles(triangulated_object,density,use_constant_mass);
    Set_Initial_Particle_Configuration(particles,initial_state,relative_to_box_center);
    typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT& copy=Copy_And_Add_Structure(triangulated_object);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Segmented_Curve
//#####################################################################
template<class TV> typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Segmented_Curve(const std::string& filename,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center,const bool use_constant_mass)
{
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    T_SEGMENTED_CURVE& segmented_curve=*T_SEGMENTED_CURVE::Create(particles);
    FILE_UTILITIES::Read_From_File(example.stream_type,filename,segmented_curve);
    std::stringstream ss;ss<<"Adding Segmented Curve - Segments = "<<segmented_curve.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    Set_Mass_Of_Particles(segmented_curve,density,use_constant_mass);
    Set_Initial_Particle_Configuration(particles,initial_state,relative_to_box_center);
    typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE& copy=Copy_And_Add_Structure(segmented_curve);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Mattress
//#####################################################################
template<class TV> TRIANGULATED_AREA<typename TV::SCALAR>& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Mattress(const GRID<VECTOR<T,2> >& mattress_grid,const bool use_constant_mass,const RIGID_GEOMETRY_STATE<TV>* initial_state,const T density,const bool reverse_triangles)
{
    PHYSBAM_ASSERT(density>0);
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    TRIANGULATED_AREA<T>& triangulated_area=*TRIANGULATED_AREA<T>::Create(particles);
    triangulated_area.Initialize_Square_Mesh_And_Particles(mattress_grid,reverse_triangles);
    triangulated_area.Check_Signed_Areas_And_Make_Consistent(true);
    std::stringstream ss;ss<<"Adding Mattress - Total Triangles : "<<triangulated_area.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    particles.Store_Velocity();
    Set_Mass_Of_Particles(triangulated_area,density,use_constant_mass);
    if(initial_state) Set_Initial_Particle_Configuration(particles,*initial_state,true);
    TRIANGULATED_AREA<T>& copy=Copy_And_Add_Structure(triangulated_area);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Mattress
//#####################################################################
template<class TV> TETRAHEDRALIZED_VOLUME<typename TV::SCALAR>& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Mattress(const GRID<VECTOR<T,3> >& mattress_grid,const bool use_constant_mass,const RIGID_GEOMETRY_STATE<TV>* initial_state,const T density)
{
    PHYSBAM_ASSERT(density>0);
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
    std::stringstream ss;ss<<"Adding Mattress - Total Tetrahedra : "<<tetrahedralized_volume.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
    particles.Store_Velocity();
    Set_Mass_Of_Particles(tetrahedralized_volume,density,use_constant_mass);
    if(initial_state) Set_Initial_Particle_Configuration(particles,*initial_state,true);
    TETRAHEDRALIZED_VOLUME<T>& copy=Copy_And_Add_Structure(tetrahedralized_volume);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Create_Embedded_Tetrahedralized_Volume
//#####################################################################
template<class TV> template<class T_SHAPE> EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<typename TV::SCALAR>& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Embedded_Tetrahedralized_Volume(const T_SHAPE& shape,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center)
{
    PARTICLES<TV>& particles=*(new PARTICLES<TV>());
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& embedding=*EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::Create(particles);
    EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume=embedding.embedded_object;
    embedded_tetrahedralized_volume.Set_Interpolation_Fraction_Threshold((T)1e-2);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_tetrahedralized_volume.simplicial_object;
    GRID<TV> grid(13,13,13,RANGE<TV>::Centered_Box());tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(grid);
    ARRAY<T> phi(particles.array_collection->Size());
    for(int p=1;p<=phi.m;p++) phi(p)=shape.Signed_Distance(particles.X(p));
    // compute embedded surface
    particles.Store_Velocity();
    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    Set_Mass_Of_Particles(tetrahedralized_volume,density,true);
    embedded_tetrahedralized_volume.Calculate_Boundary_From_Levelset_On_Nodes(phi);
    tetrahedralized_volume.mesh.number_nodes=particles.array_collection->Size();
    embedding.Create_Material_Surface_From_Manifold_Embedded_Surface();
    Set_Initial_Particle_Configuration(particles,initial_state,relative_to_box_center);
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& copy=Copy_And_Add_Structure(embedding);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Substitute_Soft_Bindings_For_Nodes
//#####################################################################
template<class TV> template<class T_OBJECT> void DEFORMABLES_STANDARD_TESTS<TV>::
Substitute_Soft_Bindings_For_Nodes(T_OBJECT& object,SOFT_BINDINGS<TV>& soft_bindings,HASHTABLE<int,int>* persistent_soft_bindings,const bool embedded_only,
    const bool use_impulses_for_collisions)
{
    PARTICLES<TV>& particles=soft_bindings.particles;
    ARRAY<int> nodes;object.mesh.elements.Flattened().Get_Unique(nodes);
    ARRAY<int> map_to_new_particles(IDENTITY_ARRAY<>(particles.array_collection->Size()));
    for(int i=1;i<=nodes.m;i++) if(!embedded_only || soft_bindings.binding_list.Binding_Index_From_Particle_Index(nodes(i))){
        int node=nodes(i),bound_node;
        if(persistent_soft_bindings && persistent_soft_bindings->Get(node,bound_node))
            persistent_soft_bindings->Delete(node);
        else{
            bound_node=particles.array_collection->Add_Element_From_Deletion_List();
            particles.array_collection->Copy_Element(*particles.array_collection,node,bound_node);}
        map_to_new_particles(node)=bound_node;
        soft_bindings.Add_Binding(VECTOR<int,2>(bound_node,node),use_impulses_for_collisions);}
    for(int i=1;i<=object.mesh.elements.m;i++) object.mesh.elements(i)=map_to_new_particles.Subset(object.mesh.elements(i));
}
//#####################################################################
// Function Read_Or_Initialize_Implicit_Surface
//#####################################################################
template<class TV> LEVELSET_IMPLICIT_OBJECT<TV>* DEFORMABLES_STANDARD_TESTS<TV>::
Read_Or_Initialize_Implicit_Surface(const std::string& levelset_filename,TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface) const
{
    // undeformed levelset
    LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    if(FILE_UTILITIES::File_Exists(levelset_filename)) FILE_UTILITIES::Read_From_File(example.stream_type,levelset_filename,undeformed_levelset);
    else{
        undeformed_triangulated_surface.Update_Bounding_Box();RANGE<TV>& box=*undeformed_triangulated_surface.bounding_box;
        GRID<TV>& grid=undeformed_levelset.levelset.grid;ARRAY<T,VECTOR<int,3> >& phi=undeformed_levelset.levelset.phi;
        grid=GRID<TV>::Create_Grid_Given_Cell_Size(box,(T)1e-2*box.Edge_Lengths().Max(),false,5);
        phi.Resize(grid.Domain_Indices());
        LEVELSET_MAKER_UNIFORM<T> levelset_maker;
        levelset_maker.Verbose_Mode();
        levelset_maker.Set_Surface_Padding_For_Flood_Fill((T)1e-3);
        levelset_maker.Use_Fast_Marching_Method(true,0);
        levelset_maker.Compute_Level_Set(undeformed_triangulated_surface,grid,phi);
        FILE_UTILITIES::Create_Directory(example.output_directory);FILE_UTILITIES::Write_To_File(example.stream_type,levelset_filename,undeformed_levelset);}
    undeformed_levelset.Update_Box();
    return &undeformed_levelset;
}
//#####################################################################
// Function Initialize_Tetrahedron_Collisions
//#####################################################################
template<class TV> void DEFORMABLES_STANDARD_TESTS<TV>::
Initialize_Tetrahedron_Collisions(const int id_number,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters,
    TRIANGULATED_SURFACE<T>* triangulated_surface)
{
    triangle_collision_parameters.perform_self_collision=false;
    tetrahedralized_volume.Update_Number_Nodes();
    if(!triangulated_surface){
        if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
        triangulated_surface=tetrahedralized_volume.triangulated_surface;}
    PARTICLES<TV>& undeformed_particles=*deformable_body_collection.particles.Clone();
    TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface=*(new TRIANGULATED_SURFACE<T>(triangulated_surface->mesh,undeformed_particles));
    undeformed_triangulated_surface.Update_Triangle_List();undeformed_triangulated_surface.Initialize_Hierarchy();
    std::string levelset_filename=STRING_UTILITIES::string_sprintf("%s/common/deformable_body_undeformed_levelset_%d.phi",example.output_directory.c_str(),id_number);
    LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*Read_Or_Initialize_Implicit_Surface(levelset_filename,undeformed_triangulated_surface);
    deformable_body_collection.collisions.collision_body_list.Add_Body(new TETRAHEDRON_COLLISION_BODY<T>(tetrahedralized_volume,undeformed_triangulated_surface,undeformed_levelset,triangulated_surface),0,true);
}
//#####################################################################
// Function Create_Drifted_Surface
//#####################################################################
template<class TV> TRIANGULATED_SURFACE<typename TV::SCALAR>& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Drifted_Surface(const TRIANGULATED_SURFACE<T>& triangulated_surface,SOFT_BINDINGS<TV>& soft_bindings,const bool use_impulses_for_collisions) const
{
    assert(&triangulated_surface.particles==&soft_bindings.particles);
    PARTICLES<TV>& particles=soft_bindings.particles;
    ARRAY<int> triangulated_surface_nodes;triangulated_surface.mesh.elements.Flattened().Get_Unique(triangulated_surface_nodes);
    ARRAY<int> child_particles(particles.array_collection->Size());
    particles.array_collection->Preallocate(particles.array_collection->Size()+triangulated_surface_nodes.m);
    for(int i=1;i<=triangulated_surface_nodes.m;i++){int p=triangulated_surface_nodes(i);
        child_particles(p)=particles.array_collection->Append(*particles.array_collection,p);
        soft_bindings.Add_Binding(VECTOR<int,2>(child_particles(p),p),use_impulses_for_collisions);}
    TRIANGULATED_SURFACE<T>& drifted_surface=*TRIANGULATED_SURFACE<T>::Create(particles);
    for(int i=1;i<=triangulated_surface.mesh.elements.m;i++) drifted_surface.mesh.elements.Append(VECTOR<int,3>::Map(child_particles,triangulated_surface.mesh.elements(i)));
    return drifted_surface;
}
//#####################################################################
// Function Set_Mass_Of_Particles 
//#####################################################################
template<class TV> template<class T_OBJECT> void DEFORMABLES_STANDARD_TESTS<TV>::
Set_Mass_Of_Particles(const T_OBJECT& object,const T density,const bool use_constant_mass)
{
    PHYSBAM_ASSERT(density>0);
    if(!object.mesh.elements.m) return;
    PARTICLES<TV>& particles=dynamic_cast<PARTICLES<TV>&>(object.particles);
    particles.Store_Mass();
    ARRAY<int> nodes;object.mesh.elements.Flattened().Get_Unique(nodes);
    if(use_constant_mass&&nodes.m){
        T mass_per_node=density*object.Total_Size()/nodes.m;
        INDIRECT_ARRAY<ARRAY_VIEW<T>,ARRAY<int>&> mass_subset=particles.mass.Subset(nodes);ARRAYS_COMPUTATIONS::Fill(mass_subset,mass_per_node);}
    else{
        INDIRECT_ARRAY<ARRAY_VIEW<T>,ARRAY<int>&> mass_subset=particles.mass.Subset(nodes);ARRAYS_COMPUTATIONS::Fill(mass_subset,(T)0);
        int nodes_per_element=object.mesh.elements(1).m;
        T density_scaled=density/(T)nodes_per_element;
        for(int t=1;t<=object.mesh.elements.m;t++){
            T mass_scaled=density_scaled*object.Signed_Size(t);
            assert(mass_scaled>0);
            particles.mass.Subset(object.mesh.elements(t))+=mass_scaled;}}
}
//#####################################################################
// Function Create_Cloth_Panel
//#####################################################################
template<class TV> TRIANGULATED_SURFACE<typename TV::SCALAR>& DEFORMABLES_STANDARD_TESTS<TV>::
Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_GEOMETRY_STATE<TV>* initial_state,
    TRIANGULATED_SURFACE_CLIPPING_HELPER<T> *clipping_function,ARRAY<int>* particle_indices)
{
    PARTICLES<TV>& particles=*(new PARTICLES<TV>());
    TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create(particles);
    TRIANGLE_MESH& mesh=triangulated_surface.mesh;
    particles.Store_Mass();
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    mesh.Initialize_Herring_Bone_Mesh(m,n);particles.array_collection->Add_Elements(mesh.number_nodes);
    T mass_node=aspect_ratio*sqr(side_length)/(m*n);ARRAYS_COMPUTATIONS::Fill(particles.mass,mass_node); // TODO: make this consistent with the density attribute
    T dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) particles.X(i+m*(j-1))=TV((i-1)*dx,(T).5,(j-1)*dy);
    if(initial_state) Set_Initial_Particle_Configuration(particles,*initial_state,true);
    if(clipping_function) (*clipping_function)(triangulated_surface);
    TRIANGULATED_SURFACE<T>& copy=Copy_And_Add_Structure(triangulated_surface,particle_indices);
    delete &particles;
    return copy;
}
//#####################################################################
// Function Embed_Particles_In_Tetrahedralized_Volume
//#####################################################################
template<class TV> void DEFORMABLES_STANDARD_TESTS<TV>::
Embed_Particles_In_Tetrahedralized_Volume(BINDING_LIST<VECTOR<T,3> >& binding_list,const POINT_CLOUD_SUBSET<VECTOR<T,3>,PARTICLES<VECTOR<T,3> > >& particles_to_embed,
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const T thickness_over_two)
{
    bool tetrahedron_list_initialized=tetrahedralized_volume.tetrahedron_list!=0;if(!tetrahedron_list_initialized) tetrahedralized_volume.Update_Tetrahedron_List();
    bool hierarchy_initialized=tetrahedralized_volume.hierarchy!=0;if(!hierarchy_initialized) tetrahedralized_volume.Initialize_Hierarchy();
    const ARRAY<TETRAHEDRON<T> >& tetrahedron_list=*tetrahedralized_volume.tetrahedron_list;
    ARRAY<int> candidate_tets;
    for(int p=1;p<=particles_to_embed.Size();p++){
        VECTOR<T,3> X=particles_to_embed.X(p);T current_thickness_over_two=thickness_over_two;candidate_tets.Resize(0);
        // find tetrahedron
        int embedding_tet=0;
        while(!candidate_tets.m){
            tetrahedralized_volume.hierarchy->Intersection_List(X,candidate_tets,current_thickness_over_two);
            if(candidate_tets.m || thickness_over_two==0) break;current_thickness_over_two*=2;}
        for(int t=1;t<=candidate_tets.m;++t) if(tetrahedron_list(candidate_tets(t)).Inside(X)){embedding_tet=candidate_tets(t);break;}
        if(!embedding_tet && candidate_tets.m){ // find closest tet
            T rho_min=FLT_MAX;
            for(int t=1;t<=candidate_tets.m;++t){
                T rho=(tetrahedron_list(candidate_tets(t)).Surface(X)-X).Magnitude_Squared();if(rho<rho_min){rho_min=rho;embedding_tet=candidate_tets(t);}}}
        // add binding
        VECTOR<int,4> vertices=tetrahedralized_volume.mesh.elements(embedding_tet);
        VECTOR<T,3> weights=tetrahedron_list(embedding_tet).Barycentric_Coordinates(X);
        binding_list.Add_Binding(new LINEAR_BINDING<VECTOR<T,3>,4>(binding_list.particles,particles_to_embed.active_indices(p),vertices,weights));}
    // clean up
    if(!tetrahedron_list_initialized){delete tetrahedralized_volume.tetrahedron_list;tetrahedralized_volume.tetrahedron_list=0;}
    if(!hierarchy_initialized){delete tetrahedralized_volume.hierarchy;tetrahedralized_volume.hierarchy=0;}
}
//#####################################################################
// Function Mark_Hard_Bindings_With_Free_Particles
//#####################################################################
template<class TV> void DEFORMABLES_STANDARD_TESTS<TV>::
Mark_Hard_Bindings_With_Free_Particles()
{
    FREE_PARTICLES<TV>* free_particles=deformable_body_collection.deformable_geometry.template Find_Structure<FREE_PARTICLES<TV>*>();
    if(!free_particles){
        free_particles=FREE_PARTICLES<TV>::Create();
        deformable_body_collection.deformable_geometry.Add_Structure(free_particles);}
    for(int i=1;i<=deformable_body_collection.binding_list.bindings.m;i++)
        free_particles->nodes.Append(deformable_body_collection.binding_list.bindings(i)->particle_index);
}
//#####################################################################
#define INSTANTIATION_HELPER_1D(T) \
    template DEFORMABLES_STANDARD_TESTS<VECTOR<T,1> >::DEFORMABLES_STANDARD_TESTS(EXAMPLE<VECTOR<T,1> >&,DEFORMABLE_BODY_COLLECTION<VECTOR<T,1> >&); \
    template TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,1> >::SEGMENTED_CURVE& DEFORMABLES_STANDARD_TESTS<VECTOR<T,1> >::Create_Segmented_Curve(const GRID<VECTOR<T,1> >&,const RIGID_GEOMETRY_STATE<VECTOR<T,1> >&,const T); \
    template void DEFORMABLES_STANDARD_TESTS<VECTOR<T,1> >::Add_Gravity();

#define INSTANTIATION_HELPER_2D(T) \
    INSTANTIATION_HELPER_ALL(T,2) \
    template TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,2> >::TRIANGULATED_OBJECT& DEFORMABLES_STANDARD_TESTS<VECTOR<T,2> >::Create_Triangulated_Object(const GRID<VECTOR<T,2> >&,const RIGID_GEOMETRY_STATE<VECTOR<T,2> >&,const T);

#define INSTANTIATION_HELPER_ALL(T,d) \
    template DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::DEFORMABLES_STANDARD_TESTS(EXAMPLE<VECTOR<T,d> >&,DEFORMABLE_BODY_COLLECTION<VECTOR<T,d> >&); \
    template void DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::Set_Initial_Particle_Configuration(GEOMETRY_PARTICLES<VECTOR<T,d> >&,const RIGID_GEOMETRY_STATE<VECTOR<T,d> >&,const bool); \
    template TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,d> >::SEGMENTED_CURVE& DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::Create_Segmented_Curve(const GRID<VECTOR<T,1> >&,const RIGID_GEOMETRY_STATE<VECTOR<T,d> >&,const T); \
    template TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,d> >::SEGMENTED_CURVE& DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::Create_Segmented_Curve(const int,const RIGID_GEOMETRY_STATE<VECTOR<T,d> >&,const T,const T); \
    template TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,d> >::SEGMENTED_CURVE& DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::Create_Segmented_Curve(const std::string&,const RIGID_GEOMETRY_STATE<VECTOR<T,d> >&,const bool,const bool); \
    template TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,d> >::TRIANGULATED_OBJECT& DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::Create_Triangulated_Object(const std::string&,const RIGID_GEOMETRY_STATE<VECTOR<T,d> >&,const bool,const bool,const T); \
    template void DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::Substitute_Soft_Bindings_For_Nodes(TRIANGULATED_SURFACE<T>&,SOFT_BINDINGS<VECTOR<T,d> >&,HASHTABLE<int,int>*,const bool,const bool); \
    template void DEFORMABLES_STANDARD_TESTS<VECTOR<T,d> >::Add_Gravity();

#define INSTANTIATION_HELPER(T) \
    INSTANTIATION_HELPER_1D(T) \
    INSTANTIATION_HELPER_2D(T) \
    INSTANTIATION_HELPER_ALL(T,3) \
    template TETRAHEDRALIZED_VOLUME<VECTOR<T,3>::SCALAR>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Create_Tetrahedralized_Volume(const std::string&,const RIGID_GEOMETRY_STATE<VECTOR<T,3> >&,const bool,const bool,const T,const T); \
    template TRIANGULATED_AREA<VECTOR<T,2>::SCALAR>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,2> >::Create_Mattress(const GRID<VECTOR<T,2> >&,const bool,const RIGID_GEOMETRY_STATE<VECTOR<T,2> >*,const T,const bool); \
    template TETRAHEDRALIZED_VOLUME<VECTOR<T,3>::SCALAR>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Create_Mattress(const GRID<VECTOR<T,3> >&,const bool,const RIGID_GEOMETRY_STATE<VECTOR<T,3> >*,const T); \
    template LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Read_Or_Initialize_Implicit_Surface(const std::string&,TRIANGULATED_SURFACE<T>&) const; \
    template void DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Initialize_Tetrahedron_Collisions(const int,TETRAHEDRALIZED_VOLUME<T>&,TRIANGLE_COLLISION_PARAMETERS<VECTOR<T,3> >&,TRIANGULATED_SURFACE<T>*); \
    template TRIANGULATED_SURFACE<VECTOR<T,3>::SCALAR>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Create_Drifted_Surface(const TRIANGULATED_SURFACE<T>&,SOFT_BINDINGS<VECTOR<T,3> >&,const bool) const; \
    template EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<VECTOR<T,3>::SCALAR>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Create_Embedded_Tetrahedralized_Volume(const SPHERE<VECTOR<T,3> >&,const RIGID_GEOMETRY_STATE<VECTOR<T,3> >&,const bool); \
    template EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<VECTOR<T,3>::SCALAR>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Create_Embedded_Tetrahedralized_Volume(const TORUS<T>&,const RIGID_GEOMETRY_STATE<VECTOR<T,3> >&,const bool); \
    template TETRAHEDRALIZED_VOLUME<T>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Copy_And_Add_Structure(TETRAHEDRALIZED_VOLUME<T>&,ARRAY<int>*); \
    template TRIANGULATED_SURFACE<T>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Copy_And_Add_Structure(TRIANGULATED_SURFACE<T>&,ARRAY<int>*); \
    template TRIANGULATED_AREA<T>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,2> >::Copy_And_Add_Structure(TRIANGULATED_AREA<T>&,ARRAY<int>*); \
    template SEGMENTED_CURVE<VECTOR<T,2> >& DEFORMABLES_STANDARD_TESTS<VECTOR<T,2> >::Copy_And_Add_Structure(SEGMENTED_CURVE<VECTOR<T,2> >&,ARRAY<int>*); \
    template TRIANGULATED_SURFACE<T>& DEFORMABLES_STANDARD_TESTS<VECTOR<T,3> >::Create_Cloth_Panel(const int,const T,const T,const RIGID_GEOMETRY_STATE<VECTOR<T,3> >*, \
        TRIANGULATED_SURFACE_CLIPPING_HELPER<T>*,ARRAY<int>*); \
    template void DEFORMABLES_STANDARD_TESTS<VECTOR<T,2> >::Set_Mass_Of_Particles<SEGMENTED_CURVE<VECTOR<T,2> > >(SEGMENTED_CURVE<VECTOR<T,2> > const&,T,bool);

INSTANTIATION_HELPER(float);
template void DEFORMABLES_STANDARD_TESTS<VECTOR<float,3> >::Mark_Hard_Bindings_With_Free_Particles();
template SEGMENTED_CURVE_2D<float>& DEFORMABLES_STANDARD_TESTS<VECTOR<float,2> >::Copy_And_Add_Structure<SEGMENTED_CURVE_2D<float> >(SEGMENTED_CURVE_2D<float>&,ARRAY<int,int>*);
template void DEFORMABLES_STANDARD_TESTS<VECTOR<float,3> >::Set_Mass_Of_Particles<TRIANGULATED_SURFACE<float> >(TRIANGULATED_SURFACE<float> const&,float,bool);
template void DEFORMABLES_STANDARD_TESTS<VECTOR<float,1> >::Set_Mass_Of_Particles<SEGMENTED_CURVE<VECTOR<float,1> > >(SEGMENTED_CURVE<VECTOR<float,1> > const&,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
template void DEFORMABLES_STANDARD_TESTS<VECTOR<double,3> >::Mark_Hard_Bindings_With_Free_Particles();
template SEGMENTED_CURVE_2D<double>& DEFORMABLES_STANDARD_TESTS<VECTOR<double,2> >::Copy_And_Add_Structure<SEGMENTED_CURVE_2D<double> >(SEGMENTED_CURVE_2D<double>&,ARRAY<int,int>*);
template void DEFORMABLES_STANDARD_TESTS<VECTOR<double,3> >::Set_Mass_Of_Particles<TRIANGULATED_SURFACE<double> >(TRIANGULATED_SURFACE<double> const&,double,bool);
template void DEFORMABLES_STANDARD_TESTS<VECTOR<double,1> >::Set_Mass_Of_Particles<SEGMENTED_CURVE<VECTOR<double,1> > >(SEGMENTED_CURVE<VECTOR<double,1> > const&,double,bool);
#endif

