#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_EXAMPLE
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_SIMPLICIAL_OBJECT.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/MELTING_EXAMPLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> MELTING_EXAMPLE<TV,d>::
MELTING_EXAMPLE(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<GRID<TV> >::TYPE type)
    :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,number_of_regions,type),initialized(false)
{
    maximum_velocity=0;
    solids_parameters.write_static_variables_every_frame=true;
    use_melting=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> MELTING_EXAMPLE<TV,d>::
~MELTING_EXAMPLE()
{}
//#####################################################################
// Function Refinement_Criteria_Precomputation
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Refinement_Criteria_Precomputation(const int object,const T dt)
{
    T max_velocity;
    if(maximum_velocity) max_velocity=maximum_velocity;
    else if(melting_parameters.use_constant_melting_speed) max_velocity=melting_parameters.constant_melting_speed;
    else max_velocity=melting_parameters.levelsets(object)->V.Maximum_Magnitude();
    near_interface_threshold=-melting_parameters.interface_refinement_multiplier*max_velocity*dt;
    // pull out Fe_hat array for fast lookup during refinement criteria
    int index=melting_parameters.body_index(object);
    FINITE_VOLUME<TV,d>& fvm=solids_parameters.solid_body_collection.template Find_Force<FINITE_VOLUME<TV,d>&>(index);
    refinement_Fe_hat=&fvm.Fe_hat;
}
//#####################################################################
// Function Refinement_Criteria
//#####################################################################
template<class TV,int d> bool MELTING_EXAMPLE<TV,d>::
Refinement_Criteria(const int object,const T_RED_SIMPLEX* simplex)
{
    if(simplex->Depth()>=melting_parameters.maximum_depth) return false;
    if(melting_parameters.refine_for_high_deformation && melting_parameters.body_type(object)==melting_parameters.DEFORMABLE && Deformation_Refinement_Criteria(object,simplex)) return true;
    if(melting_parameters.refine_near_interface && Interface_Refinement_Criteria(object,simplex)) return true;
    return false;
}
//#####################################################################
// Function Update_Rigid_Body_Helper
//#####################################################################
template<class T> static void Update_Rigid_Body_Helper(MELTING_EXAMPLE<VECTOR<T,3>,3>& example,const int object)
{
    typedef VECTOR<T,3> TV;
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    MELTING_PARAMETERS<TV,3>& melting_parameters=example.melting_parameters;
    int id(melting_parameters.body_index(object));
    FRAME<TV> frame;TWIST<TV> twist; // grid to world frame, and world velocities about grid origin
    if(id){
        RIGID_BODY<TV>& rigid_body=solids_parameters.solid_body_collection.rigid_body_collection.Rigid_Body(id);
        frame=rigid_body.Frame()*melting_parameters.rigid_body_grid_frames(object).Inverse();
        twist.angular=rigid_body.Twist().angular;
        twist.linear=rigid_body.Pointwise_Object_Velocity(frame.t);
        solids_parameters.solid_body_collection.rigid_body_collection.rigid_body_particle.Remove_Body(id);}
    EMBEDDED_OBJECT<TV,3>& embedded_object=melting_parameters.levelsets(object)->embedding.embedded_object;
    embedded_object.Update_Embedded_Particle_Positions();
    TRIANGULATED_SURFACE<T>& embedded_surface=embedded_object.embedded_object;
    if(!embedded_surface.particles.array_collection->Size()) melting_parameters.body_index(object)=0;
    else{
        T density=melting_parameters.rigid_body_coupling_density;
        TRIANGULATED_SURFACE<T>* triangulated_surface=TRIANGULATED_SURFACE<T>::Create();
        triangulated_surface->mesh.Initialize_Mesh(embedded_surface.mesh);triangulated_surface->particles.array_collection->Initialize(*embedded_surface.particles.array_collection);
        RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(solids_parameters.solid_body_collection.rigid_body_collection,true);
        MASS_PROPERTIES<TV> mass_properties(*triangulated_surface,true);mass_properties.Set_Mass(rigid_body->Mass().mass);

        if(mass_properties.Volume()<1e-6){ // to prevent the rigid body simulation from going crazy
            LOG::cout << "BODY " << object << ": Volume below 1e-6, deleting body" << std::endl;
            delete rigid_body;delete triangulated_surface;melting_parameters.body_index(object)=0;}
        else{
            // generate the level set for the rigid bodies so that collisions can be performed
            FRAME<TV> frame_local;
            mass_properties.Transform_To_Object_Frame(frame_local,rigid_body->Mass().inertia_tensor,dynamic_cast<PARTICLES<TV>&>(triangulated_surface->particles));
            rigid_body->Set_Frame(frame_local);
            triangulated_surface->Update_Bounding_Box();triangulated_surface->Refresh_Auxiliary_Structures();
            LEVELSET_IMPLICIT_OBJECT<TV>* implicit_surface=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
            GRID<TV>& grid=implicit_surface->levelset.grid;ARRAY<T,VECTOR<int,3> >& phi=implicit_surface->levelset.phi;
            grid=melting_parameters.levelsets(object)->grid.uniform_grid;int multiple=(int)ceil((T)20/grid.counts.x);
            RANGE<TV> bounding_box=*triangulated_surface->bounding_box;
            bounding_box.Scale_About_Center((T)1.5);
            grid.Initialize(multiple*grid.Domain_Indices().Maximum_Corner(),bounding_box);
            LOG::cout<<grid<<std::endl;
            phi.Resize(grid.Domain_Indices());
            LEVELSET_MAKER_UNIFORM<T> levelset_maker;
            PHYSBAM_NOT_IMPLEMENTED(); //levelset_maker.Compute_Level_Set(*triangulated_surface,grid,phi);
            FRAME<TV> grid_frame=rigid_body->Frame();
            for(UNIFORM_GRID_ITERATOR_NODE<TV> iterator(grid);iterator.Valid();iterator.Next())
                phi(iterator.Node_Index())=melting_parameters.levelsets(object)->levelset.Phi(grid_frame*iterator.Location());
            rigid_body->Add_Structure(*implicit_surface);
            rigid_body->Set_Mass(max((T)0.005,mass_properties.Volume()*density));
            rigid_body->Add_Structure(*triangulated_surface);
            melting_parameters.rigid_body_grid_frames(object)=grid_frame;
            rigid_body->Set_Frame(frame*grid_frame);
            rigid_body->Twist().angular=twist.angular;rigid_body->Update_Angular_Momentum();
            rigid_body->Twist().linear=twist.linear-TV::Cross_Product(twist.angular,frame.t-rigid_body->X());
            melting_parameters.body_index(object)=Value(solids_parameters.solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body));
            if(!id) example.Initialize_Particle_Positions_And_Velocities(object);}}
}
template<class TV,int d> static void Update_Rigid_Body_Helper(MELTING_EXAMPLE<TV,d>& example,const int object)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Initialize_Bodies()
{
    if(!initialized){Initialize_Deformable_And_Rigid_Bodies();initialized=true;}

    for(int i=1;i<=melting_parameters.levelsets.m;i++)if(melting_parameters.body_type(i)==melting_parameters.RIGID)
        Update_Rigid_Body_Helper(*this,i);
    solids_parameters.solid_body_collection.rigid_body_collection.rigid_geometry_collection.Destroy_Unreferenced_Geometry();

    ARRAY<TV> X_save(solids_parameters.solid_body_collection.deformable_body_collection.particles.X);
    for(int i=1;i<=melting_parameters.levelsets.m;i++)if(melting_parameters.body_type(i)==melting_parameters.DEFORMABLE){
        LEVELSET_SIMPLICIAL_OBJECT<TV,d>& levelset=*melting_parameters.levelsets(i);
        //X_save=levelset.particles.X;
        ARRAY<TV> Xm=levelset.Get_Material_Coordinates();
        Xm+=(T)100*TV::Axis_Vector(1); // shift the objects to some other place so that the triangle collisions won't freak out (TODO: clean up)
    }//levelset.particles.X=Xm;}

    Initialize_Forces();
    for(int i=1;i<=melting_parameters.levelsets.m;i++)if(melting_parameters.body_type(i)==melting_parameters.DEFORMABLE)
        melting_parameters.levelsets(i)->embedding.Create_Material_Surface();
    solids_parameters.Initialize_Triangle_Collisions();
    solids_parameters.solid_body_collection.deformable_body_collection.particles.X=X_save;

    solids_parameters.solid_body_collection.deformable_body_collection.collisions.Initialize_Object_Collisions(solids_parameters.collide_with_interior,solids_parameters.collision_tolerance);
    solids_parameters.solid_body_collection.Update_Time_Varying_Material_Properties(0); // inaccurate time
    solids_parameters.solid_body_collection.Update_Position_Based_State(0); // inaccurate time

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Add_Deformable_Melting_Object
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Add_Deformable_Melting_Object(const int index)
{
    int object=Add_Melting_Object(melting_parameters.DEFORMABLE,index,solids_parameters.solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,d>&>());
    Initialize_Particle_Positions_And_Velocities(object);
}
//#####################################################################
// Function Add_Rigid_Melting_Object
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Add_Rigid_Melting_Object(const int index)
{
    Add_Melting_Object(melting_parameters.RIGID,index,*EMBEDDED_MATERIAL_SURFACE<TV,d>::Create());
}
//#####################################################################
// Function Add_Melting_Object
//#####################################################################
template<class TV,int d> int MELTING_EXAMPLE<TV,d>::
Add_Melting_Object(const typename MELTING_PARAMETERS<TV,d>::BODY_TYPE type,const int index,EMBEDDED_MATERIAL_SURFACE<TV,d>& embedding)
{
    melting_parameters.body_type.Append(type);
    melting_parameters.body_index.Append(index);
    melting_parameters.levelsets.Append(new LEVELSET_SIMPLICIAL_OBJECT<TV,d>(embedding)); // TODO: fix memory leaks
    melting_parameters.rigid_body_grid_frames.Append(FRAME<TV>());
    melting_parameters.temperature.Append(0);melting_parameters.reaction.Append(0);

    int object=melting_parameters.levelsets.m;
    LEVELSET_SIMPLICIAL_OBJECT<TV,d>& levelset=*melting_parameters.levelsets(object);
    T_RED_GREEN_GRID& grid=levelset.grid;

    Initialize_Grid(object,grid);
    Setup_Initial_Refinement(object);

    levelset.Initialize();
    Initialize_Phi(object,levelset.phi);
    Initialize_Levelset_Velocity(object,levelset.V);

    levelset.Build_Embedded_Object(true);
    //levelset.particles.X=levelset.Get_Material_Coordinates();
    return object;
}
//#####################################################################
// Function Mark_Nodes_Inside_Objects // only in 2d
//#####################################################################
template<class T> void Mark_Nodes_Inside_Objects_Helper(MELTING_EXAMPLE<VECTOR<T,2>,2>& example,ARRAY<bool,VECTOR<int,2> >& inside_objects)
{
    GRID<VECTOR<T,2> >& grid=*example.fluids_parameters.grid;
    for(int object=1;object<=example.melting_parameters.levelsets.m;object++){
        TRIANGULATED_AREA<T>& material_surface=example.melting_parameters.levelsets(object)->embedding.material_surface;
        material_surface.Initialize_Hierarchy();material_surface.Update_Bounding_Box();
        for(int i=1;i<=grid.counts.x;i++)for(int j=1;j<=grid.counts.y;j++)inside_objects(i,j)=material_surface.Inside(grid.X(i,j))!=0;}
}
template<class T_EXAMPLE> void Mark_Nodes_Inside_Objects_Helper(T_EXAMPLE&,ARRAY<bool,VECTOR<int,2> >&)
{
    PHYSBAM_FATAL_ERROR();
}
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Mark_Nodes_Inside_Objects(ARRAY<bool,VECTOR<int,2> >& inside_objects)
{
    Mark_Nodes_Inside_Objects_Helper(*this,inside_objects);
}
//#####################################################################
// Function Melting_Levelset_Substep
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Melting_Levelset_Substep(const int object,const T dt,const T time)
{
    LEVELSET_SIMPLICIAL_OBJECT<TV,d>& levelset=*melting_parameters.levelsets(object);
    if(melting_parameters.use_constant_melting_speed)levelset.phi+=melting_parameters.constant_melting_speed*dt;
    else levelset.Euler_Step(dt,time);
}
//#####################################################################
// Function Update_Solids_Topology_For_Melting
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Update_Solids_Topology_For_Melting(const T dt,const T time,const bool reinitialize)
{
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        LEVELSET_SIMPLICIAL_OBJECT<TV,d>& levelset=*melting_parameters.levelsets(object);
        ARRAY<T> phi_save=levelset.phi;Melting_Levelset_Substep(object,dt,time);
        for(int i=1;i<=phi_save.m;i++)levelset.phi(i)=max(levelset.phi(i),phi_save(i)); // phi should always increase
        if(reinitialize) levelset.levelset.Fast_Marching_Method();}
    Update_Topology(dt);
}
//#####################################################################
// Function Update_Topology
//#####################################################################
template<class TV,int d> bool MELTING_EXAMPLE<TV,d>::
Refinement_Criteria_Helper(PAIR<MELTING_EXAMPLE<TV,d>*,int>* helper,const typename RED_GREEN_POLICY<VECTOR<T,d> >::RED_SIMPLEX* simplex)
{
    return helper->x->Refinement_Criteria(helper->y,simplex);
}
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Update_Topology(const T dt)
{
    PAIR<MELTING_EXAMPLE<TV,d>*,int> helper;helper.x=this;
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        LEVELSET_SIMPLICIAL_OBJECT<TV,d>& levelset=*melting_parameters.levelsets(object);
        T_RED_GREEN_GRID& grid=levelset.grid;
        POINT_CLOUD_SUBSET<TV,GEOMETRY_PARTICLES<TV> >& particles=levelset.particles;
        ARRAY<T>* temperature=melting_parameters.temperature(object);
        ARRAY<T>* reaction=melting_parameters.reaction(object);

        // map data to the cell based indices
        ARRAY<TV> cell_based_X(grid.number_of_nodes),cell_based_V(grid.number_of_nodes);
        for(int i=1;i<=grid.number_of_nodes;i++)if(levelset.node_to_particle_mapping(i))
            cell_based_X(i)=particles.X(levelset.node_to_particle_mapping(i));

        // refine
        int old_number_of_nodes=grid.number_of_nodes;
        Refinement_Criteria_Precomputation(object,dt);
        helper.y=object;grid.Refine_One_Level(&helper,Refinement_Criteria_Helper);
        LOG::cout<<"Refinement adding "<<grid.number_of_nodes-old_number_of_nodes<<" nodes."<<std::endl;

        // interpolate position, velocity to new particles
        cell_based_X.Resize(grid.number_of_nodes);grid.Interpolate_Node_Values_To_New_Nodes(cell_based_X,old_number_of_nodes);
        cell_based_V.Resize(grid.number_of_nodes);grid.Interpolate_Node_Values_To_New_Nodes(cell_based_V,old_number_of_nodes);
        levelset.phi.Resize(grid.number_of_nodes);grid.Interpolate_Node_Values_To_New_Nodes(levelset.phi,old_number_of_nodes);
        levelset.V.Resize(grid.number_of_nodes);grid.Interpolate_Node_Values_To_New_Nodes(levelset.V,old_number_of_nodes);
        if(temperature){temperature->Resize(grid.number_of_nodes);grid.Interpolate_Node_Values_To_New_Nodes(*temperature,old_number_of_nodes);}
        if(reaction){reaction->Resize(grid.number_of_nodes);grid.Interpolate_Node_Values_To_New_Nodes(*reaction,old_number_of_nodes);}

        // compact red_green_grid
        {ARRAY<int> node_mapping_array;grid.Compact_Array_Indices(0,&node_mapping_array);
        ARRAY<T>::Compact_Array_Using_Compaction_Array(levelset.phi,node_mapping_array);
        ARRAY<T>::Compact_Array_Using_Compaction_Array(levelset.V,node_mapping_array);
        if(temperature) ARRAY<T>::Compact_Array_Using_Compaction_Array(*temperature,node_mapping_array);
        if(reaction) ARRAY<T>::Compact_Array_Using_Compaction_Array(*reaction,node_mapping_array);
        ARRAY<TV> temp(node_mapping_array.m);
        ARRAY<TV>::Compact_Array_Using_Compaction_Array(cell_based_X,node_mapping_array,&temp);
        ARRAY<TV>::Compact_Array_Using_Compaction_Array(cell_based_V,node_mapping_array,&temp);}

        grid.Tree_Topology_Changed();
        levelset.levelset.Tree_Topology_Changed();
        levelset.Build_Embedded_Object(true);

        // map data back to the new particles
        for(int i=1;i<=grid.number_of_nodes;i++) if(levelset.node_to_particle_mapping(i))
            particles.X(levelset.node_to_particle_mapping(i))=cell_based_X(i);}

    Initialize_Bodies();
}
//#####################################################################
// Function Interface_Refinement_Criteria
//#####################################################################
template<class TV,int d> bool MELTING_EXAMPLE<TV,d>::
Interface_Refinement_Criteria(const int object,const T_RED_SIMPLEX* simplex)
{
    for(int i=0;i<d+1;i++)
        if(melting_parameters.levelsets(object)->phi(simplex->Node(i))<=0 && melting_parameters.levelsets(object)->phi(simplex->Node(i))>=near_interface_threshold)return true;
    return false;
}
//#####################################################################
// Function Deformation_Refinement_Criteria
//#####################################################################
template<class TV,int d> bool MELTING_EXAMPLE<TV,d>::
Deformation_Refinement_Criteria(const int object,const T_RED_SIMPLEX* simplex)
{
    T inverse_scale=pow(melting_parameters.successive_refinement_multiplier,simplex->Depth());
    T expansion=melting_parameters.expansive_refinement_threshold*inverse_scale+1,compression=1-melting_parameters.compressive_refinement_threshold*inverse_scale;
    int t=melting_parameters.levelsets(object)->cell_to_simplex_mapping(simplex->cell);
    if(t) return (*refinement_Fe_hat)(t).x11>expansion||(*refinement_Fe_hat)(t).Last()<compression;
    else if(simplex->Has_Green_Children())
        for(int cell=1;cell<=simplex->green_children->cells.m;cell++){
            int t=melting_parameters.levelsets(object)->cell_to_simplex_mapping(simplex->green_children->cells(cell));if(!t)continue;
            if((*refinement_Fe_hat)(t).x11>expansion||(*refinement_Fe_hat)(t).Last()<compression) return true;}
    return false;
}
//#####################################################################
// Function Initial_Refinement_Criteria_Helper
//#####################################################################
template<class TV,int d> bool MELTING_EXAMPLE<TV,d>::
Initial_Refinement_Criteria_Helper(INITIAL_REFINEMENT_CRITERIA_HELPER* helper,const typename RED_GREEN_POLICY<VECTOR<T,d> >::RED_SIMPLEX* simplex)
{
    if(simplex->Depth()>=*helper->maximum_depth)return false;
    bool positive=false,negative=false;
    for(int i=0;i<d+1;i++){
        if((*helper->phi)(simplex->Node(i))>0)positive=true;
        else negative=true;}
    return positive&&negative;
}
//#####################################################################
// Function Setup_Initial_Refinement
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Setup_Initial_Refinement(const int object)
{
    melting_parameters.levelsets(object)->Initialize();
    T_RED_GREEN_GRID& grid=melting_parameters.levelsets(object)->grid;
    ARRAY<T>& phi=melting_parameters.levelsets(object)->phi;
    INITIAL_REFINEMENT_CRITERIA_HELPER helper;helper.phi=&phi;helper.maximum_depth=&melting_parameters.maximum_depth;
    Initialize_Phi(object,phi);
    if(!melting_parameters.refine_near_interface)return;
    // refine the tree
    LOG::cout<<"Starting with "<<grid.number_of_cells<<" cells, "<<grid.number_of_nodes<<" nodes"<<std::endl;
    for(;;){
        int old_number_of_cells=grid.number_of_cells;
        Refinement_Criteria_Precomputation(object,0);
        grid.Refine_One_Level(&helper,Initial_Refinement_Criteria_Helper);
        if(grid.number_of_cells==old_number_of_cells) break;
        LOG::cout<<"Refined to "<<grid.number_of_cells<<" cells, "<<grid.number_of_nodes<<" nodes"<<std::endl;
        grid.Tree_Topology_Changed();
        phi.Resize(grid.number_of_nodes);
        Initialize_Phi(object,phi);}
    LOG::cout<<"Initial refinement complete"<<std::endl;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Read_Output_Files_Solids(const int frame)
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Read_Output_Files_Solids(frame);
    std::string prefix=STRING_UTILITIES::string_sprintf("%s/%d/",output_directory.c_str(),frame);
    ARRAY<int> body_type;FILE_UTILITIES::Read_From_File(stream_type,prefix+"melting_body_type",body_type);
    melting_parameters.body_type.Resize(body_type.m);for(int i=1;i<=body_type.m;i++)melting_parameters.body_type(i)=(typename MELTING_PARAMETERS<TV,d>::BODY_TYPE)body_type(i);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"melting_body_index",melting_parameters.body_index);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"melting_rigid_body_grid_frames",melting_parameters.rigid_body_grid_frames);
    melting_parameters.temperature.Delete_Pointers_And_Clean_Memory();melting_parameters.reaction.Delete_Pointers_And_Clean_Memory();
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"melting_temperature",melting_parameters.temperature);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"melting_reaction",melting_parameters.reaction);
    for(int i=1;i<=melting_parameters.levelsets.m;i++){
        std::string o=STRING_UTILITIES::string_sprintf("%smelting_%d_",prefix.c_str(),i);
        delete melting_parameters.levelsets(i);
        if(melting_parameters.body_type(i)==melting_parameters.DEFORMABLE)
            melting_parameters.levelsets(i)=new LEVELSET_SIMPLICIAL_OBJECT<TV,d>(solids_parameters.solid_body_collection.deformable_body_collection.deformable_geometry
                .template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,d>&>(melting_parameters.body_index(i)));
        else{
            melting_parameters.levelsets(i)=LEVELSET_SIMPLICIAL_OBJECT<TV,d>::Create();
            FILE_UTILITIES::Read_From_File(stream_type,o+"levelset_embedded_object",melting_parameters.levelsets(i)->embedding.embedded_object);
            melting_parameters.levelsets(i)->embedding.particles.Store_Velocity(true);} // behold the ugliness of geometry i/o
        FILE_UTILITIES::Read_From_File(stream_type,o+"levelset",*melting_parameters.levelsets(i));
        FILE_UTILITIES::Read_From_File(stream_type,o+"velocity",melting_parameters.levelsets(i)->V);}
    Initialize_Bodies();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV,int d> void MELTING_EXAMPLE<TV,d>::
Write_Output_Files(const int frame) const
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Write_Output_Files(frame);
    std::string prefix=STRING_UTILITIES::string_sprintf("%s/%d/",output_directory.c_str(),frame);
    ARRAY<int> body_type(melting_parameters.body_type.m);for(int i=1;i<=body_type.m;i++)body_type(i)=melting_parameters.body_type(i);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"melting_body_type",body_type);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"melting_body_index",melting_parameters.body_index);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"melting_rigid_body_grid_frames",melting_parameters.rigid_body_grid_frames);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"melting_temperature",melting_parameters.temperature);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"melting_reaction",melting_parameters.reaction);
    if(fluids_parameters.write_debug_data || (fluids_parameters.write_restart_data && frame%fluids_parameters.restart_data_write_rate==0))
        for(int i=1;i<=melting_parameters.levelsets.m;i++){
            std::string o=STRING_UTILITIES::string_sprintf("%smelting_%d_",prefix.c_str(),i);
            FILE_UTILITIES::Write_To_File(stream_type,o+"levelset",*melting_parameters.levelsets(i));
            FILE_UTILITIES::Write_To_File(stream_type,o+"velocity",melting_parameters.levelsets(i)->V);
            if(melting_parameters.body_type(i)==melting_parameters.RIGID)
                FILE_UTILITIES::Write_To_File(stream_type,o+"levelset_embedded_object",melting_parameters.levelsets(i)->embedding.embedded_object);
            melting_parameters.levelsets(i)->levelset.Lazy_Update_Overlay_Levelset(&melting_parameters.levelsets(i)->V);
            FILE_UTILITIES::Write_To_File(stream_type,o+"dyadic_grid",*melting_parameters.levelsets(i)->levelset.overlay_grid);
            FILE_UTILITIES::Write_To_File(stream_type,o+"dyadic_levelset",*melting_parameters.levelsets(i)->levelset.overlay_phi);
            FILE_UTILITIES::Write_To_File(stream_type,o+"dyadic_velocities",melting_parameters.levelsets(i)->levelset.overlay_velocity);}
}
//#####################################################################
template class MELTING_EXAMPLE<VECTOR<float,2>,2>;
template class MELTING_EXAMPLE<VECTOR<float,3>,2>;
template class MELTING_EXAMPLE<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MELTING_EXAMPLE<VECTOR<double,2>,2>;
template class MELTING_EXAMPLE<VECTOR<double,3>,2>;
template class MELTING_EXAMPLE<VECTOR<double,3>,3>;
#endif
#endif
