#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Level_Sets/OCTREE_REMOVED_PARTICLES_PROCESSING.h>
#include <PhysBAM_Dynamics/Level_Sets/REMOVED_PARTICLES_BLENDER_3D.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;

namespace{
//#####################################################################
// Function List_Octree_Cell_To_Box
//#####################################################################
template<class T> static void
List_Octree_Cell_To_Box(OCTREE_GRID<T>& octree_grid, OCTREE_CELL<T>* cell, const ORIENTED_BOX<VECTOR<T,3> >& box,ARRAY<OCTREE_CELL<T>*>& intersecting_cells)
{
    if(!box.Intersection(cell->Bounding_Box())) return;
    if(cell->Depth_Of_This_Cell()>=octree_grid.maximum_depth){intersecting_cells.Append(cell);return;}
    assert(cell->Has_Children());
    for(int i=0;i<8;i++) List_Octree_Cell_To_Box(octree_grid,cell->Child(i),box,intersecting_cells);
}
//#####################################################################
// Function Refine_Octree_Cell_To_Box
//#####################################################################
template<class T> static void
Refine_Octree_Cell_To_Box(OCTREE_GRID<T>& octree_grid, OCTREE_CELL<T>* cell, const ORIENTED_BOX<VECTOR<T,3> >& box,ARRAY<OCTREE_CELL<T>*>& refined_cells,ARRAY<OCTREE_CELL<T>*>* intersecting_cells)
{
    if(!box.Intersection(cell->Bounding_Box())) return;
    if(cell->Depth_Of_This_Cell()>=octree_grid.maximum_depth){if(intersecting_cells)intersecting_cells->Append(cell);return;}
    if(!cell->Has_Children()){refined_cells.Append(cell);cell->Create_Children(octree_grid.number_of_cells,0,octree_grid.number_of_nodes,0,octree_grid.number_of_faces,0,&octree_grid);}
    for(int i=0;i<8;i++) Refine_Octree_Cell_To_Box(octree_grid,cell->Child(i),box,refined_cells,intersecting_cells);
}
//#####################################################################
// Function List_Octree_Cell_To_Box
//#####################################################################
template<class T> static void
List_Octree_Cell_To_Box(OCTREE_GRID<T>& octree_grid, OCTREE_CELL<T>* cell, const RANGE<VECTOR<T,3> >& box,ARRAY<OCTREE_CELL<T>*>& intersecting_cells)
{
    if(!cell->Intersection(box)) return;
    if(cell->Depth_Of_This_Cell()>=octree_grid.maximum_depth){intersecting_cells.Append(cell);return;}
    assert(cell->Has_Children());
    for(int i=0;i<8;i++) List_Octree_Cell_To_Box(octree_grid,cell->Child(i),box,intersecting_cells);
}
//#####################################################################
// Function Refine_Octree_Cell_To_Box
//#####################################################################
template<class T> static void
Refine_Octree_Cell_To_Box(OCTREE_GRID<T>& octree_grid, OCTREE_CELL<T>* cell, const RANGE<VECTOR<T,3> >& box,ARRAY<OCTREE_CELL<T>*>& refined_cells,ARRAY<OCTREE_CELL<T>*>* intersecting_cells)
{
    if(!cell->Intersection(box)) return;
    if(cell->Depth_Of_This_Cell()>=octree_grid.maximum_depth){if(intersecting_cells)intersecting_cells->Append(cell);return;}
    if(!cell->Has_Children()){refined_cells.Append(cell);cell->Create_Children(octree_grid.number_of_cells,0,octree_grid.number_of_nodes,0,octree_grid.number_of_faces,0,&octree_grid);}
    for(int i=0;i<8;i++) Refine_Octree_Cell_To_Box(octree_grid,cell->Child(i),box,refined_cells,intersecting_cells);
}
//#####################################################################
// Function Coarsen_Cell
//#####################################################################
template<class T> static bool
Coarsen_Cell(const OCTREE_CELL<T>& cell,const ARRAY<T>& phi,const T band_width_times_small_dx)
{
    T dx=(T)1.01*band_width_times_small_dx+(T)root_three*cell.DX().x;
    for(int i=0;i<8;i++) if(abs(phi(cell.Node(i))) > dx) return true;
    return false;
}
//#####################################################################
// Function Coarsen_Tree
//#####################################################################
template<class T> static bool
Coarsen_Tree(OCTREE_CELL<T>* cell,const ARRAY<T>& phi,const T refinement_distance)
{
    if(!cell->Has_Children()){return Coarsen_Cell(*cell,phi,refinement_distance);}
    bool can_coarsen=true;for(int i=0;i<8;i++) can_coarsen&=Coarsen_Tree(cell->Child(i),phi,refinement_distance);
    if(can_coarsen){cell->Delete_Children();return Coarsen_Cell(*cell,phi,refinement_distance);}
    else return false;
}
//#####################################################################
// Function Get_Cell_Bounds
//#####################################################################
template<class T> static void
Get_Cell_Bounds(OCTREE_GRID<T>& octree_grid,const RANGE<VECTOR<T,3> >& particle_bounding_box,int i,int j,int ij,int& m_start,int& m_end,int& n_start,int& n_end,int& mn_start,int& mn_end)
{
    VECTOR<T,3> particle_box_size(particle_bounding_box.Edge_Lengths());
    int cells_x=(int)ceil(particle_box_size.x/2/octree_grid.uniform_grid.dX.x);
    int cells_y=(int)ceil(particle_box_size.y/2/octree_grid.uniform_grid.dX.y);
    int cells_z=(int)ceil(particle_box_size.z/2/octree_grid.uniform_grid.dX.z);
    // THIS DOESN'T WORK (for large particle sizes) -- some problem with refining at octree boundaries
    //int m_start=max(octree_grid.cells.domain.min_corner.x,i-cells_x);int m_end=min(octree_grid.cells.domain.max_corner.x,i+cells_x);
    //int n_start=max(octree_grid.cells.domain.min_corner.y,j-cells_y);int n_end=min(octree_grid.cells.domain.max_corner.y,j+cells_y);
    //int mn_start=max(octree_grid.cells.domain.min_corner.z,ij-cells_z);int mn_end=min(octree_grid.cells.domain.max_corner.z,ij+cells_z);
    m_start=max(1,i-cells_x);m_end=min(octree_grid.uniform_grid.counts.x,i+cells_x);
    n_start=max(1,j-cells_y);n_end=min(octree_grid.uniform_grid.counts.y,j+cells_y);
    mn_start=max(1,ij-cells_z);mn_end=min(octree_grid.uniform_grid.counts.z,ij+cells_z);
}
//#####################################################################
// Function Average_Nodes
//#####################################################################
template<class T,class TV> 
void Average_Nodes(OCTREE_GRID<T>& octree_grid,GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* collision_body_list,ARRAY<TV>& node_values,const TV& default_value,
    int to_node,int from_node1,int from_node2,int from_node3,int from_node4,int from_node5,int from_node6,int from_node7,int from_node8)
{
    COLLISION_GEOMETRY_ID body_id;int triangle_id;VECTOR<T,3> intersection_point;
    int from_nodes[]={from_node1,from_node2,from_node3,from_node4,from_node5,from_node6,from_node7,from_node8};
    int count=0;TV sum=TV();
    for(int i=1;i<=8;i++){
        VECTOR<T,3> to_location=octree_grid.Node_Location(to_node);
        VECTOR<T,3> from_location=octree_grid.Node_Location(from_nodes[i-1]);
        if(!collision_body_list->collision_geometry_collection.Intersection_Between_Points(to_location,from_location,body_id,triangle_id,intersection_point)){
            sum+=node_values(from_nodes[i-1]);count++;}}
    if(count){node_values(to_node)=sum/(T)count;}
    else{node_values(to_node)=default_value;}
}
//#####################################################################
// Function Average_Nodes
//#####################################################################
template<class T,class TV> 
void Average_Nodes(OCTREE_GRID<T>& octree_grid,GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* collision_body_list,ARRAY<TV>& node_values,const TV& default_value,
    int to_node,int from_node1,int from_node2,int from_node3,int from_node4)
{
    COLLISION_GEOMETRY_ID body_id;int triangle_id;VECTOR<T,3> intersection_point;
    int from_nodes[]={from_node1,from_node2,from_node3,from_node4};
    int count=0;TV sum=TV();
    for(int i=1;i<=4;i++){
        VECTOR<T,3> to_location=octree_grid.Node_Location(to_node);
        VECTOR<T,3> from_location=octree_grid.Node_Location(from_nodes[i-1]);
        if(!collision_body_list->collision_geometry_collection.Intersection_Between_Points(to_location,from_location,body_id,triangle_id,intersection_point)){
            sum+=node_values(from_nodes[i-1]);count++;}}
    if(count){node_values(to_node)=sum/(T)count;}
    else{node_values(to_node)=default_value;}
}
//#####################################################################
// Function Average_Nodes
//#####################################################################
template<class T,class TV> 
void Average_Nodes(OCTREE_GRID<T>& octree_grid,GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* collision_body_list,ARRAY<TV>& node_values,const TV& default_value,
    int to_node,int from_node1,int from_node2)
{
    COLLISION_GEOMETRY_ID body_id;int triangle_id;VECTOR<T,3> intersection_point;
    int from_nodes[]={from_node1,from_node2};
    int count=0;TV sum=TV();
    for(int i=1;i<=2;i++){
        VECTOR<T,3> to_location=octree_grid.Node_Location(to_node);
        VECTOR<T,3> from_location=octree_grid.Node_Location(from_nodes[i-1]);
        if(!collision_body_list->collision_geometry_collection.Intersection_Between_Points(to_location,from_location,body_id,triangle_id,intersection_point)){
            sum+=node_values(from_nodes[i-1]);count++;}}
    if(count){node_values(to_node)=sum/(T)count;}
    else{node_values(to_node)=default_value;}
}
//#####################################################################
// Function Interpolate_Node_Values_To_Direct_Children
//#####################################################################
template<class T,class TV> 
void Interpolate_Node_Values_To_Direct_Children(OCTREE_GRID<T>& octree_grid,GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* collision_body_list,OCTREE_CELL<T>* cell,
    ARRAY<TV>& node_values,
    const TV& default_value)
{
    assert(cell->Has_Children());

    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(0,7),cell->Node(0),cell->Node(1),cell->Node(2),cell->Node(3),cell->Node(4),cell->Node(5),cell->Node(6),cell->Node(7)); // middle
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(0,6),cell->Node(0),cell->Node(2),cell->Node(4),cell->Node(6)); // left
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(1,7),cell->Node(1),cell->Node(3),cell->Node(5),cell->Node(7)); // right
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(0,5),cell->Node(0),cell->Node(1),cell->Node(4),cell->Node(5)); // bottom
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(2,7),cell->Node(2),cell->Node(3),cell->Node(6),cell->Node(7)); // top
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(0,3),cell->Node(0),cell->Node(1),cell->Node(2),cell->Node(3)); // front
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(4,7),cell->Node(4),cell->Node(5),cell->Node(6),cell->Node(7)); // back
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(0,4),cell->Node(0),cell->Node(4));Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(1,5),cell->Node(1),cell->Node(5));
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(2,6),cell->Node(2),cell->Node(6));Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(3,7),cell->Node(3),cell->Node(7));
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(0,1),cell->Node(0),cell->Node(1));Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(4,5),cell->Node(4),cell->Node(5));
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(2,3),cell->Node(2),cell->Node(3));Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(6,7),cell->Node(6),cell->Node(7));
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(0,2),cell->Node(0),cell->Node(2));Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(1,3),cell->Node(1),cell->Node(3));
    Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(4,6),cell->Node(4),cell->Node(6));Average_Nodes(octree_grid,collision_body_list,node_values,default_value,cell->children->Node(5,7),cell->Node(5),cell->Node(7));
}
//#####################################################################
// Function Intersection_Between_Points
//#####################################################################
template<class T>
bool Intersection_Between_Points(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,TRIANGULATED_SURFACE<T>& triangulated_surface,const ARRAY<int>& triangle_list,const T thickness_over_two)
{
    VECTOR<T,3> ray_vector=x2-x1;T ray_vector_length_squared=ray_vector.Magnitude_Squared();
    if(ray_vector_length_squared==0){
        for(int k=1;k<=triangle_list.m;k++){
            TRIANGLE_3D<T>& triangle=(*triangulated_surface.triangle_list)(triangle_list(k));
            if(triangle.Point_Inside_Triangle(x1,thickness_over_two)){return true;}}}
    else{
        T ray_t_max=sqrt(ray_vector_length_squared);RAY<VECTOR<T,3> > ray(x1,ray_vector/ray_t_max,true);ray.semi_infinite=false;ray.t_max=ray_t_max;
        for(int k=1;k<=triangle_list.m;k++){
            TRIANGLE_3D<T>& triangle=(*triangulated_surface.triangle_list)(triangle_list(k));
            if(triangle.Intersection(ray,thickness_over_two))return true;}}
    return false;}
}
//#####################################################################
// Function Set_Collision_Aware
//#####################################################################
template<class T> void OCTREE_REMOVED_PARTICLES_PROCESSING<T>::
Set_Collision_Aware(GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* collision_body_list_input)
{
    collision_body_list=collision_body_list_input;
}
//#####################################################################
// Function Get_Ellipsoid
//#####################################################################
template<class T> void OCTREE_REMOVED_PARTICLES_PROCESSING<T>::
Get_Ellipsoid(const PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles,int p,T& radius_x,T& radius_yz,TV& major_axis)
{
    T radius=scale*particles->radius(p);
    T velocity_magnitude_squared=particles->V(p).Magnitude_Squared();
    if(velocity_magnitude_squared>1e-8){ // ellipsoid
        T speed=sqrt(velocity_magnitude_squared);
        major_axis=particles->V(p)/speed;
        if(use_velocity_scaling){
            radius_x=radius+(T).5*dt*speed;
            if(preserve_volume){radius_yz=sqrt(cube(radius)/radius_x);}
            else{radius_yz=radius;}}
        else{radius_x=3*radius;radius_yz=radius;}}
    else{ // sphere
        major_axis=TV(1,0,0); // arbitrary axis
        radius_x=radius;radius_yz=radius;}
}
//#####################################################################
// Function Get_Particle_Bounding_Boxes
//#####################################################################
template<class T> void OCTREE_REMOVED_PARTICLES_PROCESSING<T>::
Get_Particle_Bounding_Boxes(ARRAY<ARRAY<ORIENTED_BOX<TV> > >& particles_oriented_bounding_box,REMOVED_PARTICLES_BLENDER_3D<T>& particle_blender,int& number_of_particles)
{
    for(int i=1;i<=particles_array.m;i++){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles = particles_array(i);
        if(particles){
            particles_oriented_bounding_box(i).Resize(particles->array_collection->Size());
            for(int p=1;p<=particles->array_collection->Size();p++) {
                T radius_x,radius_yz;TV major_axis;Get_Ellipsoid(particles,p,radius_x,radius_yz,major_axis);
                particles_oriented_bounding_box(i)(p)=particle_blender.Get_Oriented_Bounding_Box(radius_x,radius_yz,particles->X(p),major_axis);
                number_of_particles++;}}}
}
//#####################################################################
// Function Refine_And_Create_Particle_Phi
//#####################################################################
template<class T> void OCTREE_REMOVED_PARTICLES_PROCESSING<T>::
Refine_And_Create_Particle_Phi()
{
    REMOVED_PARTICLES_BLENDER_3D<T> particle_blender(blending_parameter);
    LOG::cout<<"PARAMETERS:\n\tscale:"<<scale<<"\tblending:"<<blending_parameter<<"\toctree depth:"<<octree_maximum_depth<<"\tR:"<<particle_blender.R<<std::endl;
    LOG::cout<<"\tparticle_power:"<<particle_power<<"\tuse_velocity_scaling:"<<use_velocity_scaling<<"\tdt:"<<dt<<"\tpreserve_volume:"<<preserve_volume<<std::endl;

    LOG::cout << "Setting new octree depth to " << octree_maximum_depth << std::endl;
    octree_grid.Set_Maximum_Depth(octree_maximum_depth);

    // compute particle bounding boxes
    ARRAY<ARRAY<ORIENTED_BOX<TV> > > particles_oriented_bounding_box(particles_array.m);
    ARRAY<ARRAY<RANGE<TV> > > particles_bounding_box(particles_array.m);
    int number_of_particles=0;
    Get_Particle_Bounding_Boxes(particles_oriented_bounding_box,particle_blender,number_of_particles);
    for(int i=1;i<=particles_array.m;i++)if(particles_array(i)){
        particles_bounding_box(i).Resize(particles_array(i)->array_collection->Size());
        for(int p=1;p<=particles_array(i)->array_collection->Size();p++) particles_bounding_box(i)(p)=particles_oriented_bounding_box(i)(p).Axis_Aligned_Bounding_Box();}

    LOG::cout<<"Refining octree to particles"<<std::endl;
    LOG::cout<<"Number of cells="<<octree_grid.number_of_cells<<std::endl;
    PROGRESS_INDICATOR refinement_progress(number_of_particles);
    int refinement_timer=TIMER::Singleton()->Register_Timer();
    LOG::cout <<"Number of particles:" << number_of_particles<<std::endl;

    // precompute here to be safe
    ARRAY<VECTOR<int,3> > particle_uniform_cells(particles_array.m);
    for(int t=1;t<=particles_array.m;t++) if(particles_array(t)){ 
        VECTOR<int,3> index=octree_grid.uniform_grid.Clamp_To_Cell(octree_grid.Cell_Pointer_From_Index()(t)->Center());
        particle_uniform_cells(t)=index;}

    ARRAY<OCTREE_CELL<T>*> refined_cells;refined_cells.Preallocate(50);
    for(int t=1;t<=particles_array.m;t++){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles = particles_array(t);
        if(particles){
            VECTOR<int,3> index=particle_uniform_cells(t);int i=index.x,j=index.y,ij=index.z;
            for(int p=1;p<=particles->array_collection->Size();p++){
                ORIENTED_BOX<TV>& particle_oriented_bounding_box=particles_oriented_bounding_box(t)(p);
                RANGE<TV>& particle_bounding_box=particles_bounding_box(t)(p);
                int m_start,m_end,n_start,n_end,mn_start,mn_end;
                Get_Cell_Bounds(octree_grid,particle_bounding_box,i,j,ij,m_start,m_end,n_start,n_end,mn_start,mn_end);
                for(int m=m_start;m<=m_end;m++) for(int n=n_start;n<=n_end;n++) for(int mn=mn_start;mn<=mn_end;mn++) 
                    Refine_Octree_Cell_To_Box<T>(octree_grid,octree_grid.cells(m,n,mn),particle_oriented_bounding_box,refined_cells,0);
                refinement_progress.Progress();}}}
    octree_grid.Tree_Topology_Changed();
    LOG::cout<<"time for refinement:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(refinement_timer)<<std::endl;
    LOG::cout<<"Number of cells="<<octree_grid.number_of_cells<<std::endl;

    // Fills in leaf cell values only; rest are propagated up later
    ARRAY<bool> occupied_cell;
    if(collision_body_list){
        // TODO: refresh the occupied cells in the collision body list here
        //collision_body_list->Compute_Occupied_Cells(octree_grid,occupied_cell,false,octree_grid.Minimum_Edge_Length(),2);
    }

    // MAJOR HACK -- assumes a single deformable body
    // TODO triangulated_surface isn't used anywhere...?
    TRIANGULATED_SURFACE<T>* triangulated_surface=0;
    if(collision_body_list && collision_body_list->collision_geometry_collection.bodies.m){
        triangulated_surface=&dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>&>(*collision_body_list->collision_geometry_collection.bodies(COLLISION_GEOMETRY_ID(1))).object;
        assert(triangulated_surface->hierarchy && triangulated_surface->triangle_list);}

    // TODO: collision aware refinement here
    LOG::cout << "Refining original phi" << std::endl;
    water_phi.Resize(octree_grid.number_of_nodes);
    if(collision_body_list){
        // transfer occupied cell information to coarse cells
        collision_body_list->Transfer_Occupied_Cell_Values_To_Parents();
        for(int i=1;i<=refined_cells.m;i++){
            if(collision_body_list&&collision_body_list->Occupied_Cell_Center(refined_cells(i)->Cell())) // TODO: make sure cell center is sufficient
                Interpolate_Node_Values_To_Direct_Children(octree_grid,collision_body_list,refined_cells(i),water_phi,(T)1e-5);
            else refined_cells(i)->Interpolate_Node_Values_To_Direct_Children(water_phi);}}
    else{
        for(int i=1;i<=refined_cells.m;i++){refined_cells(i)->Interpolate_Node_Values_To_Direct_Children(water_phi);}}

    LOG::cout<<"Initializing particle phi on octree"<<std::endl;
    particle_octree_phi.Resize(octree_grid.number_of_nodes,true);
    PROGRESS_INDICATOR particle_phi_progress(number_of_particles);
    int particle_field_timer=TIMER::Singleton()->Register_Timer();
    OPERATION_HASH<> node_operations(octree_grid.number_of_nodes);
    int skipped_cells=0;
    ARRAY<OCTREE_CELL<T>*> intersecting_cells;intersecting_cells.Preallocate(50);
    for(int t=1;t<=particles_array.m;t++){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles = particles_array(t);
        if(particles){
            VECTOR<int,3> index=particle_uniform_cells(t);int i=index.x,j=index.y,ij=index.z;
            for(int p=1;p<=particles->array_collection->Size();p++){
                ORIENTED_BOX<TV>& particle_oriented_bounding_box=particles_oriented_bounding_box(t)(p);
                RANGE<TV>& particle_bounding_box=particles_bounding_box(t)(p);
                T radius_x,radius_yz;TV major_axis;Get_Ellipsoid(particles,p,radius_x,radius_yz,major_axis);
                T one_over_radius_x_squared=(T)1/sqr(radius_x),one_over_radius_yz_squared=(T)1/sqr(radius_yz);

                int m_start,m_end,n_start,n_end,mn_start,mn_end;
                Get_Cell_Bounds(octree_grid,particle_bounding_box,i,j,ij,m_start,m_end,n_start,n_end,mn_start,mn_end);
                intersecting_cells.Remove_All();
                for(int m=m_start;m<=m_end;m++) for(int n=n_start;n<=n_end;n++) for(int mn=mn_start;mn<=mn_end;mn++) 
                    List_Octree_Cell_To_Box(octree_grid,octree_grid.cells(m,n,mn),particle_oriented_bounding_box,intersecting_cells);
    
                bool check_collision_aware=false;
                if(collision_body_list) for(int c=1;!check_collision_aware && c<=intersecting_cells.m;c++) 
                    if(collision_body_list&&collision_body_list->Occupied_Cell_Center(intersecting_cells(c)->Cell())) check_collision_aware=true; // TODO: make sure cell center is sufficient

                TV particle_position=particles->X(p);COLLISION_GEOMETRY_ID body_id;int triangle_id;TV intersection_point;
                for(int c=1;c<=intersecting_cells.m;c++){
                    if(check_collision_aware){
                        for(int k=0;k<8;k++){
                            int node_index=intersecting_cells(c)->Node(k); if(node_operations.Is_Marked_Current(node_index)) continue;
                            T distance=particle_blender.Get_Distance(one_over_radius_x_squared,one_over_radius_yz_squared,particles->X(p),major_axis,octree_grid.Node_Location(node_index));
                            if(distance<=particle_blender.R){
                                if(!collision_body_list->collision_geometry_collection.Intersection_Between_Points(particle_position,octree_grid.Node_Location(node_index),body_id,triangle_id,intersection_point)){
                                    particle_octree_phi(node_index)-=octree_grid.minimum_cell_size*particle_blender.C(distance)*particle_power;}}
                            node_operations.Mark(node_index);}}
                    else{
                        skipped_cells++;
                        for(int k=0;k<8;k++){
                            int node_index=intersecting_cells(c)->Node(k); if(node_operations.Is_Marked_Current(node_index)) continue;
                            T distance=particle_blender.Get_Distance(one_over_radius_x_squared,one_over_radius_yz_squared,particles->X(p),major_axis,octree_grid.Node_Location(node_index));
                            particle_octree_phi(node_index)-=octree_grid.minimum_cell_size*particle_blender.C(distance)*particle_power;
                            node_operations.Mark(node_index);}}}
                particle_phi_progress.Progress();node_operations.Next_Operation();}}}
    if(collision_body_list){LOG::cout << "Skipped " << skipped_cells << " cells using occupied cell" << std::endl;}
    intersecting_cells.Clean_Memory();
    particles_array.Delete_Pointers_And_Clean_Memory();
    LOG::cout<<"time for particle field:"<<TIMER::Singleton()->Get_Total_Time_Since_Registration(particle_field_timer)<<std::endl;
}
//#####################################################################
// Function Merge_Phi
//#####################################################################
template<class T> void OCTREE_REMOVED_PARTICLES_PROCESSING<T>::
Merge_Phi(ARRAY<T>& result)
{
    result.Resize(octree_grid.number_of_nodes,false);
    for(int i=1;i<=octree_grid.number_of_nodes;i++) result(i)=water_phi(i)+particle_octree_phi(i);
}
//#####################################################################
// Function Union_Phi
//#####################################################################
template<class T> void OCTREE_REMOVED_PARTICLES_PROCESSING<T>::
Union_Phi(ARRAY<T>& result)
{
    result.Resize(octree_grid.number_of_nodes,false);
    for(int i=1;i<=octree_grid.number_of_nodes;i++) result(i)=min(particle_octree_phi(i)+octree_grid.minimum_cell_size*blending_parameter*particle_power,water_phi(i));
}
//#####################################################################
// Function Blend_Phi
//#####################################################################
template<class T> void OCTREE_REMOVED_PARTICLES_PROCESSING<T>::
Blend_Phi(ARRAY<T>& result,const T blend_cells)
{
    result.Resize(octree_grid.number_of_nodes,false);
    for(int i=1;i<=octree_grid.number_of_nodes;i++){
        T particles_only_phi=particle_octree_phi(i)+octree_grid.minimum_cell_size*blending_parameter*particle_power;
        T alpha=clamp(water_phi(i)/(blend_cells*octree_grid.Minimum_Edge_Length()),(T)0,(T)1);
        result(i)=(1-alpha)*(water_phi(i)+particle_octree_phi(i))+alpha*particles_only_phi;}
}
//#####################################################################
template class OCTREE_REMOVED_PARTICLES_PROCESSING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OCTREE_REMOVED_PARTICLES_PROCESSING<double>;
#endif
#endif
