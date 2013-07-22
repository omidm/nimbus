//#####################################################################
// Copyright 2008, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_REGION
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_REGION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> FRACTURE_REGION<T>::
FRACTURE_REGION(TRIANGULATED_SURFACE<T>* triangulated_surface_input,LEVELSET_IMPLICIT_OBJECT<TV>* implicit_object_input,const bool create_particle_partition)
    :triangulated_surface(triangulated_surface_input),implicit_object(implicit_object_input),particle_intersection_thickness((T).05),need_destroy_data(true),particle_partition(0)
{
    if(create_particle_partition) Initialize_Particle_Partition();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FRACTURE_REGION<T>::
~FRACTURE_REGION()
{
    if(need_destroy_data){delete triangulated_surface;delete implicit_object;}
    if(particle_partition) delete particle_partition;
}
//#####################################################################
// Function Intersect_With_Rigid_Body
//#####################################################################
template<class T> ARRAY<FRACTURE_REGION<T>*> FRACTURE_REGION<T>::
Intersect_With_Rigid_Body(const FRACTURE_REGION<T>& body,const bool use_particle_optimization,const bool tessellate_region)
{
    ARRAY<FRACTURE_REGION<T>*> new_regions;
    TV_INT offset=body.fracture_offset-fracture_offset;
    T phi_scale=body.implicit_object->levelset.grid.dX.x/implicit_object->levelset.grid.dX.x;

    // Initialize the new region's levelset.  The levelset locations are in the old region's object space.
    LEVELSET_IMPLICIT_OBJECT<TV>* fragment_implicit_object=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    RANGE<TV_INT> intersection_domain=RANGE<TV_INT>::Intersect(body.implicit_object->levelset.grid.Domain_Indices()-offset,implicit_object->levelset.grid.Domain_Indices());
    RANGE<TV_INT> fragment_counts=RANGE<TV_INT>::Zero_Box().Thickened(-INT_MAX); // Empty_Box() is not implemented properly for integers.
    for(NODE_ITERATOR iterator(implicit_object->levelset.grid,intersection_domain);iterator.Valid();iterator.Next()){
        if(implicit_object->levelset.phi(iterator.index)>0) continue;
        if(body.implicit_object->levelset.phi.Valid_Index(iterator.index+offset) && body.implicit_object->levelset.phi(iterator.index+offset)>0) continue;
        fragment_counts.Enlarge_To_Include_Point(iterator.index);}
    if(fragment_counts.Empty()){delete fragment_implicit_object;return new_regions;}
    fragment_counts=RANGE<TV_INT>::Intersect(intersection_domain,fragment_counts.Thickened(2));
    TV min_corner=implicit_object->levelset.grid.Node(fragment_counts.min_corner);
    TV max_corner=implicit_object->levelset.grid.Node(fragment_counts.max_corner);
    RANGE<TV> fragment_domain=RANGE<TV>(min_corner,max_corner);
    fragment_implicit_object->levelset.grid.Initialize(fragment_counts.Edge_Lengths()+TV_INT::All_Ones_Vector(),fragment_domain,implicit_object->levelset.grid.Is_MAC_Grid());
    fragment_implicit_object->levelset.phi.Resize(fragment_implicit_object->levelset.grid.Domain_Indices());
    fragment_implicit_object->Update_Box();fragment_implicit_object->Update_Minimum_Cell_Size();
    ARRAY<int,VECTOR<int,3> > colors(fragment_implicit_object->levelset.grid.Domain_Indices());
    for(NODE_ITERATOR iterator(fragment_implicit_object->levelset.grid);iterator.Valid();iterator.Next()){
        TV_INT my_cell_index=iterator.index+fragment_counts.min_corner-TV_INT::All_Ones_Vector();
        T body_phi=body.implicit_object->levelset.phi.Valid_Index(my_cell_index+offset)?body.implicit_object->levelset.phi(my_cell_index+offset):1;
        fragment_implicit_object->levelset.phi(iterator.index)=max(implicit_object->levelset.phi(my_cell_index),body_phi/phi_scale);
        if(fragment_implicit_object->levelset.phi(iterator.index)<=0) colors(iterator.index)=0;
        else colors(iterator.index)=-1;}
    // Separate regions
    FLOOD_FILL_3D flood_fill;ARRAY<bool,FACE_INDEX<3> > edge_is_blocked(fragment_implicit_object->levelset.grid);edge_is_blocked.Fill(false);
    int num_colors=flood_fill.Flood_Fill(colors,edge_is_blocked);
    ARRAY<RANGE<TV_INT> > region_counts;
    region_counts.Resize(num_colors);
    for(int region=1;region<=num_colors;region++) region_counts(region)=RANGE<TV_INT>::Zero_Box().Thickened(-INT_MAX);
    if(num_colors>1)
        for(NODE_ITERATOR iterator(fragment_implicit_object->levelset.grid);iterator.Valid();iterator.Next()){
            int color=colors(iterator.index);
            if(color>0) region_counts(color).Enlarge_To_Include_Point(iterator.index);}
    else region_counts(1)=fragment_implicit_object->levelset.grid.Domain_Indices();
    // for each color, construct the particles and fracture region correctly
    for(int color=1;color<=num_colors;color++){
        LEVELSET_IMPLICIT_OBJECT<TV>* region_implicit_object;
        if(num_colors==1) region_implicit_object=fragment_implicit_object;
        else{
            // Make a new implicit object for the color
            if(!region_counts(color).Size()) continue;
            region_implicit_object=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
            region_counts(color)=RANGE<TV_INT>::Intersect(fragment_implicit_object->levelset.grid.Domain_Indices(),region_counts(color).Thickened(2));
            TV min_corner=fragment_implicit_object->levelset.grid.Node(region_counts(color).min_corner);
            TV max_corner=fragment_implicit_object->levelset.grid.Node(region_counts(color).max_corner);
            RANGE<TV> region_domain=RANGE<TV>(min_corner,max_corner);
            region_implicit_object->levelset.grid.Initialize(region_counts(color).Edge_Lengths()+TV_INT::All_Ones_Vector(),region_domain,implicit_object->levelset.grid.Is_MAC_Grid());
            region_implicit_object->levelset.phi.Resize(region_implicit_object->levelset.grid.Domain_Indices());
            region_implicit_object->Update_Box();region_implicit_object->Update_Minimum_Cell_Size();
            for(NODE_ITERATOR iterator(region_implicit_object->levelset.grid);iterator.Valid();iterator.Next()){
                TV_INT old_index=iterator.index+region_counts(color).min_corner-TV_INT::All_Ones_Vector();
                if(colors(old_index)==color) region_implicit_object->levelset.phi(iterator.index)=fragment_implicit_object->levelset.phi(old_index);
                else region_implicit_object->levelset.phi(iterator.index)=abs(fragment_implicit_object->levelset.phi(old_index));}}

        // Initialize the new region's particles.  The particles.X are in this region's object space
        TRIANGULATED_SURFACE<T>* region_triangulated_surface;
        if(tessellate_region) region_triangulated_surface=TESSELLATION::Generate_Triangles(*region_implicit_object);
        else region_triangulated_surface=TRIANGULATED_SURFACE<T>::Create();
        int initial_number_of_particles=region_triangulated_surface->particles.array_collection->Size();

        particle_intersection_thickness=implicit_object->levelset.grid.dX.Min();
        MATRIX<T,3> o2bl_RS=body.levelset_RS.Inverse()*object_RS;
        TV o2bl_T=body.levelset_RS.Solve_Linear_System(object_T-body.levelset_T);

        MATRIX<T,3> bo2l_RS=levelset_RS.Inverse()*body.object_RS;
        TV bo2l_T=levelset_RS.Solve_Linear_System(body.object_T-levelset_T);
        MATRIX<T,3> bo2o_RS=object_RS.Inverse()*body.object_RS;
        TV bo2o_T=object_RS.Solve_Linear_System(body.object_T-object_T);

        if(!use_particle_optimization){
            for(int p=1;p<=triangulated_surface->particles.array_collection->Size();p++){ // clamp our particles
                TV test_point=o2bl_RS*triangulated_surface->particles.X(p)+o2bl_T;
                if((body.implicit_object->levelset.Extended_Phi(test_point)/phi_scale-particle_intersection_thickness)<=0 && 
                   region_implicit_object->levelset.grid.domain.Lazy_Inside(extra_levelset_frame.Inverse_Times(triangulated_surface->particles.X(p)))){
                    int added_particle=region_triangulated_surface->particles.array_collection->Add_Element();
                    region_triangulated_surface->particles.X(added_particle)=triangulated_surface->particles.X(p);}}
            for(int p=1;p<=body.triangulated_surface->particles.array_collection->Size();p++){ // clamp body particles
                TV test_point=bo2l_RS*body.triangulated_surface->particles.X(p)+bo2l_T;
                if(implicit_object->levelset.Extended_Phi(test_point)-particle_intersection_thickness<=0){
                    TV point_to_add=bo2o_RS*body.triangulated_surface->particles.X(p)+bo2o_T;
                    if(region_implicit_object->levelset.grid.domain.Lazy_Inside(extra_levelset_frame.Inverse_Times(point_to_add))){
                        int added_particle=region_triangulated_surface->particles.array_collection->Add_Element();
                        region_triangulated_surface->particles.X(added_particle)=point_to_add;}}}}
        else{
            // Iterate over the new levelset's cells, and find the corresponding cell in both the levelsets. Grab all particles from those cells.
            for(NODE_ITERATOR iterator(region_implicit_object->levelset.grid);iterator.Valid();iterator.Next()){
                TV_INT old_index=iterator.index+region_counts(color).min_corner-TV_INT::All_Ones_Vector();
                if(region_implicit_object->levelset.phi(iterator.index)-particle_intersection_thickness<=0){
                    TV_INT my_cell_index=old_index+fragment_counts.min_corner-TV_INT::All_Ones_Vector();
                    for(int p=1;p<=particle_partition->partition(my_cell_index).m;p++){
                        int added_particle=region_triangulated_surface->particles.array_collection->Add_Element();
                        region_triangulated_surface->particles.X(added_particle)=triangulated_surface->particles.X(particle_partition->partition(my_cell_index)(p));}
                    for(int p=1;p<=body.particle_partition->partition(my_cell_index+offset).m;p++){
                        TV point_to_add=bo2o_RS*body.triangulated_surface->particles.X(body.particle_partition->partition(my_cell_index+offset)(p))+bo2o_T;
                        int added_particle=region_triangulated_surface->particles.array_collection->Add_Element();
                        region_triangulated_surface->particles.X(added_particle)=point_to_add;}}}}

        int number_of_particles_to_swap=min(initial_number_of_particles,region_triangulated_surface->particles.array_collection->Size()-initial_number_of_particles);
        ARRAY<int> particle_map(region_triangulated_surface->particles.array_collection->Size());
        for(int i=1;i<=particle_map.m;++i) particle_map(i)=i;
        for(int index_to_swap=1;index_to_swap<=number_of_particles_to_swap;++index_to_swap)
            exchange(particle_map(index_to_swap),particle_map(particle_map.m-index_to_swap+1));
        for(int index=1;index<=number_of_particles_to_swap;++index)
            exchange(region_triangulated_surface->particles.X(index),region_triangulated_surface->particles.X(particle_map(index)));
        for(int t=1;t<=region_triangulated_surface->mesh.elements.m;t++)
            region_triangulated_surface->mesh.elements(t)=particle_map.Subset(region_triangulated_surface->mesh.elements(t));
        
        region_triangulated_surface->Update_Number_Nodes();
        if(tessellate_region) region_triangulated_surface->mesh.Initialize_Adjacent_Elements();
        if(!initial_number_of_particles) {delete region_implicit_object;delete region_triangulated_surface;continue;}
        new_regions.Append(new FRACTURE_REGION<T>(region_triangulated_surface,region_implicit_object,use_particle_optimization));}

    if(num_colors>1) delete fragment_implicit_object;
    return new_regions;
}
//#####################################################################
// Function Compute_Volume
//#####################################################################
template<class T> T FRACTURE_REGION<T>::
Compute_Volume() const
{
    ARRAY_VIEW<T>& phi_array=implicit_object->levelset.phi.array;
    int num_inside=0;
    for(int i=1;i<=phi_array.Size();i++) if(phi_array(i)<=0) num_inside++;
    return num_inside*implicit_object->levelset.grid.Cell_Size();
}
//#####################################################################
// Function Compute_Inertial_Properties
//#####################################################################
template<class T> void FRACTURE_REGION<T>::
Compute_Inertial_Properties(const T density,TV& com,T& mass,T_WORLD_SPACE_INERTIA_TENSOR& inertia) const
{
    int num_inside=0;
    com=TV();
    for(NODE_ITERATOR iterator(implicit_object->levelset.grid);iterator.Valid();iterator.Next()){
        if(implicit_object->levelset.phi(iterator.index)>0) continue; // ignore if outside
        num_inside++;
        com+=iterator.Location();}
    com/=(T)num_inside;com=extra_levelset_frame*com;
    T cell_mass=implicit_object->levelset.grid.Cell_Size()*density;
    mass=cell_mass*num_inside;

    SYMMETRIC_MATRIX<T,3> moments;
    for(NODE_ITERATOR iterator(implicit_object->levelset.grid);iterator.Valid();iterator.Next()){
        if(implicit_object->levelset.phi(iterator.index)>0) continue; // ignore if outside
        moments+=SYMMETRIC_MATRIX<T,3>::Outer_Product(extra_levelset_frame*iterator.Location()-com);}
    inertia=cell_mass*((moments.Trace()+num_inside*sqr(implicit_object->levelset.grid.dX.x)/6)-moments);
}
template<class T> void FRACTURE_REGION<T>::
Read(TYPED_ISTREAM& input)
{
    Read_Binary(input,fracture_offset);
    triangulated_surface=new TRIANGULATED_SURFACE<T>(*new TRIANGLE_MESH(),*new GEOMETRY_PARTICLES<TV>());
    implicit_object=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    Read_Binary(input,*triangulated_surface);
    Read_Binary(input,*implicit_object);
}
template<class T> void FRACTURE_REGION<T>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,fracture_offset);
    Write_Binary(output,*triangulated_surface);
    Write_Binary(output,*implicit_object);
}
template<class T> void FRACTURE_REGION<T>::
Initialize_Particle_Partition()
{
    if(particle_partition) delete particle_partition;
    particle_partition=new PARTICLE_PARTITION<TV>(implicit_object->levelset.grid.domain,implicit_object->levelset.grid.counts,triangulated_surface->particles,false,false);
}
//#####################################################################
template class FRACTURE_REGION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_REGION<double>;
#endif
