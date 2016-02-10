#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_RLE.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_RLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_RLE<T_GRID>::
PARTICLE_LEVELSET_RLE(T_GRID& grid_input)
    :PARTICLE_LEVELSET<T_GRID>(grid_input,phi),grid(grid_input),use_removed_negative_particles_in_long_cells(false),use_removed_positive_particles_in_long_cells(false),
    removed_negative_particles_in_long_cells(0),removed_positive_particles_in_long_cells(0),mpi_grid(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_RLE<T_GRID>::
~PARTICLE_LEVELSET_RLE()
{
    delete removed_negative_particles_in_long_cells;delete removed_positive_particles_in_long_cells;
}
//#####################################################################
// Function Initialize_Particle_Levelset_Grid_Values
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Initialize_Particle_Levelset_Grid_Values()
{
    PARTICLE_LEVELSET<T_GRID>::Initialize_Particle_Levelset_Grid_Values();phi.Resize(grid.number_of_cells);
    if(!removed_negative_particles_in_long_cells) removed_negative_particles_in_long_cells=template_removed_particles.Clone();
    if(!removed_positive_particles_in_long_cells) removed_positive_particles_in_long_cells=template_removed_particles.Clone();
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Seed_Particles(const T time,const bool verbose)
{
    bool normals_defined=levelset.normals!=0;levelset.Compute_Normals(); // make sure normals are accurate

    int n=Reseed_Add_Particles(negative_particles,positive_particles,-1,PARTICLE_LEVELSET_NEGATIVE,time,0);
    for(int b=1;b<=negative_particles.m;b++)if(negative_particles(b)) negative_particles(b)->array_collection->Compact();
    int p=Reseed_Add_Particles(positive_particles,negative_particles,1,PARTICLE_LEVELSET_POSITIVE,time,0);
    for(int b=1;b<=negative_particles.m;b++)if(positive_particles(b)) positive_particles(b)->array_collection->Compact(); // TODO: revisit
    if(verbose) LOG::cout << n  << " negative particles & " << p << " positive particles" << std::endl;

    if(!normals_defined){delete levelset.normals;levelset.normals=0;}
}
//#####################################################################
// Function Attract_Individual_Particle_To_Interface
//#####################################################################
template<class T_GRID> bool PARTICLE_LEVELSET_RLE<T_GRID>::
Attract_Individual_Particle_To_Interface_And_Adjust_Radius(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const T phi_min,const T phi_max,const BLOCK_ITERATOR& block,
    const int index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(block.Block());
    T phi=levelset.Phi(block,cell_particles.X(index));bool inside=(phi >= phi_min && phi <= phi_max);
    RANGE<TV> box=block.Bounding_Box();
    if(!inside){
        T phi_goal=random.Get_Uniform_Number(phi_min,phi_max);int iteration=0;
        while(!inside && iteration < maximum_iterations_for_attraction){
            TV N=levelset.Normal(block,cell_particles.X(index)),X;T distance=phi_goal-phi,dt=1;bool inside_domain=false;
            while(!inside_domain && iteration<maximum_iterations_for_attraction){
                X=cell_particles.X(index)+dt*distance*N;
                if(grid.domain.Lazy_Outside(X)){dt*=.5;iteration++;}else inside_domain=true;}
            if(!inside_domain) break; // ran out of iterations
            phi=levelset.Phi(block,X);inside=(phi >= phi_min && phi <= phi_max);
            if(!inside){
                dt*=.5;iteration++;X=cell_particles.X(index)+dt*distance*N;
                phi=levelset.Phi(block,X);inside=(phi >= phi_min && phi <= phi_max);}
            cell_particles.X(index)=X;}}
    if(!inside || box.Lazy_Outside(cell_particles.X(index))){cell_particles.array_collection->Delete_Element(index);return false;}
    TV V_temp;levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,index,V_temp,particle_type,0,time);
    if(box.Lazy_Outside(cell_particles.X(index))){cell_particles.array_collection->Delete_Element(index);return false;}
    phi=levelset.Collision_Aware_Phi(block,cell_particles.X(index));
    if(phi<phi_min || phi>phi_max){Delete_Particle(cell_particles,index);return false;}
    cell_particles.radius(index)=clamp(abs(phi),minimum_particle_radius,maximum_particle_radius);
    return true;
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Adjust_Particle_Radii()
{
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(negative_particles(b))
        Adjust_Particle_Radii(block,*negative_particles(b),-1);}
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(positive_particles(b))
        Adjust_Particle_Radii(block,*positive_particles(b),1);}
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Adjust_Particle_Radii(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign)
{
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    // new radius is negative if the particle is on the wrong side of the interface
    if(sign==1) for(int k=1;k<=particles.array_collection->Size();k++)particles.radius(k)=clamp(levelset.Phi(block,particles.X(k)),minimum_particle_radius,maximum_particle_radius);
    else for(int k=1;k<=particles.array_collection->Size();k++) particles.radius(k)=clamp(-levelset.Phi(block,particles.X(k)),minimum_particle_radius,maximum_particle_radius);
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Modify_Levelset_Using_Escaped_Particles()
{
    levelset.Compute_Block_Minimum_And_Maximum();
    if(bias_towards_negative_particles){
        Modify_Levelset_Using_Escaped_Particles(phi,positive_particles,1);Modify_Levelset_Using_Escaped_Particles(phi,negative_particles,-1);}
    else{
        ARRAY<T> phi_minus(phi),phi_plus(phi);
        Modify_Levelset_Using_Escaped_Particles(phi_minus,negative_particles,-1);Modify_Levelset_Using_Escaped_Particles(phi_plus,positive_particles,1);
        for(int c=1;c<=grid.number_of_cells;c++)if(abs(levelset.phi(c)) <= half_band_width) levelset.phi(c)=minmag(phi_minus(c),phi_plus(c));}
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Modify_Levelset_Using_Escaped_Particles(ARRAY<T>& phi,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const int sign)
{
    T radius_multiplier=-sign*outside_particle_distance_multiplier,one_over_radius_multiplier=1/radius_multiplier;
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particles(b)){
        bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
        for(PARTICLE_LEVELSET_PARTICLES<TV>* link=particles(b);link;link=link->next){PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*link;
            for(int k=1;k<=cell_particles.array_collection->Size();k++)if(levelset.Phi(block,cell_particles.X(k))*one_over_radius_multiplier > cell_particles.radius(k)){
                for(int ii=0;ii<T_GRID::number_of_cells_per_block;ii++){int c=block.Cell(ii);
                    T radius_minus_sign_phi=cell_particles.radius(k)-sign*phi(c);
                    if(radius_minus_sign_phi > 0){TV center=block.Cell_X(ii);
                        T distance_squared=(center-cell_particles.X(k)).Magnitude_Squared();
                        if(distance_squared < sqr(radius_minus_sign_phi)){
                            static COLLISION_GEOMETRY_ID body_id;static int aggregate_id;static TV intersection_point;
                            if(near_objects && levelset.collision_body_list->collision_geometry_collection.Intersection_Between_Points(center,cell_particles.X(k),body_id,aggregate_id,intersection_point)) continue;
                            phi(c)=sign*(cell_particles.radius(k)-sqrt(distance_squared));}}}}}
        if(near_objects) levelset.Disable_Collision_Aware_Interpolation();}}
}
//#####################################################################
// Function Update_Particle_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Update_Particle_Cells(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles)
{
    if(mpi_grid) Exchange_Boundary_Particles(*mpi_grid,*this,particles,(PARTICLE_LEVELSET_PARTICLES<TV>*)0,(int)(cfl_number+(T)1.5));

    ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*> particle_links;particle_links.Preallocate(10);
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particles(b)){
        RANGE<TV> box=block.Bounding_Box();PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(b);
        particle_links.Remove_All();for(PARTICLE_LEVELSET_PARTICLES<TV>* link=particles(b);link;link=link->next)particle_links.Append(link);
        for(int link=particle_links.m;link>=1;link--)for(int k=particle_links(link)->array_collection->Size();k>=1;k--)if(box.Lazy_Outside(particle_links(link)->X(k))){
            int final_b=grid.Clamped_Block_Index(particle_links(link)->X(k),1);if(!final_b){Delete_Particle(*particle_links(link),k);continue;}
            if(!particles(final_b)) particles(final_b)=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
            int absolute_index=(link-1)*particle_pool.number_particles_per_cell+k;
            Move_Particle(cell_particles,*particles(final_b),absolute_index);}}}
    for(int b=1;b<=particles.m;b++)if(particles(b) && !particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));
}
//#####################################################################
// Function Update_Particle_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Update_Particle_Cells(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles_in_long_cells,const T dt)
{
    if(mpi_grid){
        T maximum_speed_squared=0;
        for(int b=1;b<=grid.number_of_blocks;b++) if(particles(b)) maximum_speed_squared=max(maximum_speed_squared,particles(b)->V.Magnitude_Squared());
        T maximum_distance_in_cells=dt*sqrt(maximum_speed_squared)/grid.Minimum_Edge_Length();
        Exchange_Boundary_Particles(*mpi_grid,*this,particles,particles_in_long_cells,(int)(maximum_distance_in_cells+(T)1.5));}

    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particles(b)){
        RANGE<TV> box=block.Bounding_Box();PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& cell_particles=*particles(b);
        for(int k=cell_particles.array_collection->Size();k>=1;k--)if(box.Lazy_Outside(cell_particles.X(k))){
            int final_b=grid.Clamped_Block_Index(cell_particles.X(k),1);
            if(!final_b){
                if(particles_in_long_cells) Move_Particle(cell_particles,*particles_in_long_cells,k);
                else cell_particles.array_collection->Delete_Element(k);}
            else{
                if(!particles(final_b)) particles(final_b)=Allocate_Particle<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >();
                Move_Particle(cell_particles,*particles(final_b),k);}}}}
    if(particles_in_long_cells)
        for(int k=particles_in_long_cells->array_collection->Size();k>=1;k--){
            int final_b=grid.Clamped_Block_Index(particles_in_long_cells->X(k),1);if(!final_b) continue;
            if(!particles(final_b)) particles(final_b)=Allocate_Particle<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >();
            Move_Particle(*particles_in_long_cells,*particles(final_b),k);}
    LOG::cout<<"particles in long cells = "<<particles_in_long_cells<<std::endl;
    for(int b=1;b<=particles.m;b++)if(particles(b) && !particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));
}
//#####################################################################
// Function Compact_Particles_Into_Single_Particle_Bin
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Compact_Particles_Into_Single_Particle_Bin()
{
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();
        if(negative_particles(b)) Compact_Particles_Into_Single_Particle_Bin(block,*negative_particles(b),-1);
        if(positive_particles(b)) Compact_Particles_Into_Single_Particle_Bin(block,*positive_particles(b),1);}
}
//#####################################################################
// Function Compact_Particles_Into_Single_Particle_Bin
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Compact_Particles_Into_Single_Particle_Bin(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign)
{
    if(!particles.next) return;
    assert(particles.array_collection->Size()==particle_pool.number_particles_per_cell);
    int heap_size=0;ARRAY<int> heap_particle_indices(particle_pool.number_particles_per_cell);ARRAY<T> heap_phi_minus_radius(particle_pool.number_particles_per_cell);
    for(int k=particles.array_collection->Size();k>=1;k--){
        heap_size++;heap_particle_indices(heap_size)=k;heap_phi_minus_radius(heap_size)=sign*levelset.Phi(block,particles.X(k))-particles.radius(k);}
    ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices); // when heap is full, order values with largest on top 
    PARTICLE_LEVELSET_PARTICLES<TV>* particles_link=&particles;
    while(particles_link->next){
        particles_link=particles_link->next;assert(!particles_link->next||particles_link->array_collection->Size()==particle_pool.number_particles_per_cell);
        for(int k=particles_link->array_collection->Size();k>=1;k--){
            T phi_minus_radius=sign*levelset.Phi(block,particles_link->X(k))-particles_link->radius(k);
            if(phi_minus_radius < heap_phi_minus_radius(1)){ // delete particle on top of heap & add new particle
                particles.array_collection->Copy_Element(*particles_link->array_collection,k,heap_particle_indices(1));
                heap_phi_minus_radius(1)=phi_minus_radius;
                ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices,1,heap_phi_minus_radius.m);}}}
    Free_Particle_And_Clear_Pointer(particles.next);
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Euler_Step_Particles(const ARRAY<T>& V,const T dt,const T time,const bool use_second_order_for_nonremoved_particles,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    if(use_second_order_for_nonremoved_particles){
        Second_Order_Runge_Kutta_Step_Particles(V,negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,time,update_particle_cells_after_euler_step,verbose);
        Second_Order_Runge_Kutta_Step_Particles(V,positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,time,update_particle_cells_after_euler_step,verbose);}
    else{
        Euler_Step_Particles(V,negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,time,update_particle_cells_after_euler_step);
        Euler_Step_Particles(V,positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,time,update_particle_cells_after_euler_step);}

    Euler_Step_Removed_Particles(dt,time,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Euler_Step_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    if(use_removed_negative_particles)
        Euler_Step_Removed_Particles(removed_negative_particles,use_removed_negative_particles_in_long_cells?removed_negative_particles_in_long_cells:0,
            PARTICLE_LEVELSET_REMOVED_NEGATIVE,dt,time,update_particle_cells_after_euler_step,verbose);
    if(use_removed_positive_particles)
        Euler_Step_Removed_Particles(removed_positive_particles,use_removed_positive_particles_in_long_cells?removed_positive_particles_in_long_cells:0,
            PARTICLE_LEVELSET_REMOVED_POSITIVE,dt,time,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Euler_Step_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Euler_Step_Removed_Particles(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles_in_long_cells,
    const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    int number_of_deleted_particles=0,number_of_non_occupied_cells=0;
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particles(b)){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& cell_particles=*particles(b);
        if(!levelset.collision_body_list->Swept_Occupied_Block(block)){
            number_of_non_occupied_cells++;
            for(int k=cell_particles.array_collection->Size();k>=1;k--) // since not an occupied cell, don't need to adjust for objects
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,cell_particles.V(k),particle_type,dt,time);}
        else{
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,cell_particles.V(k),particle_type,dt,time);
                T collision_distance=Particle_Collision_Distance(cell_particles.quantized_collision_distance(k));
                if(!Adjust_Particle_For_Objects(cell_particles.X(k),cell_particles.V(k),cell_particles.radius(k),collision_distance,particle_type,dt,time)){
                    cell_particles.array_collection->Delete_Element(k);number_of_deleted_particles++;}}}
        cell_particles.Euler_Step_Position(dt);}}
    // handle particles in long cells
    if(particles_in_long_cells){
        for(int k=particles_in_long_cells->array_collection->Size();k>=1;k--){
            levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(*particles_in_long_cells,k,particles_in_long_cells->V(k),particle_type,dt,time);
            T collision_distance=Particle_Collision_Distance(particles_in_long_cells->quantized_collision_distance(k));
            // TODO: make this always check nonrasterized objects
            if(!Adjust_Particle_For_Objects(particles_in_long_cells->X(k),particles_in_long_cells->V(k),particles_in_long_cells->radius(k),collision_distance,particle_type,dt,time)){
                particles_in_long_cells->array_collection->Delete_Element(k);number_of_deleted_particles++;}}
        particles_in_long_cells->Euler_Step_Position(dt);}
    if(verbose){
        if(number_of_deleted_particles) LOG::cout<<"Deleted "<<number_of_deleted_particles<<" "<<PARTICLE_LEVELSET<T_GRID>::Particle_Type_Name(particle_type)<<" due to crossover"<<std::endl;
        if(number_of_non_occupied_cells) LOG::cout<<"Skipped "<<number_of_non_occupied_cells<<" non occupied cells"<<std::endl;}
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles,particles_in_long_cells,dt);
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Euler_Step_Particles(const ARRAY<T>& V,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
    const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool assume_particles_in_correct_cells)
{
    PHYSBAM_NOT_IMPLEMENTED();
/*
    assert(assume_particles_in_correct_cells);
    //T object_contour=T_GRID::ITERATOR_CELL::Short_Max_Length(grid)*(this->min_collision_distance_factor+cfl_number);
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particles(b)){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(b);
        typename T_GRID::LINEAR_INTERPOLATION_MAC_HELPER block_velocity(block,V);
        //bool check_objects=(*levelset_object.block_minimum)(b)<=object_contour;
        for(int k=1;k<=cell_particles.array_collection->Size();k++){
            TV velocity=block_velocity.Interpolate_Face(cell_particles.X(k));
            PHYSBAM_NOT_IMPLEMENTED();//if(check_objects) levelset.levelset_callbacks->Adjust_Particle_For_Objects(cell_particles,k,velocity,particle_type,dt,time);
            cell_particles.X(k)+=dt*velocity;}}}
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
*/
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Second_Order_Runge_Kutta_Step_Particles(const ARRAY<T>& V,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
    const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    T_FACE_LOOKUP V_lookup(V);
    T_FACE_LOOKUP_COLLIDABLE V_lookup_collidable(V_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP V_lookup_collidable_lookup(V_lookup_collidable,V_lookup);
    T_LINEAR_INTERPOLATION_SCALAR linear_interpolation; // use for second step since particle may not be in initial block
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    COLLISION_GEOMETRY_ID body_id;int aggregate_id;T start_phi,end_phi;TV body_normal,body_velocity;int number_of_deleted_particles=0,number_of_non_occupied_cells=0;
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particles(b)){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(b);
        T_LINEAR_INTERPOLATION_MAC_HELPER block_velocity(block,V);RANGE<TV> box=block.Bounding_Box();
        if(!levelset.collision_body_list->Swept_Occupied_Block(block)){ // not occupied, so advect ignoring objects
            number_of_non_occupied_cells++;
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                TV velocity=block_velocity.Interpolate_Face(cell_particles.X(k));
                TV X_new=cell_particles.X(k)+dt*velocity;
                if(box.Lazy_Inside(X_new)) velocity=(T).5*(velocity+block_velocity.Interpolate_Face(X_new));
                else{
                    const BLOCK_ITERATOR block(grid,X_new);if(!block){cell_particles.X(k)=X_new;continue;} // particle will be deleted in cell update
                    velocity=(T).5*(velocity+T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face(block,V,X_new));}
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
                cell_particles.X(k)+=dt*velocity;}}
        else{ // collision aware advection
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                T collision_distance=Particle_Collision_Distance(cell_particles.quantized_collision_distance(k));
                BLOCK_ITERATOR block_new(block);RANGE<TV> box_new=box;
                // first push particle out
                if(particle_type==PARTICLE_LEVELSET_NEGATIVE){bool particle_crossover;
                    if(levelset.collision_body_list->Push_Out_Point(cell_particles.X(k),collision_distance,true,particle_crossover)&&particle_crossover){
                        cell_particles.array_collection->Delete_Element(k);number_of_deleted_particles++;continue;} // delete due to crossover in push out
                    if(box_new.Lazy_Outside(cell_particles.X(k))){
                        block_new.Initialize(cell_particles.X(k));if(!block_new) continue; // particle will be deleted in cell update
                        box_new=block_new.Bounding_Box();}}
                TV velocity=linear_interpolation_collidable.From_Block_Face(block_new,V_lookup_collidable_lookup,cell_particles.X(k));
                // adjust normal component of velocity to not move into nearest body
                bool got_interaction_with_body=false;
                if(particle_type==PARTICLE_LEVELSET_NEGATIVE)
                    got_interaction_with_body=levelset.collision_body_list->Get_Body_Penetration(cell_particles.X(k),cell_particles.X(k),collision_distance,dt,body_id,aggregate_id,
                        start_phi,end_phi,body_normal,body_velocity);
                if(got_interaction_with_body){
                    T relative_normal_velocity=TV::Dot_Product(velocity-body_velocity,body_normal);
                    if(relative_normal_velocity<0) velocity-=relative_normal_velocity*body_normal;}
                // compute position using velocity but clamp if intersects any body
                TV X_new=cell_particles.X(k)+dt*velocity;
                RAY<TV> ray;
                if(RAY<TV>::Create_Non_Degenerate_Ray(cell_particles.X(k),X_new-cell_particles.X(k),ray))
                    if(levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id))
                        X_new=ray.Point(ray.t_max);
                // get average velocity
                block_new.Initialize(X_new);if(!block_new) continue; // particle will be deleted in cell update
                velocity=(T).5*(velocity+linear_interpolation_collidable.From_Block_Face(block_new,V_lookup_collidable_lookup,X_new));
                // adjust normal component of velocity again
                if(got_interaction_with_body){
                    T relative_normal_velocity=TV::Dot_Product(velocity-body_velocity,body_normal);
                    if(relative_normal_velocity<0) velocity-=relative_normal_velocity*body_normal;}
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
                if(!Adjust_Particle_For_Objects(cell_particles.X(k),velocity,cell_particles.radius(k),collision_distance,particle_type,dt,time)){
                    cell_particles.array_collection->Delete_Element(k);number_of_deleted_particles++;}
                else cell_particles.X(k)+=dt*velocity;}}}}
    if(verbose){
        if(number_of_deleted_particles) LOG::cout<<"Deleted "<<number_of_deleted_particles<<" "<<PARTICLE_LEVELSET<T_GRID>::Particle_Type_Name(particle_type)<<" due to crossover"<<std::endl;
        if(number_of_non_occupied_cells) LOG::cout<<"Skipped "<<number_of_non_occupied_cells<<" non occupied cells"<<std::endl;}
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_RLE<T_GRID>::
Reseed_Particles(const T time,ARRAY<bool>* cell_centered_mask)
{
    int new_particles=0;
    bool normals_defined=levelset.normals!=0;levelset.Compute_Normals(); // make sure normals are accurate

    if(!cell_centered_mask) new_particles-=Reseed_Delete_Particles(negative_particles,-1);
    new_particles+=Reseed_Add_Particles(negative_particles,positive_particles,-1,PARTICLE_LEVELSET_NEGATIVE,time,cell_centered_mask);
    if(!cell_centered_mask) for(int b=1;b<=negative_particles.m;b++)if(negative_particles(b)) negative_particles(b)->array_collection->Compact();

    if(!cell_centered_mask) new_particles-=Reseed_Delete_Particles(positive_particles,1);
    new_particles+=Reseed_Add_Particles(positive_particles,negative_particles,1,PARTICLE_LEVELSET_POSITIVE,time,cell_centered_mask);
    if(!cell_centered_mask) for(int b=1;b<=positive_particles.m;b++)if(positive_particles(b)) positive_particles(b)->array_collection->Compact();

    if(!normals_defined){delete levelset.normals;levelset.normals=0;}
    return new_particles;
}
//#####################################################################
// Function Reseed_Delete_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_RLE<T_GRID>::
Reseed_Delete_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const int sign)
{
    levelset.Compute_Block_Minimum_And_Maximum();
    int number_deleted=0;ARRAY<int> heap_particle_indices(number_particles_per_cell);ARRAY<T> heap_phi_minus_radius(number_particles_per_cell);
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particles(b)){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(b);assert(!cell_particles.next);
        for(int i=0;i<T_GRID::number_of_cells_per_block;i++){
            T unsigned_phi=sign*levelset.phi(block.Cell(i));if(0<=unsigned_phi && unsigned_phi<=half_band_width) goto NEAR_THE_INTERFACE;}

        if(cell_particles.array_collection->Size() == 0) Free_Particle_And_Clear_Pointer(particles(b));
        else{ // delete all non-escaped particles
            ARRAY<bool> escaped;Identify_Escaped_Particles(block,cell_particles,escaped,sign);
            for(int index=cell_particles.array_collection->Size();index>=1;index--)if(!escaped(index)) cell_particles.array_collection->Delete_Element(index);
            if(cell_particles.array_collection->Size() == 0) Free_Particle_And_Clear_Pointer(particles(b));}
        continue;

        NEAR_THE_INTERFACE:; // can only get in here via the goto
        if(cell_particles.array_collection->Size() <= number_particles_per_cell) continue;
        ARRAY<bool> escaped;Identify_Escaped_Particles(block,cell_particles,escaped,sign);int number_of_escaped_particles=escaped.Number_True();
        int total_particles=cell_particles.array_collection->Size()-number_of_escaped_particles;
        if(total_particles > number_particles_per_cell){ // too many particles - delete particles with a heap sort
            number_deleted+=total_particles-number_particles_per_cell;int heap_size=0;
            for(int index=1;index<=cell_particles.array_collection->Size();index++)if(!escaped.m || !escaped(index)){
                T phi_minus_radius=sign*levelset.Phi(cell_particles.X(index))-cell_particles.radius(index);
                if(heap_size < number_particles_per_cell){ // add particle to heap
                    heap_size++;heap_particle_indices(heap_size)=index;heap_phi_minus_radius(heap_size)=phi_minus_radius;
                    if(heap_size==number_particles_per_cell) ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices);} // when heap is full, order values with largest on top
                else{ // excess particles don't fit in the heap
                    if(phi_minus_radius < heap_phi_minus_radius(1)){ // delete particle on top of heap & add new particle
                        cell_particles.array_collection->Add_To_Deletion_List(heap_particle_indices(1));
                        heap_phi_minus_radius(1)=phi_minus_radius;heap_particle_indices(1)=index;
                        ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices,1,heap_phi_minus_radius.m);}
                    else cell_particles.array_collection->Add_To_Deletion_List(index);}} // delete new particle, larger than top of heap
            cell_particles.array_collection->Delete_Elements_On_Deletion_List();}}}
    return number_deleted;
}
//#####################################################################
// Function Reseed_Add_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_RLE<T_GRID>::
Reseed_Add_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& other_particles,
    const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,ARRAY<bool>* cell_centered_mask)
{
    int number_added=0;
    for(BLOCK_ITERATOR block(grid,1);block;block++){
        for(int i=0;i<T_GRID::number_of_cells_per_block;i++){
            T unsigned_phi=sign*levelset.phi(block.Cell(i));
            if(0<=unsigned_phi && unsigned_phi<=half_band_width && (!cell_centered_mask||(*cell_centered_mask)(block.Cell(i)))) goto NEAR_THE_INTERFACE;}
        continue;

        NEAR_THE_INTERFACE:; // can only get in here via the goto
        int b=block.Block();
        if(!particles(b)) particles(b)=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
        int total_particles=particles(b)->array_collection->Size(),total_other_particles=0;if(other_particles(b)) total_other_particles=other_particles(b)->array_collection->Size();
        if(total_particles+total_other_particles >= number_particles_per_cell) continue;
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(b);assert(!cell_particles.next);
        int number_of_particles_to_add=number_particles_per_cell-total_particles-total_other_particles;
        if(sign == -1){ // we add the negative particles first, and don't want to add too many of them...
            if(total_other_particles) number_of_particles_to_add=(int)((T)number_of_particles_to_add*(T)total_particles/(T)(total_particles+total_other_particles)+1);
            else if(!total_particles) number_of_particles_to_add=(number_of_particles_to_add+1)/2;}
        number_added+=number_of_particles_to_add;
        cell_particles.array_collection->Preallocate(total_particles+number_of_particles_to_add);
        T phi_min=sign*minimum_particle_radius,phi_max=sign*half_band_width;if(phi_min > phi_max) exchange(phi_min,phi_max);
        RANGE<TV> box=block.Bounding_Box();
        int attempts=0;
        ARRAY_VIEW<int>* id=store_unique_particle_id?cell_particles.array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID):0;
        for(int k=1;k<=number_of_particles_to_add;k++){
            int index=cell_particles.array_collection->Add_Element();
            if(id) (*id)(index)=++last_unique_particle_id;
            cell_particles.quantized_collision_distance(index)=(unsigned short)(random.Get_Number()*USHRT_MAX);
            cell_particles.X(index)=random.Get_Uniform_Vector(box);
            if(!Attract_Individual_Particle_To_Interface_And_Adjust_Radius(particles,phi_min,phi_max,block,index,particle_type,time) && ++attempts <= 5) k--;}}
    return number_added;
}
//#####################################################################
// Function Identify_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Identify_Escaped_Particles(const int sign)
{
    levelset.Compute_Block_Minimum_And_Maximum();
    if(!sign || sign == 1){
        escaped_positive_particles.Resize(grid.number_of_blocks);
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();
            if(positive_particles(b))
                Identify_Escaped_Particles(block,*positive_particles(b),escaped_positive_particles(b),1);
            else escaped_positive_particles(b).Clean_Memory();}}
    if(!sign || sign == -1){
        escaped_negative_particles.Resize(grid.number_of_blocks);
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();
            if(negative_particles(b))
                Identify_Escaped_Particles(block,*negative_particles(b),escaped_negative_particles(b),-1);
            else escaped_negative_particles(b).Clean_Memory();}}
}
//#####################################################################
// Function Identify_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Identify_Escaped_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign)
{
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    T inverse_radius_multiplier=-sign/outside_particle_distance_multiplier;
    escaped.Resize(particles.array_collection->Size(),false,false);ARRAYS_COMPUTATIONS::Fill(escaped,false);
    for(int k=1;k<=particles.array_collection->Size();k++)if(inverse_radius_multiplier*levelset.Phi(block,particles.X(k)) > particles.radius(k)) escaped(k)=true;
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Delete_Particles_Outside_Grid()
{
    Delete_Particles_Outside_Grid(positive_particles,grid);
    Delete_Particles_Outside_Grid(negative_particles,grid);
    if(use_removed_positive_particles){
        Delete_Particles_Outside_Grid(removed_positive_particles,grid);
        Delete_Particles_Outside_Grid(*removed_positive_particles_in_long_cells);}
    if(use_removed_negative_particles){
        Delete_Particles_Outside_Grid(removed_negative_particles,grid);
        Delete_Particles_Outside_Grid(*removed_negative_particles_in_long_cells);}
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Delete_Particles_Outside_Grid(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_in_long_cells)
{
    const T_BOX_HORIZONTAL& box=grid.horizontal_grid.domain;
    for(int k=particles_in_long_cells.array_collection->Size();k>=1;k--)if(!box.Lazy_Inside(particles_in_long_cells.X(k).Horizontal_Vector())) particles_in_long_cells.array_collection->Delete_Element(k);
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void PARTICLE_LEVELSET_RLE<T_GRID>::
Delete_Particles_Outside_Grid(ARRAY<T_PARTICLES*>& particles,const RLE_GRID_2D<T>&)
{
    for(RLE_GRID_ITERATOR_BLOCK_2D<T> block(grid,RANGE<VECTOR<int,1> >(0,0));block;block++){int b=block.Block();if(particles(b)){
        for(int k=particles(b)->array_collection->Size();k>=1;k--) if(particles(b)->X(k).x < grid.domain.min_corner.x) particles(b)->array_collection->Delete_Element(k);
        if(!particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));}}
    for(RLE_GRID_ITERATOR_BLOCK_2D<T> block(grid,RANGE<VECTOR<int,1> >(grid.uniform_grid.counts.x-1,grid.uniform_grid.counts.x-1));block;block++){int b=block.Block();if(particles(b)){
        for(int k=particles(b)->array_collection->Size();k>=1;k--) if(particles(b)->X(k).x > grid.domain.max_corner.x) particles(b)->array_collection->Delete_Element(k);
        if(!particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));}}
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void PARTICLE_LEVELSET_RLE<T_GRID>::
Delete_Particles_Outside_Grid(ARRAY<T_PARTICLES*>& particles,const RLE_GRID_3D<T>&)
{
    for(RLE_GRID_ITERATOR_BLOCK_3D<T> block(grid,RANGE<VECTOR<int,2> >(0,0,0,grid.uniform_grid.counts.z-1));block;block++){int b=block.Block();if(particles(b)){
        for(int k=particles(b)->array_collection->Size();k>=1;k--) if(particles(b)->X(k).x < grid.domain.min_corner.x) particles(b)->array_collection->Delete_Element(k);
        if(!particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));}}
    for(RLE_GRID_ITERATOR_BLOCK_3D<T> block(grid,RANGE<VECTOR<int,2> >(grid.uniform_grid.counts.x-1,grid.uniform_grid.counts.x-1,0,grid.uniform_grid.counts.z-1));block;block++){int b=block.Block();if(particles(b)){
        for(int k=particles(b)->array_collection->Size();k>=1;k--) if(particles(b)->X(k).x > grid.domain.max_corner.x) particles(b)->array_collection->Delete_Element(k);
        if(!particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));}}
    for(RLE_GRID_ITERATOR_BLOCK_3D<T> block(grid,RANGE<VECTOR<int,2> >(0,grid.uniform_grid.counts.x-1,0,0));block;block++){int b=block.Block();if(particles(b)){
        for(int k=particles(b)->array_collection->Size();k>=1;k--) if(particles(b)->X(k).z < grid.domain.min_corner.z) particles(b)->array_collection->Delete_Element(k);
        if(!particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));}}
    for(RLE_GRID_ITERATOR_BLOCK_3D<T> block(grid,RANGE<VECTOR<int,2> >(0,grid.uniform_grid.counts.x-1,grid.uniform_grid.counts.z-1,grid.uniform_grid.counts.z-1));block;block++){int b=block.Block();if(particles(b)){
        for(int k=particles(b)->array_collection->Size();k>=1;k--) if(particles(b)->X(k).z > grid.domain.max_corner.z) particles(b)->array_collection->Delete_Element(k);
        if(!particles(b)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(b));}}
}
//#####################################################################
// Function Identify_And_Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Identify_And_Remove_Escaped_Particles(const ARRAY<T>& V,const T radius_fraction,const bool verbose)
{
    T_FACE_LOOKUP V_lookup(V);
    T_FACE_LOOKUP_COLLIDABLE V_lookup_collidable(V_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP V_lookup_collidable_lookup(V_lookup_collidable,V_lookup);

    Identify_Escaped_Particles();int p=0,n=0;
    if(use_removed_positive_particles)
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(escaped_positive_particles(b).m)
            p+=Remove_Escaped_Particles(block,*positive_particles(b),escaped_positive_particles(b),1,removed_positive_particles(b),V_lookup_collidable_lookup,radius_fraction);}
    else
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(escaped_positive_particles(b).m)
            p+=Remove_Escaped_Particles(block,*positive_particles(b),escaped_positive_particles(b),1,radius_fraction);}
    if(use_removed_negative_particles)
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(escaped_negative_particles(b).m)
            n+=Remove_Escaped_Particles(block,*negative_particles(b),escaped_negative_particles(b),-1,removed_negative_particles(b),V_lookup_collidable_lookup,radius_fraction);}
    else
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(escaped_negative_particles(b).m)
            n+=Remove_Escaped_Particles(block,*negative_particles(b),escaped_negative_particles(b),-1,radius_fraction);}
    LOG::cout<<"new removed positive particles = "<<p<<std::endl<<"new removed negative particles = "<<n<<std::endl;
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> template<class T_FACE_LOOKUP_LOOKUP> int PARTICLE_LEVELSET_RLE<T_GRID>::
Remove_Escaped_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles,const T_FACE_LOOKUP_LOOKUP& V,const T radius_fraction)
{
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    int deleted=0;T minus_sign_over_radius_fraction=-sign/radius_fraction;
    for(int k=particles.array_collection->Size();k>=1;k--)if(escaped(k) && minus_sign_over_radius_fraction*levelset.Phi(block,particles.X(k)) > particles.radius(k)){
        if(!removed_particles) removed_particles=Allocate_Particle<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >();
        removed_particles->array_collection->Take(*particles.array_collection,k);
        removed_particles->V.Last()=linear_interpolation_collidable.From_Block_Face(block,V,removed_particles->X.Last());}
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
    return deleted;
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_RLE<T_GRID>::
Remove_Escaped_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,const T radius_fraction)
{
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    int deleted=0;
    for(int k=particles.array_collection->Size();k>=1;k--)if(escaped(k) && -sign*levelset.Phi(block,particles.X(k)) > radius_fraction*particles.radius(k)){
        particles.array_collection->Delete_Element(k);deleted++;}
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
    return deleted;
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Reincorporate_Removed_Particles(const T radius_fraction)
{
    if(use_removed_positive_particles)
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(removed_positive_particles(b))
            Reincorporate_Removed_Particles(block,positive_particles(b),1,*removed_positive_particles(b),radius_fraction);}
    if(use_removed_negative_particles)
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(removed_negative_particles(b))
            Reincorporate_Removed_Particles(block,negative_particles(b),-1,*removed_negative_particles(b),radius_fraction);}
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Reincorporate_Removed_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>*& particles,const int sign,
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& removed_particles,const T radius_fraction)
{
    assert(!particles || !particles->next || particles->array_collection->Size()==particle_pool.number_particles_per_cell);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    T one_over_radius_multiplier=-sign/radius_fraction;
    for(int k=removed_particles.array_collection->Size();k>=1;k--)if(one_over_radius_multiplier*levelset.Phi(block,removed_particles.X(k)) < removed_particles.radius(k)){
        if(!particles) particles=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
        Move_Particle(removed_particles,*particles,k);} // maybe recalculate radius to get rid of raining effect?
    if(particles) Compact_Particles_Into_Single_Particle_Bin(block,*particles,sign);
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Delete_Particles_Far_From_Interface
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Delete_Particles_Far_From_Interface(const int discrete_band)
{
    assert(discrete_band<8);
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& cell_neighbors=grid.Short_Cell_Neighbors();
    const ARRAY<VECTOR<bool,T_GRID::dimension> >& cell_neighbors_visible=levelset.collision_body_list->cell_neighbors_visible;
    ARRAY<char> near_interface(grid.number_of_cells);
    for(int cell=1;cell<=grid.number_of_cells;cell++)for(int axis=1;axis<=T_GRID::dimension;axis++){int neighbor=cell_neighbors(cell)(2*axis);
        if(neighbor && cell_neighbors_visible(cell)(axis) && LEVELSET_UTILITIES<T>::Interface(phi(cell),phi(neighbor))){
            near_interface(cell)=near_interface(neighbor)=1;}}
    for(int distance=1;distance<=discrete_band;distance++){
        char new_mask=1<<distance,old_mask=new_mask-1;
        for(int cell=1;cell<=grid.number_of_cells;cell++)for(int axis=1;axis<=T_GRID::dimension;axis++){int neighbor=cell_neighbors(cell)(2*axis);
            if(neighbor && cell_neighbors_visible(cell)(axis) && (near_interface(cell)|near_interface(neighbor))&old_mask){
                near_interface(cell)|=new_mask;near_interface(neighbor)|=new_mask;}}}
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(negative_particles(b) || positive_particles(b)){
        for(int i=0;i<T_GRID::number_of_cells_per_block;i++)if(near_interface(block.Cell(i))) goto NEAR_INTERFACE;
        // if not near interface, delete particles 
        Free_Particle_And_Clear_Pointer(negative_particles(b));
        Free_Particle_And_Clear_Pointer(positive_particles(b));
        NEAR_INTERFACE:;}}
}
//#####################################################################
// Function Delete_Particles_In_Local_Maximum_Phi_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Delete_Particles_In_Local_Maximum_Phi_Cells()
{
    T tolerance=levelset.small_number*grid.Minimum_Edge_Length();
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& cell_neighbors=grid.Short_Cell_Neighbors();
    for(CELL_ITERATOR cell(grid,1);cell;cell++){int c=cell.Cell();
        if(phi(c)<=0 || phi(c)>=grid.Minimum_Edge_Length()) continue;
        bool local_maximum=true;
        for(int n=1;n<=T_GRID::number_of_neighbors_per_cell;n++){int neighbor=cell_neighbors(c)(n);
            if(!neighbor || phi(c)<phi(neighbor)-tolerance){local_maximum=false;break;}}
        if(!local_maximum) continue;
        TV X=cell.X();TV_INT I=cell.I();
        for(int i=1;i<=T_GRID::number_of_cells_per_block;i++){
            T_BLOCK block(grid,grid.uniform_grid.Node_Cell_Index(I,i)); // TODO: don't call a function with a name that makes no sense in this context
            if(!block) continue;int b=block.Block();
            if(positive_particles(b)){PARTICLE_LEVELSET_PARTICLES<TV>& particles=*positive_particles(b);assert(!particles.next);
                for(int k=particles.array_collection->Size();k>=1;k--)if((particles.X(k)-X).Magnitude_Squared()<=sqr(particles.radius(k))) particles.array_collection->Delete_Element(k);}}}
}
//#####################################################################
// Function Transfer_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Transfer_Particles(const T_GRID& new_grid)
{
    ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*> new_positive_particles(new_grid.number_of_blocks),new_negative_particles(new_grid.number_of_blocks);
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*> new_removed_positive_particles(use_removed_positive_particles?new_grid.number_of_blocks:0),
        new_removed_negative_particles(use_removed_negative_particles?new_grid.number_of_blocks:0);
    for(BLOCK_ITERATOR old_block(grid,1);old_block;old_block++){
        BLOCK_ITERATOR new_block(new_grid,old_block.I());if(!new_block) continue;
        int b1=old_block.Block(),b2=new_block.Block();
        new_negative_particles(b2)=negative_particles(b1);negative_particles(b1)=0;
        new_positive_particles(b2)=positive_particles(b1);positive_particles(b1)=0;
        if(use_removed_negative_particles){new_removed_negative_particles(b2)=removed_negative_particles(b1);removed_negative_particles(b1)=0;}
        if(use_removed_positive_particles){new_removed_positive_particles(b2)=removed_positive_particles(b1);removed_positive_particles(b1)=0;}}
    if(use_removed_negative_particles_in_long_cells)
        Transfer_Particles_In_Long_Cells(new_grid,removed_negative_particles,new_removed_negative_particles,*removed_negative_particles_in_long_cells);
    if(use_removed_positive_particles_in_long_cells)
        Transfer_Particles_In_Long_Cells(new_grid,removed_positive_particles,new_removed_positive_particles,*removed_positive_particles_in_long_cells);
    ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>::Exchange_Arrays(negative_particles,new_negative_particles);
    ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>::Exchange_Arrays(positive_particles,new_positive_particles);
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::Exchange_Arrays(removed_negative_particles,new_removed_negative_particles);
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::Exchange_Arrays(removed_positive_particles,new_removed_positive_particles);
    for(int b=1;b<=grid.number_of_blocks;b++){Free_Particle_And_Clear_Pointer(new_negative_particles(b));Free_Particle_And_Clear_Pointer(new_positive_particles(b));}
    new_removed_negative_particles.Delete_Pointers_And_Clean_Memory();new_removed_positive_particles.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Transfer_Particles_In_Long_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_RLE<T_GRID>::
Transfer_Particles_In_Long_Cells(const T_GRID& new_grid,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& new_particles,
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_in_long_cells)
{
    for(int k=particles_in_long_cells.array_collection->Size();k>=1;k--){
        T_BLOCK new_block(new_grid,particles_in_long_cells.X(k));if(!new_block) continue;int b=new_block.Block();
        if(!new_particles(b)) new_particles(b)=Allocate_Particle<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >();
        Move_Particle(particles_in_long_cells,*new_particles(b),k);}
    for(int b=1;b<=grid.number_of_blocks;b++)if(particles(b)){
        particles_in_long_cells.array_collection->Take(*particles(b)->array_collection);
        Free_Particle_And_Clear_Pointer(particles(b));}
}
//#####################################################################
template class PARTICLE_LEVELSET_RLE<RLE_GRID_2D<float> >;
template class PARTICLE_LEVELSET_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_RLE<RLE_GRID_2D<double> >;
template class PARTICLE_LEVELSET_RLE<RLE_GRID_3D<double> >;
#endif
#endif
