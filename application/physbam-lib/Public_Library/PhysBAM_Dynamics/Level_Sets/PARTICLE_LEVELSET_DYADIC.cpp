#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//##################################################################### 
template<class T_GRID> PARTICLE_LEVELSET_DYADIC<T_GRID>::
PARTICLE_LEVELSET_DYADIC(T_GRID& grid_input,ARRAY<T>& phi_input)
    :PARTICLE_LEVELSET<T_GRID>(grid_input,phi_input)
{}
//#####################################################################
// Destructor
//##################################################################### 
template<class T_GRID> PARTICLE_LEVELSET_DYADIC<T_GRID>::
~PARTICLE_LEVELSET_DYADIC()
{}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Seed_Particles(const T time,const bool verbose)
{
    Reseed_Add_Particles(negative_particles,positive_particles,-1,PARTICLE_LEVELSET_NEGATIVE,time,0);
    Reseed_Add_Particles(positive_particles,negative_particles,1,PARTICLE_LEVELSET_POSITIVE,time,0);
    Compact_Particles_Into_Single_Particle_Bin();
    int p=0;for(int i=1;i<=levelset.grid.number_of_cells;i++) if(positive_particles(i)) p+=positive_particles(i)->array_collection->Size();
    int n=0;for(int i=1;i<=levelset.grid.number_of_cells;i++) if(negative_particles(i)) n+=negative_particles(i)->array_collection->Size();
    if(verbose) LOG::cout<<n <<" negative particles & "<<p<<" positive particles"<<std::endl;
}
//#####################################################################
// Function Attract_Individual_Particle_To_Interface
//#####################################################################
template<class T_GRID> bool PARTICLE_LEVELSET_DYADIC<T_GRID>::
Attract_Individual_Particle_To_Interface_And_Adjust_Radius(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const T phi_min,const T phi_max,const BLOCK_DYADIC<T_GRID>& block,
    const int absolute_index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const bool delete_particles_that_leave_original_cell)
{
    int index;PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=Get_Particle_Link(*particles(block.Base_Cell()->Cell()),absolute_index,index);
    T phi=levelset.Phi(block,cell_particles.X(index));bool inside=(phi>=phi_min && phi<=phi_max);
    const CELL* current_cell=block.Base_Cell();BLOCK_DYADIC<T_GRID> current_block=block;
    if(!inside){
        T phi_goal=random.Get_Uniform_Number(phi_min,phi_max);TV X,N;int iteration=0;
        while(!inside && iteration<maximum_iterations_for_attraction){
            N=levelset.Normal(current_block,cell_particles.X(index));T distance=phi_goal-phi,dt=1;bool inside_domain=false;
            while(!inside_domain && iteration<maximum_iterations_for_attraction){
                X=cell_particles.X(index)+dt*distance*N;
                if(levelset.grid.Outside(X)){dt*=.5;iteration++;}else inside_domain=true;}
            if(!inside_domain)break; // ran out of iterations
            current_cell=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,X);
            if(!current_cell) break; // out of the band
            current_block.Initialize(current_cell);
            phi=levelset.Phi(current_block,X);inside=(phi>=phi_min && phi<=phi_max);
            if(!inside){
                dt*=.5;iteration++;X=cell_particles.X(index)+dt*distance*N;
                current_cell=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,X);
                if(!current_cell) break; // out of the band
                current_block.Initialize(current_cell);
                phi=levelset.Phi(current_block,X);inside=(phi>=phi_min && phi<=phi_max);}
            cell_particles.X(index)=X;}}
    if(!inside){Delete_Particle(cell_particles,index);return false;}
    else{
        assert(levelset.collision_body_list);
        // TODO: the old collision aware code did not do this callback for fer that it would hamper collisions
        TV V_temp;levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,index,V_temp,particle_type,0,time);
        phi=levelset.Collision_Aware_Phi(cell_particles.X(index));
        if(phi<phi_min || phi>phi_max){Delete_Particle(cell_particles,index);return false;}
        cell_particles.radius(index)=clamp(abs(phi),minimum_particle_radius,maximum_particle_radius);
        if(delete_particles_that_leave_original_cell && block.Base_Cell()!=current_cell){Delete_Particle(cell_particles,index);return false;}
        if(!particles(current_cell->Cell())) particles(current_cell->Cell())=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
        if(block.Base_Cell()!=current_cell) Move_Particle(cell_particles,*particles(current_cell->Cell()),index);
        return true;}
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Adjust_Particle_Radii()
{
    ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    for(int cell_index=1;cell_index<=negative_particles.m;cell_index++) if(negative_particles(cell_index))
        Adjust_Particle_Radii(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*negative_particles(cell_index),-1);
    for(int cell_index=1;cell_index<=positive_particles.m;cell_index++) if(positive_particles(cell_index))
        Adjust_Particle_Radii(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*positive_particles(cell_index),1);
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Adjust_Particle_Radii(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign)
{
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    // new radius is negative if the particle is on the wrong side of the interface
    if(sign==1) for(int k=1;k<=particles.array_collection->Size();k++) particles.radius(k)=max(minimum_particle_radius,min(maximum_particle_radius,levelset.Phi(block,particles.X(k))));
    else for(int k=1;k<=particles.array_collection->Size();k++) particles.radius(k)=max(minimum_particle_radius,min(maximum_particle_radius,-levelset.Phi(block,particles.X(k))));
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Modify_Levelset_Using_Escaped_Particles()
{
    if(bias_towards_negative_particles){
        Modify_Levelset_Using_Escaped_Particles(levelset.phi,positive_particles,1);Modify_Levelset_Using_Escaped_Particles(levelset.phi,negative_particles,-1);
        if(reincorporate_removed_particles_everywhere) Modify_Levelset_Using_Escaped_Particles(levelset.phi,removed_negative_particles,-1);}
    else{
        ARRAY<T> phi_minus(levelset.phi),phi_plus(levelset.phi);
        Modify_Levelset_Using_Escaped_Particles(phi_minus,negative_particles,-1);Modify_Levelset_Using_Escaped_Particles(phi_plus,positive_particles,1);
        for(int i=1;i<=levelset.grid.number_of_cells;i++) if(abs(levelset.phi(i))<=half_band_width) levelset.phi(i)=minmag(phi_minus(i),phi_plus(i));}
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Modify_Levelset_Using_Escaped_Particles(ARRAY<T>& phi,ARRAY<T_PARTICLES*>& particles,const int sign)
{
    T one_over_radius_multiplier=-sign/outside_particle_distance_multiplier;
    for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_PARTICLES<TV>* particles_link=particles(iterator.Cell_Index());BLOCK_DYADIC<T_GRID> block(levelset.grid,iterator.Cell_Pointer());
        bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
        while(particles_link){
            PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles_link;
            for(int k=1;k<=cell_particles.array_collection->Size();k++) if(one_over_radius_multiplier*levelset.Phi(block,cell_particles.X(k))>cell_particles.radius(k)){
                for(int ii=0;ii<T_GRID::number_of_cells_per_node;ii++){
                    T radius_minus_sign_phi=cell_particles.radius(k)-sign*phi(block.cells[ii]->Cell());
                    if(radius_minus_sign_phi>0){TV center=block.cells[ii]->Center();
                        T distance_squared=(center-cell_particles.X(k)).Magnitude_Squared();
                        if(distance_squared<sqr(radius_minus_sign_phi)){
                            static COLLISION_GEOMETRY_ID body_id;static int aggregate_id;static TV intersection_point;
                            if(near_objects && levelset.collision_body_list->collision_geometry_collection.Intersection_Between_Points(center,cell_particles.X(k),body_id,aggregate_id,intersection_point)) continue;
                            phi(block.cells[ii]->Cell())=sign*(cell_particles.radius(k)-sqrt(distance_squared));}}}}
            particles_link=particles_link->next;}
        if(near_objects) levelset.Disable_Collision_Aware_Interpolation();}
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Euler_Step_Particles(const ARRAY<T>& face_velocities,const T dt,const T time,const bool use_second_order_for_nonremoved_particles,const bool update_particle_cells_after_euler_step,
    const bool verbose)
{
    if(use_second_order_for_nonremoved_particles){
        Second_Order_Runge_Kutta_Step_Particles(face_velocities,negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,time,update_particle_cells_after_euler_step,verbose);
        Second_Order_Runge_Kutta_Step_Particles(face_velocities,positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,time,update_particle_cells_after_euler_step,verbose);}
    else{
        Euler_Step_Particles(face_velocities,negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,time,update_particle_cells_after_euler_step);
        Euler_Step_Particles(face_velocities,positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,time,update_particle_cells_after_euler_step);}

    Euler_Step_Removed_Particles(dt,time,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Euler_Step_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    if(use_removed_negative_particles)
        Euler_Step_Removed_Particles(removed_negative_particles,PARTICLE_LEVELSET_REMOVED_NEGATIVE,dt,time,update_particle_cells_after_euler_step,verbose);
    if(use_removed_positive_particles)
        Euler_Step_Removed_Particles(removed_positive_particles,PARTICLE_LEVELSET_REMOVED_POSITIVE,dt,time,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Euler_Step_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Euler_Step_Removed_Particles(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
    const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    int number_of_deleted_particles=0,number_of_non_occupied_cells=0;
    for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
        if(!levelset.collision_body_list->Swept_Occupied_Block(BLOCK_DYADIC<T_GRID>(levelset.grid,iterator.Cell_Pointer()))){
            number_of_non_occupied_cells++;
            for(int k=cell_particles.array_collection->Size();k>=1;k--) // since not an occupied cell, don't need to adjust for objects
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,cell_particles.V(k),particle_type,dt,time);}
        else{
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,cell_particles.V(k),particle_type,dt,time);
                T collision_distance=Particle_Collision_Distance(cell_particles.quantized_collision_distance(k));
                if(!Adjust_Particle_For_Objects(cell_particles.X(k),cell_particles.V(k),cell_particles.radius(k),collision_distance,particle_type,dt,time)){
                    cell_particles.array_collection->Delete_Element(k);number_of_deleted_particles++;}}}
        cell_particles.Euler_Step_Position(dt);}
    if(verbose){
        if(number_of_deleted_particles) LOG::cout<<"Deleted "<<number_of_deleted_particles<<" "<<PARTICLE_LEVELSET<T_GRID>::Particle_Type_Name(particle_type)<<" due to crossover"<<std::endl;
        if(number_of_non_occupied_cells) LOG::cout<<"Skipped "<<number_of_non_occupied_cells<<" non occupied cells"<<std::endl;}
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Euler_Step_Particles(const ARRAY<T>& face_velocities,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
    const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool assume_particles_in_correct_cells)
{
    LINEAR_INTERPOLATION_DYADIC<T_GRID,TV> linear_interpolation;
    if(assume_particles_in_correct_cells)
        for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
            PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
            for(int k=1;k<=cell_particles.array_collection->Size();k++){
                TV velocity=linear_interpolation.From_Block_Face(levelset.grid,BLOCK_DYADIC<T_GRID>(levelset.grid,iterator.Cell_Pointer()),face_velocities,cell_particles.X(k));
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
                T collision_distance=Particle_Collision_Distance(cell_particles.quantized_collision_distance(k));
                Adjust_Particle_For_Objects(cell_particles.X(k),velocity,cell_particles.radius(k),collision_distance,particle_type,dt,time);
                cell_particles.X(k)+=dt*velocity;}}
    else
        for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
            PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
            for(int k=1;k<=cell_particles.array_collection->Size();k++){
                TV velocity=linear_interpolation.Clamped_To_Array_Face(levelset.grid,face_velocities,0,cell_particles.X(k));
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
                T collision_distance=Particle_Collision_Distance(cell_particles.quantized_collision_distance(k));
                Adjust_Particle_For_Objects(cell_particles.X(k),velocity,cell_particles.radius(k),collision_distance,particle_type,dt,time);
                cell_particles.X(k)+=dt*velocity;}}
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Second_Order_Runge_Kutta_Step_Particles(const ARRAY<T>& face_velocities,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,
    const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    T_FACE_LOOKUP face_velocities_lookup(face_velocities);
    T_FACE_LOOKUP_COLLIDABLE face_velocities_lookup_collidable(face_velocities_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP face_velocities_lookup_collidable_lookup(face_velocities_lookup_collidable,face_velocities_lookup);
    T_LINEAR_INTERPOLATION_SCALAR linear_interpolation; // use for second step since particle may not be in initial block
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    COLLISION_GEOMETRY_ID body_id;int aggregate_id;T start_phi,end_phi;TV body_normal,body_velocity;int number_of_deleted_particles=0,number_of_non_occupied_cells=0;
    for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
        CELL* current_cell=iterator.Cell_Pointer();BLOCK_DYADIC<T_GRID> current_block(levelset.grid,current_cell);
        if(!levelset.collision_body_list->Swept_Occupied_Block(current_block)){ // not occupied, so advect ignoring objects
            number_of_non_occupied_cells++;
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                TV velocity=linear_interpolation.From_Block_Face(levelset.grid,current_block,face_velocities,cell_particles.X(k));
                TV X_new=cell_particles.X(k)+dt*velocity;
                CELL* cell_new=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,X_new);
                if(!cell_new){cell_particles.array_collection->Delete_Element(k);continue;}
                velocity=(T).5*(velocity+linear_interpolation.From_Block_Face(levelset.grid,BLOCK_DYADIC<T_GRID>(levelset.grid,cell_new),face_velocities,X_new));
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
                cell_particles.X(k)+=dt*velocity;}}
        else{ // collision aware advection
            assert(!cell_particles.array_collection->Size() || !current_cell->Has_Children()); // ensure no particles in interior cells
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                T collision_distance=Particle_Collision_Distance(cell_particles.quantized_collision_distance(k));
                BLOCK_DYADIC<T_GRID> block_new=current_block;
                // first push particle out
                if(particle_type==PARTICLE_LEVELSET_NEGATIVE){bool particle_crossover;
                    if(levelset.collision_body_list->Push_Out_Point(cell_particles.X(k),collision_distance,true,particle_crossover)&&particle_crossover){
                        cell_particles.array_collection->Delete_Element(k);number_of_deleted_particles++;continue;} // delete due to crossover in push out
                    CELL* cell_new=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,cell_particles.X(k));
                    if(!cell_new){cell_particles.array_collection->Delete_Element(k);continue;}
                    block_new.Initialize(cell_new);}
                TV velocity=linear_interpolation_collidable.From_Block_Face(levelset.grid,block_new,face_velocities_lookup_collidable_lookup,cell_particles.X(k));
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
                CELL* cell_new=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,X_new);
                if(!cell_new){cell_particles.array_collection->Delete_Element(k);continue;}
                block_new.Initialize(cell_new);
                if(!cell_new){cell_particles.array_collection->Add_To_Deletion_List(k);continue;}
                velocity=(T).5*(velocity+linear_interpolation_collidable.From_Block_Face(levelset.grid,block_new,face_velocities_lookup_collidable_lookup,X_new));
                // adjust normal component of velocity again
                if(got_interaction_with_body){
                    T relative_normal_velocity=TV::Dot_Product(velocity-body_velocity,body_normal);
                    if(relative_normal_velocity<0) velocity-=relative_normal_velocity*body_normal;}
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
                if(!Adjust_Particle_For_Objects(cell_particles.X(k),velocity,cell_particles.radius(k),collision_distance,particle_type,dt,time)){
                    cell_particles.array_collection->Delete_Element(k);number_of_deleted_particles++;}
                else cell_particles.X(k)+=dt*velocity;}}}
    if(verbose){
        if(number_of_deleted_particles) LOG::cout<<"Deleted "<<number_of_deleted_particles<<" "<<PARTICLE_LEVELSET<T_GRID>::Particle_Type_Name(particle_type)<<" due to crossover"<<std::endl;
        if(number_of_non_occupied_cells) LOG::cout<<"Skipped "<<number_of_non_occupied_cells<<" non occupied cells"<<std::endl;}
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_DYADIC<T_GRID>::
Reseed_Particles(const T time,ARRAY<bool>* cell_centered_mask)
{
    int new_particles=0;
    bool normals_defined=(levelset.normals!=0);levelset.Compute_Normals(); // make sure normals are accurate

    if(!cell_centered_mask) new_particles-=Reseed_Delete_Particles(negative_particles,-1);
    new_particles+=Reseed_Add_Particles(negative_particles,positive_particles,-1,PARTICLE_LEVELSET_NEGATIVE,time,cell_centered_mask);

    if(!cell_centered_mask) new_particles-=Reseed_Delete_Particles(positive_particles,1);
    new_particles+=Reseed_Add_Particles(positive_particles,negative_particles,1,PARTICLE_LEVELSET_POSITIVE,time,cell_centered_mask);

    if(!normals_defined){delete levelset.normals;levelset.normals=0;}
    return new_particles;
}
//#####################################################################
// Function Reseed_Delete_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_DYADIC<T_GRID>::
Reseed_Delete_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const int sign)
{
    int number_deleted=0;ARRAY<int> heap_particle_indices(number_particles_per_cell);ARRAY<T> heap_phi_minus_radius(number_particles_per_cell);
    for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        int cell_index=iterator.Cell_Index();CELL* base_cell=iterator.Cell_Pointer();BLOCK_DYADIC<T_GRID> base_block(levelset.grid,base_cell);
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(cell_index);
        for(int i=0;i<T_GRID::number_of_cells_per_block;i++){
            T unsigned_phi=sign*levelset.phi(base_block.cells[i]->Cell());if(0<=unsigned_phi && unsigned_phi<=half_band_width) goto NEAR_THE_INTERFACE;}
        
        if(cell_particles.array_collection->Size()==0) Free_Particle_And_Clear_Pointer(particles(cell_index));
        else{ // delete all non-escaped particles
            ARRAY<bool> escaped;Identify_Escaped_Particles(base_block,cell_particles,escaped,sign);
            for(int index=cell_particles.array_collection->Size();index>=1;index--) if(!escaped(index)) cell_particles.array_collection->Delete_Element(index);
            if(cell_particles.array_collection->Size()==0) Free_Particle_And_Clear_Pointer(particles(cell_index));}
        continue;
        
        NEAR_THE_INTERFACE:; // can only get here via the goto
        if(!particles(cell_index) || particles(cell_index)->array_collection->Size()<=number_particles_per_cell) continue;
        ARRAY<bool> escaped;Identify_Escaped_Particles(base_block,*particles(cell_index),escaped,sign);int number_of_escaped_particles=escaped.Number_True();
        int total_particles=particles(cell_index)->array_collection->Size()-number_of_escaped_particles;
        if(total_particles>number_particles_per_cell){ // too many particles - delete particles with a heap sort
            number_deleted+=number_particles_per_cell-total_particles;int heap_size=0;
            for(int index=1;index<=cell_particles.array_collection->Size();index++) if(!escaped(index)){
                T phi_minus_radius=sign*levelset.Phi(base_block,cell_particles.X(index))-cell_particles.radius(index);
                if(heap_size<number_particles_per_cell){ // add particle to heap
                    heap_size++;heap_particle_indices(heap_size)=index;heap_phi_minus_radius(heap_size)=phi_minus_radius;
                    if(heap_size==number_particles_per_cell) ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices);} // when heap is full, order values with largest on top
                else{ // excess particles don't fit in the heap
                    if(phi_minus_radius<heap_phi_minus_radius(1)){ // delete particle on top of heap & add new particle
                        cell_particles.array_collection->Add_To_Deletion_List(heap_particle_indices(1));
                        heap_phi_minus_radius(1)=phi_minus_radius;heap_particle_indices(1)=index;
                        ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices,1,heap_phi_minus_radius.m);}
                    else cell_particles.array_collection->Add_To_Deletion_List(index);}} // delete new particle, larger than top of heap
            cell_particles.array_collection->Delete_Elements_On_Deletion_List();}}
    return number_deleted;
}
//#####################################################################
// Function Reseed_Add_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_DYADIC<T_GRID>::
Reseed_Add_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& other_particles,const int sign,
    const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,ARRAY<bool>* cell_centered_mask)
{
    int number_added=0;
    for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()){
        int cell_index=iterator.Cell_Index();CELL* cell=iterator.Cell_Pointer();
        if(!levelset.grid.fully_refined_block(cell->Cell()))continue;BLOCK_DYADIC<T_GRID> block(levelset.grid,cell);
        if(!levelset.grid.Domain().Inside(block.Center(),(T)-1e-6))continue;
        for(int i=0;i<T_GRID::number_of_cells_per_block;i++){
            T unsigned_phi=sign*levelset.phi(block.cells[i]->Cell());
            if(0<=unsigned_phi && unsigned_phi<=half_band_width && (!cell_centered_mask||(*cell_centered_mask)(block.cells[i]->Cell()))) goto NEAR_THE_INTERFACE;}
        continue;
        
        NEAR_THE_INTERFACE:; // can only get here via the goto
        if(!particles(cell_index)) particles(cell_index)=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
        int total_particles=particles(cell_index)->array_collection->Size(),total_other_particles=0;if(other_particles(cell_index)) total_other_particles=other_particles(cell_index)->array_collection->Size();
        if(total_particles+total_other_particles>=number_particles_per_cell) continue;
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(cell_index);
        int number_of_particles_to_add=number_particles_per_cell-total_particles-total_other_particles;
        if(sign==-1){ // we add the negative particles first, and don't want to add too many of them...
            if(total_other_particles) number_of_particles_to_add=(int)((T)number_of_particles_to_add*(T)total_particles/(T)(total_particles+total_other_particles)+1);
            else if(!total_particles) number_of_particles_to_add=(number_of_particles_to_add+1)/2;}
        number_added+=number_of_particles_to_add;
        T phi_min=sign*minimum_particle_radius,phi_max=sign*half_band_width;if(phi_min>phi_max) exchange(phi_min,phi_max);
        RANGE<TV> block_bounding_box=block.Bounding_Box();
        int attempts=0;
        ARRAY_VIEW<int>* id=store_unique_particle_id?cell_particles.array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID):0;
        for(int k=1;k<=number_of_particles_to_add;k++){
            int index=cell_particles.array_collection->Add_Element();
            if(id) (*id)(index)=++last_unique_particle_id;
            cell_particles.quantized_collision_distance(index)=(unsigned short)(random.Get_Number()*USHRT_MAX);
            cell_particles.X(index)=random.Get_Uniform_Vector(block_bounding_box);
            if(!Attract_Individual_Particle_To_Interface_And_Adjust_Radius(particles,phi_min,phi_max,block,index,particle_type,time,true) && ++attempts<=5) k--;}}
    return number_added;
}
//#####################################################################
// Function Identify_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Identify_Escaped_Particles(const int sign)
{
    ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    if(!sign || sign==1){
        escaped_positive_particles.Resize(levelset.grid.number_of_cells);
        for(int cell_index=1;cell_index<=positive_particles.m;cell_index++)
            if(positive_particles(cell_index)) 
                Identify_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*positive_particles(cell_index),escaped_positive_particles(cell_index),1);
            else escaped_positive_particles(cell_index).Resize(0);}
    if(!sign || sign==-1){
        escaped_negative_particles.Resize(levelset.grid.number_of_cells);
        for(int cell_index=1;cell_index<=negative_particles.m;cell_index++)
            if(negative_particles(cell_index)) 
                Identify_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*negative_particles(cell_index),escaped_negative_particles(cell_index),-1);
            else escaped_negative_particles(cell_index).Resize(0);}
}
//#####################################################################
// Function Identify_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Identify_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign)
{
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    escaped.Resize(particles.array_collection->Size());ARRAYS_COMPUTATIONS::Fill(escaped,false);T one_over_radius_multiplier=-sign/outside_particle_distance_multiplier;
    for(int k=1;k<=particles.array_collection->Size();k++) if(one_over_radius_multiplier*levelset.Phi(block,particles.X(k))>particles.radius(k)) escaped(k)=true;
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Delete_Deep_Escaped_Particles
//#####################################################################
// delete those farther than radius_fraction too deep
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Delete_Deep_Escaped_Particles(const T radius_fraction,const bool need_to_identify_escaped_particles,const bool verbose)
{
    ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    if(verbose) LOG::cout<<"positive particles - ";int number_of_positive_particles_deleted=0;
    for(int cell_index=1;cell_index<=positive_particles.m;cell_index++) if(positive_particles(cell_index))
        number_of_positive_particles_deleted+=Delete_Deep_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*positive_particles(cell_index),
                                                            escaped_positive_particles(cell_index),1,radius_fraction,need_to_identify_escaped_particles);
    if(verbose) LOG::cout<<"deleted "<<number_of_positive_particles_deleted<<" particles"<<std::endl;
    if(verbose) LOG::cout<<"negative particles - ";int number_of_negative_particles_deleted=0;
    for(int cell_index=1;cell_index<=negative_particles.m;cell_index++) if(negative_particles(cell_index))
        number_of_negative_particles_deleted+=Delete_Deep_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*negative_particles(cell_index),
                                                            escaped_negative_particles(cell_index),-1,radius_fraction,need_to_identify_escaped_particles);
    if(verbose) LOG::cout<<"deleted "<<number_of_negative_particles_deleted<<" particles"<<std::endl;
}
//#####################################################################
// Function Delete_Deep_Escaped_Particles
//#####################################################################
// delete those farther than radius_fraction too deep
template<class T_GRID> int PARTICLE_LEVELSET_DYADIC<T_GRID>::
Delete_Deep_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign,const T radius_fraction,const bool need_to_identify_escaped_particles)
{
    assert(!particles.next);
    if(need_to_identify_escaped_particles) Identify_Escaped_Particles(block,particles,escaped,sign);
    int deleted=0;int minus_sign=-sign;
    for(int k=particles.array_collection->Size();k>=1;k--) if(escaped(k) && minus_sign*levelset.Phi(block,particles.X(k))>radius_fraction*particles.radius(k)){
        particles.array_collection->Delete_Element(k);deleted++;}
    return deleted;
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Delete_Particles_Outside_Grid()
{
    // delete particles far outside grid
    RANGE<TV> enlarged_domain=levelset.grid.Domain();enlarged_domain.Change_Size((T).25*levelset.grid.Minimum_Edge_Length());
    for(CELL_ITERATOR iterator(levelset.grid,levelset.grid.number_of_ghost_cells,UNIFORM_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){int cell_index=iterator.Cell_Index();
        if(!enlarged_domain.Lazy_Inside(levelset.grid.Node_Location(T_GRID::number_of_nodes_per_cell-1,iterator.Cell_Pointer()))){
            Free_Particle_And_Clear_Pointer(negative_particles(cell_index));
            Free_Particle_And_Clear_Pointer(positive_particles(cell_index));
            if(use_removed_positive_particles) Free_Particle_And_Clear_Pointer(removed_positive_particles(cell_index));
            if(use_removed_negative_particles) Free_Particle_And_Clear_Pointer(removed_negative_particles(cell_index));}}
    // delete particles barely outside grid
    RANGE<TV> domain=levelset.grid.Domain();TV domain_boundaries[2]={domain.Minimum_Corner(),domain.Maximum_Corner()};
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int side=1;side<=2;side++)
        for(NODE_ITERATOR iterator(levelset.grid,levelset.grid.Map_Individual_Side_Boundary_Nodes(2*(axis-1)+side));iterator.Valid();iterator.Next()){
            CELL* cell=iterator.Deepest_Cell();if(cell->Depth_Of_This_Cell()>levelset.grid.maximum_depth)continue;
            if(negative_particles(cell->Cell()))Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*negative_particles(cell->Cell()),2*side-3);
            if(positive_particles(cell->Cell()))Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*positive_particles(cell->Cell()),2*side-3);
            if(use_removed_negative_particles) if(removed_negative_particles(cell->Cell()))
                Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*removed_negative_particles(cell->Cell()),2*side-3);
            if(use_removed_positive_particles) if(removed_positive_particles(cell->Cell()))
                Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*removed_positive_particles(cell->Cell()),2*side-3);}
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int side)
{
    for(int k=particles.array_collection->Size();k>=1;k--) if(side*particles.X(k)[axis]>side*domain_boundary) particles.array_collection->Delete_Element(k); // side should be -1 or 1
}
//#####################################################################
// Function Identify_And_Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Identify_And_Remove_Escaped_Particles(const ARRAY<T>& face_velocities,const T radius_fraction,const bool verbose)
{
    T_FACE_LOOKUP face_velocities_lookup(face_velocities);
    T_FACE_LOOKUP_COLLIDABLE face_velocities_lookup_collidable(face_velocities_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP face_velocities_lookup_collidable_lookup(face_velocities_lookup_collidable,face_velocities_lookup);

    Identify_Escaped_Particles();int p=0,n=0;
    ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    if(use_removed_positive_particles){
        for(int cell_index=1;cell_index<=positive_particles.m;cell_index++) if(positive_particles(cell_index))
            p+=Remove_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*positive_particles(cell_index),escaped_positive_particles(cell_index),1,
                removed_positive_particles(cell_index),face_velocities_lookup_collidable_lookup,radius_fraction);}
    else{
        for(int cell_index=1;cell_index<=positive_particles.m;cell_index++) if(positive_particles(cell_index))
            p+=Remove_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*positive_particles(cell_index),escaped_positive_particles(cell_index),1,
                radius_fraction);}
    if(use_removed_negative_particles){
        for(int cell_index=1;cell_index<=negative_particles.m;cell_index++) if(negative_particles(cell_index))
            n+=Remove_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*negative_particles(cell_index),escaped_negative_particles(cell_index),-1,
                removed_negative_particles(cell_index),face_velocities_lookup_collidable_lookup,radius_fraction);}
    else{
        for(int cell_index=1;cell_index<=negative_particles.m;cell_index++) if(negative_particles(cell_index))
            n+=Remove_Escaped_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),*negative_particles(cell_index),escaped_negative_particles(cell_index),-1,
                radius_fraction);}
    LOG::cout<<"new removed positive particles = "<<p<<std::endl;LOG::cout<<"new removed negative particles = "<<n<<std::endl;
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> template<class T_FACE_LOOKUP_LOOKUP> int PARTICLE_LEVELSET_DYADIC<T_GRID>::
Remove_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles,const T_FACE_LOOKUP_LOOKUP& face_velocities,const T radius_fraction)
{
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    int deleted=0;
    for(int k=particles.array_collection->Size();k>=1;k--) if(escaped(k) && -sign*levelset.Phi(block,particles.X(k))>radius_fraction*particles.radius(k)){
        if(!removed_particles) removed_particles=Allocate_Particle<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >();
        removed_particles->array_collection->Take(*particles.array_collection,k);
        removed_particles->V.Last()=linear_interpolation_collidable.From_Block_Face(levelset.grid,block,face_velocities,removed_particles->X.Last());}
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
    return deleted;
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_DYADIC<T_GRID>::
Remove_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,const T radius_fraction)
{
    assert(!particles.next);
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    int deleted=0;
    for(int k=particles.array_collection->Size();k>=1;k--) if(escaped(k) && -sign*levelset.Phi(block,particles.X(k))>radius_fraction*particles.radius(k)){
        particles.array_collection->Delete_Element(k);deleted++;}
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
    return deleted;
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Reincorporate_Removed_Particles(const T radius_fraction)
{
    ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    if(use_removed_positive_particles)
        for(int cell_index=1;cell_index<=removed_positive_particles.m;cell_index++) if(removed_positive_particles(cell_index))
            Reincorporate_Removed_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),positive_particles(cell_index),1,*removed_positive_particles(cell_index),
                radius_fraction);
    if(use_removed_negative_particles)
        for(int cell_index=1;cell_index<=removed_negative_particles.m;cell_index++) if(removed_negative_particles(cell_index))
            Reincorporate_Removed_Particles(BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(cell_index)),negative_particles(cell_index),-1,*removed_negative_particles(cell_index),
                radius_fraction);
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Reincorporate_Removed_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>*& particles,const int sign,
                                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& removed_particles,const T radius_fraction)
{
    bool near_objects=levelset.collision_body_list->Occupied_Block(block);if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    for(int k=removed_particles.array_collection->Size();k>=1;k--) if(-sign*levelset.Phi(block,removed_particles.X(k))<radius_fraction*removed_particles.radius(k)){
        if(!particles) particles=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
        Move_Particle(removed_particles,*particles,k);} // maybe recalculate radius to get rid of raining effect?
    if(particles) Compact_Particles_Into_Single_Particle_Bin(*particles,block,sign);
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Delete_Particles_Far_From_Interface
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Delete_Particles_Far_From_Interface(const int discrete_band)
{
    assert(discrete_band<8);
    const ARRAY<VECTOR<CELL*,T_GRID::number_of_neighbors_per_cell> >& cell_neighbors=levelset.grid.Neighbors();
    const ARRAY<VECTOR<bool,T_GRID::dimension> >& cell_neighbors_visible=levelset.collision_body_list->cell_neighbors_visible;
    ARRAY<char> near_interface(levelset.grid.number_of_cells);
    for(int cell=1;cell<=levelset.grid.number_of_cells;cell++){
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            CELL* cell_neighbor=cell_neighbors(cell)(2*axis);if(cell_neighbor==0)continue;
            int neighbor=cell_neighbor->Cell();
            if(neighbor && cell_neighbors_visible(cell)(axis) && LEVELSET_UTILITIES<T>::Interface(levelset.phi(cell),levelset.phi(neighbor))){
                near_interface(cell)=near_interface(neighbor)=1;}}}
    for(int distance=1;distance<=discrete_band;distance++){
        char new_mask=1<<distance,old_mask=new_mask-1;
        for(int cell=1;cell<=levelset.grid.number_of_cells;cell++)
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                CELL* cell_neighbor=cell_neighbors(cell)(2*axis);if(cell_neighbor==0)continue;
                int neighbor=cell_neighbor->Cell();
                if(neighbor && cell_neighbors_visible(cell)(axis) && (near_interface(cell)|near_interface(neighbor))&old_mask){
                    near_interface(cell)|=new_mask;near_interface(neighbor)|=new_mask;}}}
    for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()){int b=iterator.Cell_Index();if(negative_particles(b) || positive_particles(b)){
        BLOCK_DYADIC<T_GRID> block(levelset.grid,iterator.Cell_Pointer());
        for(int i=0;i<T_GRID::number_of_cells_per_block;i++) if(near_interface(block.Cell(i))) goto NEAR_INTERFACE;
        // if not near interface, delete particles 
        Free_Particle_And_Clear_Pointer(negative_particles(b));
        Free_Particle_And_Clear_Pointer(positive_particles(b));
        NEAR_INTERFACE:;}}
}
//#####################################################################
// Function Add_Negative_Particle
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_DYADIC<T_GRID>::
Add_Negative_Particle(const TV& location,const T collision_distance_percentage)
{
    PHYSBAM_NOT_IMPLEMENTED(); // TODO: figure out whether/how this should be collision aware

    const CELL* base_cell=levelset.grid.Base_Cell(location);
    if(!base_cell) return;
    int i=base_cell->Cell();
    if(!negative_particles(i))negative_particles(i)=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
    else if(negative_particles(i)->array_collection->Size()>=number_particles_per_cell)return;
    int index=negative_particles(i)->array_collection->Add_Element();
    negative_particles(i)->X(index)=location;negative_particles(i)->radius(index)=maximum_particle_radius;
    if(store_unique_particle_id){
        ARRAY_VIEW<int>* id=negative_particles(i)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);
        PHYSBAM_ASSERT(id);(*id)(index)=++last_unique_particle_id;}
    negative_particles(i)->quantized_collision_distance(index)=(unsigned short)(random.Get_Number()*USHRT_MAX*collision_distance_percentage);
}
//#####################################################################
template class PARTICLE_LEVELSET_DYADIC<QUADTREE_GRID<float> >;
template class PARTICLE_LEVELSET_DYADIC<OCTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_DYADIC<QUADTREE_GRID<double> >;
template class PARTICLE_LEVELSET_DYADIC<OCTREE_GRID<double> >;
#endif
#endif
