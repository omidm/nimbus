//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Duc Nguyen, Nick Rasmussen, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/VOF_ADVECTION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#ifndef WIN32
    #include <ext/atomicity.h>
#endif
using namespace PhysBAM;
#ifndef WIN32
    using namespace __gnu_cxx;
#endif
//#####################################################################
// Constructor
//##################################################################### 
template<class T_GRID> PARTICLE_LEVELSET_UNIFORM<T_GRID>::
PARTICLE_LEVELSET_UNIFORM(T_GRID& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells_input)
    :PARTICLE_LEVELSET<T_GRID>(grid_input,phi_input,number_of_ghost_cells_input),mpi_grid(0),vof_advection(0),thread_queue(0)
{
#ifdef USE_PTHREADS
    pthread_mutex_init(&cell_lock,0);
#endif
}
//#####################################################################
// Destructor
//##################################################################### 
template<class T_GRID> PARTICLE_LEVELSET_UNIFORM<T_GRID>::
~PARTICLE_LEVELSET_UNIFORM()
{
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&cell_lock);
#endif
    delete vof_advection;
    for(NODE_ITERATOR iterator(levelset.grid,positive_particles.Domain_Indices());iterator.Valid();iterator.Next()) Delete_All_Particles_In_Cell(iterator.Node_Index());
}
//#####################################################################
// Function Consistency_Check
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Consistency_Check(RANGE<TV_INT> domain,T_ARRAYS_PARTICLES& particles)
{
  // Remove consistency check for testing purposes.
  /*
    HASHTABLE<typename T_ARRAYS_PARTICLES::ELEMENT> hash;
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block);
        while(cell_particles){
            assert(!hash.Contains(cell_particles));hash.Insert(cell_particles);
            assert(cell_particles->array_collection->Size()!=0);
            if(cell_particles->next)assert(cell_particles->array_collection->Size()==particle_pool.number_particles_per_cell);
            for(int k=1;k<=cell_particles->array_collection->Size();k++)
                assert(cell_particles->radius(k)>0);
            cell_particles=cell_particles->next;}}
            */
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Seed_Particles(const T time,const bool verbose)
{
    int n=Reseed_Add_Particles(negative_particles,positive_particles,-1,PARTICLE_LEVELSET_NEGATIVE,time,0);

    int p=0;
    if(!only_use_negative_particles){
        p=Reseed_Add_Particles(positive_particles,negative_particles,1,PARTICLE_LEVELSET_POSITIVE,time,0);}

    if(verbose) {std::stringstream ss;ss<<n<<" negative particles and "<<p<<" positive particles"<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Attract_Individual_Particle_To_Interface
//#####################################################################
template<class T_GRID> bool PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Attract_Individual_Particle_To_Interface_And_Adjust_Radius(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles,const T phi_min,const T phi_max,const BLOCK_UNIFORM<T_GRID>& block,
    const int index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const bool delete_particles_that_leave_original_cell,RANDOM_NUMBERS<T>& local_random,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index)
{
    T phi=levelset.Phi(cell_particles.X(index));bool inside=(phi>=phi_min && phi<=phi_max);
    TV_INT current_block_index=block.block_index;
    if(!inside){
        T phi_goal=local_random.Get_Uniform_Number(phi_min,phi_max);TV X,N;int iteration=0;
        while(!inside && iteration<maximum_iterations_for_attraction){
            N=levelset.Normal(cell_particles.X(index));T distance=phi_goal-phi,dt=1;bool inside_domain=false;
            while(!inside_domain && iteration<maximum_iterations_for_attraction){
                X=cell_particles.X(index)+dt*distance*N;
                if(levelset.grid.Outside(X)){dt*=.5;iteration++;}else inside_domain=true;}
            if(!inside_domain) break; // ran out of iterations
            phi=levelset.Phi(X);inside=(phi>=phi_min && phi<=phi_max);
            if(!inside){
                dt*=.5;iteration++;X=cell_particles.X(index)+dt*distance*N;
                phi=levelset.Phi(X);inside=(phi>=phi_min && phi<=phi_max);}
            cell_particles.X(index)=X;}}
    if(!inside){
        assert(cell_particles.radius(index)>0);
        cell_particles.radius(index)=-1;
        return false;}
    else{
        // TODO: the old collision aware code did not do this callback for fear that it would hamper collisions
        TV V_temp;levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,index,V_temp,particle_type,0,time);
        if(levelset.collision_body_list) phi=levelset.Collision_Aware_Phi(cell_particles.X(index));
        else phi=levelset.Phi(cell_particles.X(index));
        if(phi<phi_min || phi>phi_max){
            assert(cell_particles.radius(index)>0);
            cell_particles.radius(index)=-1;
            return false;}
        cell_particles.radius(index)=clamp(abs(phi),minimum_particle_radius,maximum_particle_radius);
        current_block_index=levelset.grid.Block_Index(cell_particles.X(index),1);
        if(delete_particles_that_leave_original_cell && block.block_index!=current_block_index){
            assert(cell_particles.radius(index)>0);
            cell_particles.radius(index)=-1;
            return false;}
        //if(!particles(current_block_index)) particles(current_block_index)=Allocate_Particles(template_particles);
        //if(block.block_index!=current_block_index) Move_Particle(cell_particles,*particles(current_block_index),index);
        //move_particles(block.block_index).Append(TRIPLE<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT,int>(&cell_particles,current_block_index,index));
        list_to_process(current_block_index).Append(TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int>(block.block_index,&cell_particles,index));
        /*if(mutexes && thread_queue) pthread_mutex_lock(&(*mutexes)((*domain_index)(current_block_index)));
        if(!particles(current_block_index)) particles(current_block_index)=Allocate_Particles(template_particles);
        Copy_Particle(cell_particles,*particles(current_block_index),index);
        if(mutexes && thread_queue) pthread_mutex_unlock(&(*mutexes)((*domain_index)(current_block_index)));
        assert(cell_particles.radius(index)>0);
        cell_particles.radius(index)=-1;*/
        return true;}
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Adjust_Particle_Radii()
{
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(domain,thread_queue).Run(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Adjust_Particle_Radii_Threaded);
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Adjust_Particle_Radii_Threaded(RANGE<TV_INT>& domain)
{
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(negative_particles(block)) Adjust_Particle_Radii(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*negative_particles(block),-1);}
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(positive_particles(block)) Adjust_Particle_Radii(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*positive_particles(block),1);}
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Adjust_Particle_Radii(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign)
{
    bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block):false;if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    // new radius is negative if the particle is on the wrong side of the interface
    PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=&particles;
    while(cell_particles){
        if(sign==1) for(int k=1;k<=cell_particles->array_collection->Size();k++)cell_particles->radius(k)=max(minimum_particle_radius,min(maximum_particle_radius,levelset.Phi(particles.X(k))));
        else for(int k=1;k<=cell_particles->array_collection->Size();k++)cell_particles->radius(k)=max(minimum_particle_radius,min(maximum_particle_radius,-levelset.Phi(particles.X(k))));
        cell_particles=cell_particles->next;}
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Modify_Levelset_Using_Escaped_Particles(T_FACE_ARRAYS_SCALAR* V,ARRAY<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES*>* other_positive_particles)
{
    if(bias_towards_negative_particles){
        if(other_positive_particles)for(int i=1;i<=other_positive_particles->m;i++)Modify_Levelset_Using_Escaped_Particles(levelset.phi,*(*other_positive_particles)(i),0,1);
        Modify_Levelset_Using_Escaped_Particles(levelset.phi,positive_particles,0,1);Modify_Levelset_Using_Escaped_Particles(levelset.phi,negative_particles,0,-1);
        if(reincorporate_removed_particles_everywhere) Modify_Levelset_Using_Escaped_Particles(levelset.phi,removed_negative_particles,V,-1);}
    else{
        T_ARRAYS_SCALAR phi_minus(levelset.phi),phi_plus(levelset.phi);
        Modify_Levelset_Using_Escaped_Particles(phi_minus,negative_particles,0,-1);Modify_Levelset_Using_Escaped_Particles(phi_plus,positive_particles,0,1);
        if(other_positive_particles)for(int i=1;i<=other_positive_particles->m;i++)Modify_Levelset_Using_Escaped_Particles(phi_plus,*(*other_positive_particles)(i),0,1);
        if(reincorporate_removed_particles_everywhere) Modify_Levelset_Using_Escaped_Particles(phi_minus,removed_negative_particles,V,-1);
        for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();levelset.phi(cell)=minmag(phi_minus(cell),phi_plus(cell));}}
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Modify_Levelset_Using_Escaped_Particles(T_ARRAYS_SCALAR& phi,T_ARRAYS_PARTICLES& particles,T_FACE_ARRAYS_SCALAR* V,const int sign)
{
    T one_over_radius_multiplier=-sign/outside_particle_distance_multiplier;
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(domain,particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(levelset.grid.Domain_Indices(),thread_queue).template Run<T_ARRAYS_SCALAR&,T_ARRAYS_PARTICLES&,T_FACE_ARRAYS_SCALAR*,int,T>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Modify_Levelset_Using_Escaped_Particles_Threaded,phi,particles,V,sign,one_over_radius_multiplier);
    Consistency_Check(domain,particles);
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Modify_Levelset_Using_Escaped_Particles_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi,T_ARRAYS_PARTICLES& particles,T_FACE_ARRAYS_SCALAR* V,const int sign,const T one_over_radius_multiplier)
{
    RANGE<TV_INT> range_domain(domain);range_domain.max_corner+=TV_INT::All_Ones_Vector();
    for(NODE_ITERATOR iterator(levelset.grid,range_domain);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block_index);
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block):false;if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
        while(cell_particles){for(int k=1;k<=cell_particles->array_collection->Size();k++)if(one_over_radius_multiplier*levelset.Phi(cell_particles->X(k))>cell_particles->radius(k)){
            for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++){TV_INT cell=block.Cell(cell_index);if(!domain.Lazy_Inside(cell)) continue;
                T radius_minus_sign_phi=cell_particles->radius(k)-sign*phi(cell);
                if(radius_minus_sign_phi>0){TV center=levelset.grid.Center(cell);
                    T distance_squared=(center-cell_particles->X(k)).Magnitude_Squared();
                    if(distance_squared<sqr(radius_minus_sign_phi)){
                        static COLLISION_GEOMETRY_ID body_id;static int aggregate_id;static TV intersection_point;
                        if(near_objects && levelset.collision_body_list->collision_geometry_collection.Intersection_Between_Points(center,cell_particles->X(k),body_id,aggregate_id,intersection_point)) continue;
                        phi(cell)=sign*(cell_particles->radius(k)-sqrt(distance_squared));
                        /*if(V) (*V)(ii,jj,iijj)=cell_particles.V(k);*/}}}}
            cell_particles=cell_particles->next;}
        if(near_objects) levelset.Disable_Collision_Aware_Interpolation();}}
}
//#####################################################################
// Function Update_Particles_To_Reflect_Mass_Conservation
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Update_Particles_To_Reflect_Mass_Conservation(T_ARRAYS_SCALAR& phi_old,const bool update_particle_cells,const bool verbose)
{
    T_FAST_LEVELSET levelset_old(levelset.grid,phi_old);
    Update_Particles_To_Reflect_Mass_Conservation(levelset_old,negative_particles,PARTICLE_LEVELSET_NEGATIVE,update_particle_cells,verbose);
    Update_Particles_To_Reflect_Mass_Conservation(levelset_old,positive_particles,PARTICLE_LEVELSET_POSITIVE,update_particle_cells,verbose);
    Adjust_Particle_Radii();
}
//#####################################################################
// Function Update_Particles_To_Reflect_Mass_Conservation
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Update_Particles_To_Reflect_Mass_Conservation(const T_FAST_LEVELSET& levelset_old,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,
    const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const bool update_particle_cells,const bool verbose)
{
    const int maximum_iterations=5;
    const T epsilon=levelset_old.grid.Minimum_Edge_Length()*(T)1e-2;
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=particles(block_index);
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        while(cell_particles){for(int k=cell_particles->array_collection->Size();k>=1;k--){
            TV X_new=cell_particles->X(k);
            T phi_old=levelset_old.Phi(X_new);
            T phi_new=levelset.Phi(X_new);
            int iterations=0;
            while(abs(phi_old-phi_new)>epsilon && iterations<maximum_iterations){
                // figure out how much phi has changed at this point, step that far in the normal direction
                TV normal=levelset.Normal(cell_particles->X(k));
                TV destination_normal=levelset.Normal(cell_particles->X(k)+(phi_old-phi_new)*normal);
                X_new+=(T).5*(phi_old-phi_new)*(normal+destination_normal).Normalized(); // half
                TV velocity=X_new-cell_particles->X(k);
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(*cell_particles,k,velocity,particle_type,1/*=dt*/,0);
                cell_particles->X(k)+=velocity;
                phi_new=levelset.Phi(X_new);
                X_new=cell_particles->X(k);
                iterations++;}}
            cell_particles=cell_particles->next;}}}
    if(update_particle_cells) Update_Particle_Cells(particles);
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Euler_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,const T dt,const T time,const bool use_second_order_for_nonremoved_particles,const bool update_particle_cells_after_euler_step,
    const bool verbose,const bool analytic_test,const bool default_to_forward_euler_outside)
{
    if(use_second_order_for_nonremoved_particles){
        Second_Order_Runge_Kutta_Step_Particles(V,negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,time,update_particle_cells_after_euler_step,verbose,default_to_forward_euler_outside);
        Second_Order_Runge_Kutta_Step_Particles(V,positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,time,update_particle_cells_after_euler_step,verbose,default_to_forward_euler_outside);}
    else{
        Euler_Step_Particles(V,negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,time,update_particle_cells_after_euler_step);
        Euler_Step_Particles(V,positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,time,update_particle_cells_after_euler_step);}
    if(analytic_test){
        Second_Order_Runge_Kutta_Step_Particles(V,removed_negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,time,update_particle_cells_after_euler_step,default_to_forward_euler_outside);
        Second_Order_Runge_Kutta_Step_Particles(V,removed_positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,time,update_particle_cells_after_euler_step,default_to_forward_euler_outside);}
    else
        Euler_Step_Removed_Particles(dt,time,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Euler_Step_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    if(use_removed_negative_particles) Euler_Step_Removed_Particles(removed_negative_particles,PARTICLE_LEVELSET_REMOVED_NEGATIVE,dt,time,update_particle_cells_after_euler_step,verbose);
    if(use_removed_positive_particles) Euler_Step_Removed_Particles(removed_positive_particles,PARTICLE_LEVELSET_REMOVED_POSITIVE,dt,time,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Euler_Step_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Euler_Step_Removed_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
    const bool update_particle_cells_after_euler_step,const bool verbose)
{
    int number_of_deleted_particles=0,number_of_non_occupied_cells=0;
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(particles(block)){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& cell_particles=*particles(block);
        if(!levelset.collision_body_list->Swept_Occupied_Block(T_BLOCK(levelset.grid,block))){
            number_of_non_occupied_cells++;
            for(int k=cell_particles.array_collection->Size();k>=1;k--) // since not an occupied block, don't need to adjust for objects
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,cell_particles.V(k),particle_type,dt,time);}
        else{
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,cell_particles.V(k),particle_type,dt,time);
                T collision_distance=Particle_Collision_Distance(cell_particles.quantized_collision_distance(k));
                if(!Adjust_Particle_For_Objects(cell_particles.X(k),cell_particles.V(k),cell_particles.radius(k),collision_distance,particle_type,dt,time)){
                    Delete_Particle_And_Clean_Memory(particles(block),cell_particles,k);number_of_deleted_particles++;}}}
        cell_particles.Euler_Step_Position(dt);
        if(!cell_particles.array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(block));}}
    if(verbose){
        std::stringstream ss;
        if(number_of_deleted_particles) ss<<"Deleted "<<number_of_deleted_particles<<" "<<PARTICLE_LEVELSET<T_GRID>::Particle_Type_Name(particle_type)<<" due to crossover"<<std::endl;
        if(number_of_non_occupied_cells) ss<<"Skipped "<<number_of_non_occupied_cells<<" non occupied cells"<<std::endl;
        LOG::filecout(ss.str());}
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Euler_Step_Particles_Wrapper(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
    const bool update_particle_cells_after_euler_step,const bool assume_particles_in_correct_blocks,const bool enforce_domain_boundaries)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(levelset.grid,number_of_ghost_cells,false);
    levelset.boundary->Fill_Ghost_Cells_Face(levelset.grid,V,face_velocities_ghost,time,number_of_ghost_cells);
    Euler_Step_Particles(face_velocities_ghost,particles,particle_type,dt,time,update_particle_cells_after_euler_step,assume_particles_in_correct_blocks,enforce_domain_boundaries);
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Euler_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
    const bool update_particle_cells_after_euler_step,const bool assume_particles_in_correct_blocks,const bool enforce_domain_boundaries)
{
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(domain,particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(domain,thread_queue).template Run<const T_FACE_ARRAYS_SCALAR&,T_ARRAYS_PARTICLES&,PARTICLE_LEVELSET_PARTICLE_TYPE,T,T,bool,bool>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Euler_Step_Particles_Threaded,V,particles,particle_type,dt,time,assume_particles_in_correct_blocks,enforce_domain_boundaries);
    Consistency_Check(domain,particles);
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
    Consistency_Check(domain,particles);
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Euler_Step_Particles_Threaded(RANGE<TV_INT>& domain,const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
    const bool assume_particles_in_correct_blocks,const bool enforce_domain_boundaries)
{
    if(assume_particles_in_correct_blocks){
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(particles(block)){
            // TODO: optimize using slope precomputation
            typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block);
            T_LINEAR_INTERPOLATION_MAC_HELPER linear_interpolation_mac_helper(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),V);
            while(cell_particles){for(int k=1;k<=cell_particles->array_collection->Size();k++){
                TV velocity=linear_interpolation_mac_helper.Interpolate_Face(cell_particles->X(k));
                if(enforce_domain_boundaries) levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(*cell_particles,k,velocity,particle_type,dt,time);
                T collision_distance=Particle_Collision_Distance(cell_particles->quantized_collision_distance(k));
                Adjust_Particle_For_Objects(cell_particles->X(k),velocity,cell_particles->radius(k),collision_distance,particle_type,dt,time);
                cell_particles->X(k)+=dt*velocity;}
                cell_particles=cell_particles->next;}}}}
    else{
        T_LINEAR_INTERPOLATION_SCALAR linear_interpolation;T_FACE_LOOKUP V_lookup(V);
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(particles(block)){
            typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block);
            while(cell_particles){for(int k=1;k<=cell_particles->array_collection->Size();k++){
                TV velocity=linear_interpolation.Clamped_To_Array_Face(levelset.grid,V_lookup,cell_particles->X(k));
                if(enforce_domain_boundaries) levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(*cell_particles,k,velocity,particle_type,dt,time);
                T collision_distance=Particle_Collision_Distance(cell_particles->quantized_collision_distance(k));
                Adjust_Particle_For_Objects(cell_particles->X(k),velocity,cell_particles->radius(k),collision_distance,particle_type,dt,time);
                cell_particles->X(k)+=dt*velocity;}
                cell_particles=cell_particles->next;}}}}
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Second_Order_Runge_Kutta_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
    const bool update_particle_cells_after_euler_step,const bool verbose,const bool default_to_forward_euler_outside)
{
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(domain,particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(domain,thread_queue).template Run<const T_FACE_ARRAYS_SCALAR&,T_ARRAYS_PARTICLES&,const PARTICLE_LEVELSET_PARTICLE_TYPE,T,T,bool,bool>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Second_Order_Runge_Kutta_Step_Particles_Threaded,V,particles,particle_type,dt,time,verbose,default_to_forward_euler_outside);
    Consistency_Check(domain,particles);
    if(update_particle_cells_after_euler_step) Update_Particle_Cells(particles);
    Consistency_Check(domain,particles);
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Second_Order_Runge_Kutta_Step_Particles_Threaded(RANGE<TV_INT>& domain,const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool verbose,const bool default_to_forward_euler_outside)
{
    T_FACE_LOOKUP V_lookup(V);
    T_FACE_LOOKUP_COLLIDABLE V_lookup_collidable(V_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP V_lookup_collidable_lookup(V_lookup_collidable,V_lookup);
    T_LINEAR_INTERPOLATION_SCALAR linear_interpolation; // use for second step since particle may not be in initial block
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    COLLISION_GEOMETRY_ID body_id;int aggregate_id;T start_phi,end_phi;TV body_normal,body_velocity;int number_of_deleted_particles=0,number_of_non_occupied_cells=0,number_of_occupied_cells=0;
    T_ARRAYS_SCALAR phi_ghost;
    if(default_to_forward_euler_outside){
        phi_ghost.Resize(levelset.grid.Domain_Indices(number_of_ghost_cells));
        levelset.boundary->Fill_Ghost_Cells(levelset.grid,levelset.phi,phi_ghost,dt,time,number_of_ghost_cells);}

    const T one_over_dt=1/dt;

    int max_particles=0;
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block_index);
        max_particles=max(max_particles,cell_particles->array_collection->Size());
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        while(cell_particles){
            if(!levelset.collision_body_list || !levelset.collision_body_list->Swept_Occupied_Block(block)){ // not occupied, so advect ignoring objects
                number_of_non_occupied_cells++;
                T_LINEAR_INTERPOLATION_MAC_HELPER linear_interpolation_mac_helper(block,V);
                for(int k=1;k<=cell_particles->array_collection->Size();k++){
                    TV velocity=linear_interpolation_mac_helper.Interpolate_Face(cell_particles->X(k));
                    TV X_new=cell_particles->X(k)+dt*velocity;
                    if(default_to_forward_euler_outside){
                        T phi_X_new=linear_interpolation.Clamped_To_Array(levelset.grid,phi_ghost,X_new);
                        if(phi_X_new<3.0*levelset.grid.min_dX){
                            velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array_Face(levelset.grid,V_lookup,X_new));}}
                    else{
                        velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array_Face(levelset.grid,V_lookup,X_new));}
                    levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(*cell_particles,k,velocity,particle_type,dt,time);
                    cell_particles->X(k)+=dt*velocity;}}
            else{ // collision aware advection
                number_of_occupied_cells++;
                for(int k=cell_particles->array_collection->Size();k>=1;k--){
                    T collision_distance=Particle_Collision_Distance(cell_particles->quantized_collision_distance(k));
                    // first push particle out
                    if(particle_type==PARTICLE_LEVELSET_NEGATIVE){bool particle_crossover;
                        if(levelset.collision_body_list->Push_Out_Point(cell_particles->X(k),collision_distance,true,particle_crossover)&&particle_crossover){
                            Delete_Particle_And_Clean_Memory(particles(block_index),*cell_particles,k);number_of_deleted_particles++;continue;}} // delete due to crossover in push out
                    // TODO: optimize by calling From_Base_Node (or something like that) since we know what cell we're in
                    TV velocity=linear_interpolation_collidable.Clamped_To_Array_Face(levelset.grid,V_lookup_collidable_lookup,cell_particles->X(k));
                    // adjust normal component of velocity to not move into nearest body
                    bool got_interaction_with_body=false;
                    if(particle_type==PARTICLE_LEVELSET_NEGATIVE)
                        got_interaction_with_body=levelset.collision_body_list->Get_Body_Penetration(cell_particles->X(k),cell_particles->X(k),collision_distance,dt,body_id,aggregate_id,
                            start_phi,end_phi,body_normal,body_velocity);
                    if(got_interaction_with_body){
                        T relative_normal_velocity=TV::Dot_Product(velocity-body_velocity,body_normal);
                        if(relative_normal_velocity<0) velocity-=relative_normal_velocity*body_normal;}
                    // compute position using velocity but clamp if intersects any body
                    TV X_new=cell_particles->X(k)+dt*velocity;
                    RAY<TV> ray;
                    if(RAY<TV>::Create_Non_Degenerate_Ray(cell_particles->X(k),X_new-cell_particles->X(k),ray) && levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id))
                        X_new=ray.Point(ray.t_max);
                    // get average velocity
                    if(default_to_forward_euler_outside){
                        T phi_X_new=linear_interpolation.Clamped_To_Array(levelset.grid,phi_ghost,X_new);
                        if(phi_X_new<3.0*levelset.grid.min_dX){
                            velocity=(T).5*(velocity+linear_interpolation_collidable.Clamped_To_Array_Face(levelset.grid,V_lookup_collidable_lookup,X_new));}}
                    else{
                        velocity=(T).5*(velocity+linear_interpolation_collidable.Clamped_To_Array_Face(levelset.grid,V_lookup_collidable_lookup,X_new));}
                    // adjust normal component of velocity again
                    if(got_interaction_with_body){
                        T relative_normal_velocity=TV::Dot_Product(velocity-body_velocity,body_normal);
                        if(relative_normal_velocity<0) velocity-=relative_normal_velocity*body_normal;}
                    levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(*cell_particles,k,velocity,particle_type,dt,time);
                    velocity=(X_new-cell_particles->X(k))*one_over_dt;
                    if(!Adjust_Particle_For_Objects(cell_particles->X(k),velocity,cell_particles->radius(k),collision_distance,particle_type,dt,time)){
                        Delete_Particle_And_Clean_Memory(particles(block_index),*cell_particles,k);number_of_deleted_particles++;}
                    else cell_particles->X(k)+=dt*velocity;}}
            cell_particles=cell_particles->next;}
        if(particles(block_index)->array_collection->Size()==0) Free_Particle_And_Clear_Pointer(particles(block_index));
        }}
    if(verbose){
        std::stringstream ss;
        ss<<"max_particles_per_cell "<<max_particles<<std::endl;
        if(number_of_deleted_particles) ss<<"Deleted "<<number_of_deleted_particles<<" "<<PARTICLE_LEVELSET<T_GRID>::Particle_Type_Name(particle_type)<<" due to crossover"<<std::endl;
        if(number_of_non_occupied_cells) ss<<number_of_occupied_cells<<" occupied_cells "<<", "<<number_of_non_occupied_cells<<" non occupied cells"<<std::endl;
        LOG::filecout(ss.str());}
}
//#####################################################################
// Function Update_Particle_Cells
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Update_Particle_Cells(T_ARRAYS_PARTICLES& particles)
{
    typedef typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE T_PARTICLES;
    const T_PARTICLES& template_particles=choice<(2-IS_SAME<T_PARTICLES,PARTICLE_LEVELSET_PARTICLES<TV> >::value)>(this->template_particles,this->template_removed_particles);
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(domain,particles);
    if(mpi_grid){
        Exchange_Boundary_Particles(*mpi_grid,template_particles,particles,(int)(cfl_number+(T)1.5),*this);}
    Consistency_Check(domain,particles);
    if(thread_queue){
        DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV> threaded_iterator(domain,thread_queue);

        ARRAY<ARRAY<int>,TV_INT> number_of_particles_per_block(levelset.grid.Block_Indices());
        threaded_iterator.template Run<T_ARRAYS_PARTICLES&,ARRAY<ARRAY<int>,TV_INT>&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Update_Particle_Cells_Part_One_Threaded,particles,number_of_particles_per_block);
        
        ARRAY<int,TV_INT> domain_index(domain.Thickened(1));
        for(int i=1;i<=threaded_iterator.domains.m;i++)
            for(NODE_ITERATOR iterator(levelset.grid,threaded_iterator.domains(i));iterator.Valid();iterator.Next())
                domain_index(iterator.Node_Index())=i; // TODO: thread
        
        ARRAY<ARRAY<TRIPLE<TV_INT,T_PARTICLES*,int> >,TV_INT> list_to_process(domain.Thickened(1));
        Update_Particle_Cells_Find_List(domain,particles,number_of_particles_per_block,list_to_process); //This can not be threaded and be deterministic
        threaded_iterator.template Run<T_ARRAYS_PARTICLES&,const ARRAY<ARRAY<int>,TV_INT>&,ARRAY<ARRAY<TRIPLE<TV_INT,T_PARTICLES*,int> >,TV_INT>&,const ARRAY<int,TV_INT>*>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Update_Particle_Cells_Part_Two_Threaded,particles,number_of_particles_per_block,list_to_process,&domain_index);
        threaded_iterator.template Run<T_ARRAYS_PARTICLES&,const ARRAY<ARRAY<int>,TV_INT>&,ARRAY<ARRAY<TRIPLE<TV_INT,T_PARTICLES*,int> >,TV_INT>&,const ARRAY<int,TV_INT>*>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Update_Particle_Cells_Part_Three_Threaded,particles,number_of_particles_per_block,list_to_process,&domain_index);
        threaded_iterator.template Run<T_ARRAYS_PARTICLES&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Marked_Particles,particles);}
    else{
        ARRAY<ARRAY<int>,TV_INT> number_of_particles_per_block(levelset.grid.Block_Indices());
        ARRAY<ARRAY<TRIPLE<TV_INT,T_PARTICLES*,int> >,TV_INT> list_to_process(domain.Thickened(1));
        Update_Particle_Cells_Part_One_Threaded(domain,particles,number_of_particles_per_block);
        Update_Particle_Cells_Find_List(domain,particles,number_of_particles_per_block,list_to_process);
        Update_Particle_Cells_Part_Two_Threaded(domain,particles,number_of_particles_per_block,list_to_process,0);
        Update_Particle_Cells_Part_Three_Threaded(domain,particles,number_of_particles_per_block,list_to_process,0);
        Delete_Marked_Particles(domain,particles);}
    Consistency_Check(domain,particles);
}
//#####################################################################
// Function Update_Particle_Cells_Threaded
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Update_Particle_Cells_Part_One_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block)
{
    typedef typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE T_PARTICLES;
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block);
        while(cell_particles){number_of_particles_per_block(block).Append(cell_particles->array_collection->Size());cell_particles=cell_particles->next;}}
}
//#####################################################################
// Function Update_Particle_Cells_Threaded
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Update_Particle_Cells_Find_List(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process)
{
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
    typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block);
    for(int i=1;i<=number_of_particles_per_block(block).m;i++){
        for(int k=1;k<=number_of_particles_per_block(block)(i);k++){
            TV_INT final_block=levelset.grid.Block_Index(cell_particles->X(k),1);
            if(final_block!=block){
                list_to_process(final_block).Append(TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int>(block,cell_particles,k));}}
        cell_particles=cell_particles->next;}}
}
//#####################################################################
// Function Update_Particle_Cells_Threaded
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Update_Particle_Cells_Part_Two_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index)
{
    typedef typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE T_PARTICLES;
    const T_PARTICLES& template_particles=choice<(2-IS_SAME<T_PARTICLES,PARTICLE_LEVELSET_PARTICLES<TV> >::value)>(this->template_particles,this->template_removed_particles);
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT final_block=iterator.Node_Index();
        for(int i=1;i<=list_to_process(final_block).m;i++){
            T_PARTICLES* cell_particles=list_to_process(final_block)(i).y;int k=list_to_process(final_block)(i).z;
            if(!particles(final_block))particles(final_block)=Allocate_Particles(template_particles);
            Copy_Particle(*cell_particles,*particles(final_block),k);}}
}
//#####################################################################
// Function Update_Particle_Cells_Threaded
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Update_Particle_Cells_Part_Three_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index)
{
    typedef typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE T_PARTICLES;
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT final_block=iterator.Node_Index();
        for(int i=1;i<=list_to_process(final_block).m;i++){
            T_PARTICLES* cell_particles=list_to_process(final_block)(i).y;int k=list_to_process(final_block)(i).z;
            assert(cell_particles->radius(k)>0);
            cell_particles->radius(k)=-(T)1;}}
}
//#####################################################################
// Function Delete_Marked_Particles
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Marked_Particles(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles)
{
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(particles(block)){
        typename T_ARRAYS_PARTICLES::ELEMENT cell_particles=particles(block);
        while(1){
            for(int k=1;k<=cell_particles->array_collection->Size();k++) if(cell_particles->radius(k)<0){
                Delete_Particle_And_Clean_Memory(particles(block),*cell_particles,k);k--;}
            if(!cell_particles->next)break;
            cell_particles=cell_particles->next;}
        if(particles(block)->array_collection->Size()==0) Free_Particle_And_Clear_Pointer(particles(block));}}
}
//#####################################################################
// Function Advect_Particles_In_Ghost_Region
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Advect_Particles_In_Ghost_Region(T_GRID& grid,T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T dt,int number_of_ghost_cells)
{
    // use second order Runge Kutta
    T_FACE_LOOKUP V_lookup(V);
    T_LINEAR_INTERPOLATION_SCALAR linear_interpolation;

    for(NODE_ITERATOR iterator(grid,number_of_ghost_cells-1,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        TV_INT block_index=iterator.Node_Index();
        if(particles(block_index)){
            PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=particles(block_index);
            BLOCK_UNIFORM<T_GRID> block(grid,block_index);
            while(cell_particles){
                T_LINEAR_INTERPOLATION_MAC_HELPER linear_interpolation_mac_helper(block,V);
                for(int k=1;k<=cell_particles->array_collection->Size();k++){
                    TV velocity=linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,cell_particles->X(k));
                    TV X_new=cell_particles->X(k)+dt*velocity;
                    velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,X_new));
                    cell_particles->X(k)+=dt*velocity;
                    if(grid.Domain().Lazy_Inside(cell_particles->X(k))){
                        TV_INT to_block=grid.Block_Index(cell_particles->X(k),1);
                        if(!particles(to_block))particles(to_block)=Allocate_Particles(template_particles);
                        Copy_Particle(*cell_particles,*particles(to_block),k);}}
                cell_particles=cell_particles->next;}
            if(particles(block_index)->array_collection->Size()==0) Free_Particle_And_Clear_Pointer(particles(block_index));}}
}
//#####################################################################
// Function Reseed_Particles_In_Ghost_Region
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Particles_In_Ghost_Region(const T time,int number_of_ghost_cells)
{
    int new_particles=0;
    bool normals_defined=(levelset.normals!=0);levelset.Compute_Normals(); // make sure normals are accurate
    
    T_ARRAYS_BOOL cell_centered_mask(levelset.grid.Domain_Indices(number_of_ghost_cells));
    cell_centered_mask.Fill(false);
    for(CELL_ITERATOR iterator(levelset.grid,number_of_ghost_cells,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        assert(!levelset.grid.domain.Lazy_Inside(iterator.Location()));
        cell_centered_mask(iterator.Cell_Index())=true;}

    new_particles-=Reseed_Delete_Particles_In_Whole_Region(number_of_ghost_cells-1,negative_particles,-1,&cell_centered_mask);
    new_particles+=Reseed_Add_Particles_In_Whole_Region(number_of_ghost_cells-1,negative_particles,positive_particles,-1,PARTICLE_LEVELSET_NEGATIVE,time,&cell_centered_mask);

    if(!only_use_negative_particles){
        new_particles-=Reseed_Delete_Particles_In_Whole_Region(number_of_ghost_cells-1,positive_particles,1,&cell_centered_mask);
        new_particles+=Reseed_Add_Particles_In_Whole_Region(number_of_ghost_cells-1,positive_particles,negative_particles,1,PARTICLE_LEVELSET_POSITIVE,time,&cell_centered_mask);
    }

    if(normals_defined){delete levelset.normals;levelset.normals=0;}
    return new_particles;
}
//#####################################################################
// Function Reseed_Delete_Particles_In_Whole_Region
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Delete_Particles_In_Whole_Region(int number_of_ghost_cells,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign,T_ARRAYS_BOOL* cell_centered_mask)
{
    int number_deleted=0;
    for(NODE_ITERATOR iterator(levelset.grid,number_of_ghost_cells,T_GRID::WHOLE_REGION);iterator.Valid();iterator.Next()){
        TV_INT block_index=iterator.Node_Index();
        if(particles(block_index)){
            ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> > deletion_list;
            PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(block_index);
            BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
            for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++){
                TV_INT cell=block.Cell(cell_index);
                if(!levelset.phi.domain.Lazy_Inside(cell)) continue;
                T unsigned_phi=sign*levelset.phi(cell);
                if(0<=unsigned_phi && unsigned_phi<=half_band_width && (!cell_centered_mask||(*cell_centered_mask)(cell))){
                    if(cell_particles.array_collection->Size()==0){Free_Particle_And_Clear_Pointer(particles(block_index));}
                    else{
                        PARTICLE_LEVELSET_PARTICLES<TV>* local_cell_particles=&cell_particles;
                        while(local_cell_particles){
                            for(int index=1;index<=local_cell_particles->array_collection->Size();index++) 
                                Add_Particle_To_Deletion_List(deletion_list,*local_cell_particles,index);
                            local_cell_particles=local_cell_particles->next;}
                        Delete_Particles_From_Deletion_List(deletion_list,cell_particles);
                        if(cell_particles.array_collection->Size()==0) 
                            Free_Particle_And_Clear_Pointer(particles(block_index));}}}}}
    return number_deleted;
}
//#####################################################################
// Function Reseed_Add_Particles_In_Whole_Region
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Add_Particles_In_Whole_Region(int number_of_ghost_cells,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
    const T time,T_ARRAYS_BOOL* cell_centered_mask)
{
    int number_added=0;
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices(number_of_ghost_cells));domain.max_corner+=TV_INT::All_Ones_Vector();
   
    Consistency_Check(domain,particles);
    ARRAY<int,TV_INT> number_of_particles_to_add(domain);
  
    ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT> list_to_process(domain.Thickened(1));

    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){
        TV_INT block_index=iterator.Node_Index();
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++){
            TV_INT cell=block.Cell(cell_index);
            if(!levelset.phi.domain.Lazy_Inside(cell)) continue;
            T unsigned_phi=sign*levelset.phi(cell);
            if(0<=unsigned_phi && unsigned_phi<=half_band_width && (!cell_centered_mask||(*cell_centered_mask)(cell))){
                if(!particles(block_index)) particles(block_index)=Allocate_Particles(template_particles);
                int total_particles=0,total_other_particles=0;
                PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=particles(block_index);
                while(cell_particles){total_particles+=cell_particles->array_collection->Size();cell_particles=cell_particles->next;}
                cell_particles=other_particles(block_index);
                while(cell_particles){total_other_particles+=cell_particles->array_collection->Size();cell_particles=cell_particles->next;}
                if(total_particles+total_other_particles>=number_particles_per_cell){ 
                    if(particles(block_index)->array_collection->Size()==0) 
                        Free_Particle_And_Clear_Pointer(particles(block_index));
                    continue;}
                cell_particles=particles(block_index);
                number_of_particles_to_add(block_index)=number_particles_per_cell-total_particles-total_other_particles;
                if(sign==-1){ // we add the negative particles first, and don't want to add too many of them...
                    if(total_other_particles) number_of_particles_to_add(block_index)=(int)((T)number_of_particles_to_add(block_index)*(T)total_particles/(T)(total_particles+total_other_particles)+1);
                    else if(!total_particles) number_of_particles_to_add(block_index)=(number_of_particles_to_add(block_index)+1)/2;}}}}
    
    RANDOM_NUMBERS<T> local_random;
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){
        TV_INT block_index=iterator.Node_Index();
        if(!number_of_particles_to_add(block_index)) continue;
        VECTOR<T,TV::dimension+1> h;for(int axis=1;axis<=TV::dimension;++axis) h(axis)=(T)block_index(axis);h(TV::dimension+1) = time;
        local_random.Set_Seed(Hash(h));
        ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> > deletion_list;
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=particles(block_index);
        T phi_min=sign*minimum_particle_radius,phi_max=sign*half_band_width;if(phi_min>phi_max) exchange(phi_min,phi_max);
        RANGE<TV> block_bounding_box=block.Bounding_Box();
        //int attempts=0;
        ARRAY_VIEW<int>* id=store_unique_particle_id?cell_particles->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID):0;
        for(int k=1;k<=number_of_particles_to_add(block_index);k++){
            PARTICLE_LEVELSET_PARTICLES<TV>* local_cell_particles=cell_particles;
            int index=Add_Particle(cell_particles);
            if(local_cell_particles!=cell_particles){
                if(store_unique_particle_id) id=cell_particles->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);}
#ifndef WIN32
            if(id) (*id)(index)= __exchange_and_add(&last_unique_particle_id,1);
#else
            PHYSBAM_ASSERT(!thread_queue);
            if(id) (*id)(index)=++last_unique_particle_id;
#endif
            cell_particles->quantized_collision_distance(index)=(unsigned short)(local_random.Get_Number()*USHRT_MAX);
            cell_particles->X(index)=local_random.Get_Uniform_Vector(block_bounding_box);
            //attract particles to the interface
            T phi=levelset.Phi(cell_particles->X(index));bool inside=(phi>=phi_min && phi<=phi_max);
            if(!inside){
                T phi_goal=local_random.Get_Uniform_Number(phi_min,phi_max);TV X,N;int iteration=0;
                while(!inside && iteration<maximum_iterations_for_attraction){
                    N=levelset.Normal(cell_particles->X(index));T distance=phi_goal-phi,dt=1;bool inside_domain=false;
                    while(!inside_domain && iteration<maximum_iterations_for_attraction){
                        X=cell_particles->X(index)+dt*distance*N;
                        if(levelset.grid.Ghost_Domain(number_of_ghost_cells).Lazy_Outside(X)){dt*=.5;iteration++;}else inside_domain=true;}
                    if(!inside_domain) break; // ran out of iterations
                    phi=levelset.Phi(X);inside=(phi>=phi_min && phi<=phi_max);
                    if(!inside){
                        dt*=.5;iteration++;X=cell_particles->X(index)+dt*distance*N;
                        phi=levelset.Phi(X);inside=(phi>=phi_min && phi<=phi_max);}
                    cell_particles->X(index)=X;}}
            if(!inside) Add_Particle_To_Deletion_List(deletion_list,*cell_particles,index);
            else cell_particles->radius(index)=clamp(abs(phi),minimum_particle_radius,maximum_particle_radius);}
        Delete_Particles_From_Deletion_List(deletion_list,*particles(block_index));
        if(particles(block_index)->array_collection->Size()==0) Free_Particle_And_Clear_Pointer(particles(block_index));}
    
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()) number_added+=number_of_particles_to_add(iterator.Node_Index());
    return number_added;
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Particles(const T time,T_ARRAYS_BOOL* cell_centered_mask)
{
    int new_particles=0;
    bool normals_defined=(levelset.normals!=0);levelset.Compute_Normals(); // make sure normals are accurate

    if(!cell_centered_mask) new_particles-=Reseed_Delete_Particles(negative_particles,-1);
    new_particles+=Reseed_Add_Particles(negative_particles,positive_particles,-1,PARTICLE_LEVELSET_NEGATIVE,time,cell_centered_mask);

    if(!only_use_negative_particles){
        if(!cell_centered_mask) new_particles-=Reseed_Delete_Particles(positive_particles,1);
        new_particles+=Reseed_Add_Particles(positive_particles,negative_particles,1,PARTICLE_LEVELSET_POSITIVE,time,cell_centered_mask);}

    if(!normals_defined){delete levelset.normals;levelset.normals=0;}
    return new_particles;
}
//#####################################################################
// Function Reseed_Delete_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Delete_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign)
{
    int number_deleted=0;ARRAY<int> heap_particle_indices(number_particles_per_cell);ARRAY<T> heap_phi_minus_radius(number_particles_per_cell);
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> > deletion_list;
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(block_index);
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++){TV_INT cell=block.Cell(cell_index);
            T unsigned_phi=sign*levelset.phi(cell);if(0<=unsigned_phi && unsigned_phi<=half_band_width) goto NEAR_THE_INTERFACE;}
        
        //TODO (mlentine): Is this right?
        if(cell_particles.array_collection->Size()==0){Free_Particle_And_Clear_Pointer(particles(block_index));}
        else{ // delete all non-escaped particles
            PARTICLE_LEVELSET_PARTICLES<TV>* local_cell_particles=&cell_particles;
            while(local_cell_particles){
                ARRAY<bool> escaped;Identify_Escaped_Particles(block,*local_cell_particles,escaped,sign);
                for(int index=1;index<=local_cell_particles->array_collection->Size();index++) if(!escaped(index)) Add_Particle_To_Deletion_List(deletion_list,*local_cell_particles,index);
                local_cell_particles=local_cell_particles->next;}
            Delete_Particles_From_Deletion_List(deletion_list,cell_particles);
            if(cell_particles.array_collection->Size()==0){Free_Particle_And_Clear_Pointer(particles(block_index));}}
        continue;
        
        NEAR_THE_INTERFACE:; // can only get here via the goto
        if(cell_particles.array_collection->Size()<=number_particles_per_cell) continue; //This assumes particle_pool.number_particles_per_cell>number_particles_per_cell
        ARRAY<ARRAY<bool> > all_escaped;
        int number_of_escaped_particles=0,total_particles=0;
        PARTICLE_LEVELSET_PARTICLES<TV>* local_cell_particles=&cell_particles;
        while(local_cell_particles){
            ARRAY<bool> escaped;Identify_Escaped_Particles(block,cell_particles,escaped,sign);number_of_escaped_particles+=escaped.Number_True();all_escaped.Append(escaped);
            total_particles+=local_cell_particles->array_collection->Size();local_cell_particles=local_cell_particles->next;}
        total_particles-=number_of_escaped_particles;
        if(total_particles>number_particles_per_cell){ // too many particles - delete particles with a heap sort
            number_deleted+=total_particles-number_particles_per_cell;int heap_size=0;
            local_cell_particles=&cell_particles;
            for(int i=1;local_cell_particles;i++){for(int index=1;index<=local_cell_particles->array_collection->Size();index++) if(!all_escaped(i)(index)){
                T phi_minus_radius=sign*levelset.Phi(local_cell_particles->X(index))-local_cell_particles->radius(index);
                if(heap_size<number_particles_per_cell){ // add particle to heap
                    heap_size++;heap_particle_indices(heap_size)=index+particle_pool.number_particles_per_cell*(i-1);heap_phi_minus_radius(heap_size)=phi_minus_radius;
                    if(heap_size==number_particles_per_cell) ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices);} // when heap is full, order values with largest on top
                else{ // excess particles don't fit in the heap
                    if(phi_minus_radius<heap_phi_minus_radius(1)){ // delete particle on top of heap & add new particle
                        int deletion_index;
                        PARTICLE_LEVELSET_PARTICLES<TV>& deletion_particle_list=Get_Particle_Link(cell_particles,heap_particle_indices(1),deletion_index);
                        Add_Particle_To_Deletion_List(deletion_list,deletion_particle_list,deletion_index); 
                        heap_phi_minus_radius(1)=phi_minus_radius;heap_particle_indices(1)=index+particle_pool.number_particles_per_cell*(i-1);
                        ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices,1,heap_phi_minus_radius.m);}
                    else Add_Particle_To_Deletion_List(deletion_list,*local_cell_particles,index);}} // delete new particle, larger than top of heap
            local_cell_particles=local_cell_particles->next;}
            Delete_Particles_From_Deletion_List(deletion_list,cell_particles);
            if(cell_particles.array_collection->Size()==0){Free_Particle_And_Clear_Pointer(particles(block_index));}}}}
    return number_deleted;
}
//#####################################################################
// Function Reseed_Add_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Add_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
    const T time,T_ARRAYS_BOOL* cell_centered_mask)
{
    int number_added=0;
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    //ARRAY<ARRAY<TRIPLE<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT,int> >,TV_INT> move_particles(domain);
    Consistency_Check(domain,particles);
    ARRAY<int,TV_INT> number_of_particles_to_add(domain);
    if(thread_queue){
        DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV> threaded_iterator(domain,thread_queue);
        threaded_iterator.template Run<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&,int,T_ARRAYS_BOOL*,ARRAY<int,TV_INT>&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Reseed_Add_Particles_Threaded_Part_One,particles,other_particles,sign,cell_centered_mask,number_of_particles_to_add);

        ARRAY<int,TV_INT> domain_index(domain);
        for(int i=1;i<=threaded_iterator.domains.m;i++)
            for(NODE_ITERATOR iterator(levelset.grid,threaded_iterator.domains(i));iterator.Valid();iterator.Next())
                domain_index(iterator.Node_Index())=i; // TODO: thread

        ARRAY<ARRAY<int>,TV_INT> number_of_particles_per_block(levelset.grid.Block_Indices());
        ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT> list_to_process(domain.Thickened(1));
        threaded_iterator.template Run<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&,int,const PARTICLE_LEVELSET_PARTICLE_TYPE,T,const ARRAY<int,TV_INT>&,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>&,const ARRAY<int,TV_INT>*>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Reseed_Add_Particles_Threaded_Part_Two,particles,sign,particle_type,time,number_of_particles_to_add,list_to_process,&domain_index);
        threaded_iterator.template Run<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&,const ARRAY<ARRAY<int>,TV_INT>&,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>&,const ARRAY<int,TV_INT>*>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Update_Particle_Cells_Part_Two_Threaded,particles,number_of_particles_per_block,list_to_process,&domain_index);
        threaded_iterator.template Run<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&,const ARRAY<ARRAY<int>,TV_INT>&,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>&,const ARRAY<int,TV_INT>*>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Update_Particle_Cells_Part_Three_Threaded,particles,number_of_particles_per_block,list_to_process,&domain_index);
        threaded_iterator.template Run<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Marked_Particles,particles);}
    else{
        ARRAY<ARRAY<int>,TV_INT> number_of_particles_per_block(levelset.grid.Block_Indices());
        ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT> list_to_process(domain.Thickened(1));
        Reseed_Add_Particles_Threaded_Part_One(domain,particles,other_particles,sign,cell_centered_mask,number_of_particles_to_add);
        Reseed_Add_Particles_Threaded_Part_Two(domain,particles,sign,particle_type,time,number_of_particles_to_add,list_to_process,0);
        Update_Particle_Cells_Part_Two_Threaded(domain,particles,number_of_particles_per_block,list_to_process,0);
        Update_Particle_Cells_Part_Three_Threaded(domain,particles,number_of_particles_per_block,list_to_process,0);
        Delete_Marked_Particles(domain,particles);}
    Consistency_Check(domain,particles);
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()) number_added+=number_of_particles_to_add(iterator.Node_Index());
    return number_added;
}
//#####################################################################
// Function Reseed_Add_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Add_Particles_Threaded_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,T_ARRAYS_BOOL* cell_centered_mask,ARRAY<int,TV_INT>& number_of_particles_to_add)
{
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++){TV_INT cell=block.Cell(cell_index);
            T unsigned_phi=sign*levelset.phi(cell);
            if(0<=unsigned_phi && unsigned_phi<=half_band_width && (!cell_centered_mask||(*cell_centered_mask)(cell))) goto NEAR_THE_INTERFACE;}
        continue;

        NEAR_THE_INTERFACE:; // can only get here via the goto 
        if(!particles(block_index)) particles(block_index)=Allocate_Particles(template_particles);
        int total_particles=0,total_other_particles=0;
        PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=particles(block_index);
        while(cell_particles){total_particles+=cell_particles->array_collection->Size();cell_particles=cell_particles->next;}
        cell_particles=other_particles(block_index);
        while(cell_particles){total_other_particles+=cell_particles->array_collection->Size();cell_particles=cell_particles->next;}
        if(total_particles+total_other_particles>=number_particles_per_cell) continue;
        cell_particles=particles(block_index);
        number_of_particles_to_add(block_index)=number_particles_per_cell-total_particles-total_other_particles;
        if(sign==-1){ // we add the negative particles first, and don't want to add too many of them...
            if(total_other_particles) number_of_particles_to_add(block_index)=(int)((T)number_of_particles_to_add(block_index)*(T)total_particles/(T)(total_particles+total_other_particles)+1);
            else if(!total_particles) number_of_particles_to_add(block_index)=(number_of_particles_to_add(block_index)+1)/2;}}
}
//#####################################################################
// Function Reseed_Add_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reseed_Add_Particles_Threaded_Part_Two(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const ARRAY<int,TV_INT>& number_of_particles_to_add,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index)
{
    RANDOM_NUMBERS<T> local_random;
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();
        if(!number_of_particles_to_add(block_index)) continue;
        VECTOR<T,TV::dimension+1> h;for(int axis=1;axis<=TV::dimension;++axis) h(axis)=(T)block_index(axis);h(TV::dimension+1) = time;
        local_random.Set_Seed(Hash(h));
        BLOCK_UNIFORM<T_GRID> block(levelset.grid,block_index);
        PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=particles(block_index);
        T phi_min=sign*minimum_particle_radius,phi_max=sign*half_band_width;if(phi_min>phi_max) exchange(phi_min,phi_max);
        RANGE<TV> block_bounding_box=block.Bounding_Box();
        int attempts=0;
        ARRAY_VIEW<int>* id=store_unique_particle_id?cell_particles->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID):0;
        ARRAY_VIEW<T>* material_volume=vof_advection?cell_particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME):0;
        for(int k=1;k<=number_of_particles_to_add(block_index);k++){
            PARTICLE_LEVELSET_PARTICLES<TV>* local_cell_particles=cell_particles;
            int index=Add_Particle(cell_particles);
            if(local_cell_particles!=cell_particles){
                if(store_unique_particle_id) id=cell_particles->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);
                if(vof_advection) material_volume=cell_particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME);}
#ifndef WIN32
            if(id) (*id)(index)= __exchange_and_add(&last_unique_particle_id,1);
#else
            PHYSBAM_ASSERT(!thread_queue);
            if(id) (*id)(index)=++last_unique_particle_id;
#endif
            if(material_volume) (*material_volume)(index)=0;
            cell_particles->quantized_collision_distance(index)=(unsigned short)(local_random.Get_Number()*USHRT_MAX);
            cell_particles->X(index)=local_random.Get_Uniform_Vector(block_bounding_box);
            cell_particles->radius(index)=1;//Default value
            if(!Attract_Individual_Particle_To_Interface_And_Adjust_Radius(particles,*cell_particles,phi_min,phi_max,block,index,particle_type,time,true,local_random,list_to_process,domain_index) && ++attempts<=5) k--;}}
}
//#####################################################################
// Function Identify_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Identify_Escaped_Particles(RANGE<TV_INT>& domain,const int sign)
{
    if(!sign || sign==1){
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
            if(positive_particles(block)) Identify_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*positive_particles(block),escaped_positive_particles(block),1);
            else escaped_positive_particles(block).Resize(0);}}
    if(!sign || sign==-1){
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
            if(negative_particles(block)) Identify_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*negative_particles(block),escaped_negative_particles(block),-1);
            else escaped_negative_particles(block).Resize(0);}}
}
//#####################################################################
// Function Identify_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Identify_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign)
{
    bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block):false;
    bool enable_collision_aware_already = (levelset.collision_unaware_interpolation)==0?false:true;
    if(near_objects && !enable_collision_aware_already) {
        levelset.Enable_Collision_Aware_Interpolation(sign);}
    escaped.Resize(particles.array_collection->Size());ARRAYS_COMPUTATIONS::Fill(escaped,false);T one_over_radius_multiplier=-sign/outside_particle_distance_multiplier;
    for(int k=1;k<=particles.array_collection->Size();k++) if(one_over_radius_multiplier*levelset.Phi(particles.X(k))>particles.radius(k)) escaped(k)=true;
    if(near_objects && !enable_collision_aware_already) {
        levelset.Disable_Collision_Aware_Interpolation();}
}
//#####################################################################
// Function Delete_All_Particles_In_Cell
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_All_Particles_In_Cell(const TV_INT& block)
{
    if(positive_particles.Valid_Index(block)) Free_Particle_And_Clear_Pointer(positive_particles(block));
    if(negative_particles.Valid_Index(block)) Free_Particle_And_Clear_Pointer(negative_particles(block));
    if(use_removed_positive_particles && removed_positive_particles.Valid_Index(block)){delete removed_positive_particles(block);removed_positive_particles(block)=0;}
    if(use_removed_negative_particles && removed_negative_particles.Valid_Index(block)){delete removed_negative_particles(block);removed_negative_particles(block)=0;}
}
//#####################################################################
// Function Delete_Positive_Particles_In_Cell
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Positive_Particles_In_Cell(const TV_INT& block)
{
    if(positive_particles.Valid_Index(block)) Free_Particle_And_Clear_Pointer(positive_particles(block));
    if(use_removed_positive_particles && removed_positive_particles.Valid_Index(block)){delete removed_positive_particles(block);removed_positive_particles(block)=0;}
}
//#####################################################################
// Function Delete_Deep_Escaped_Particles
//#####################################################################
// delete those farther than radius_fraction too deep
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Deep_Escaped_Particles(const T radius_fraction,const bool need_to_identify_escaped_particles,const bool verbose)
{
    std::stringstream ss;
    if(verbose) ss << "positive particles - ";int number_of_positive_particles_deleted=0;
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(positive_particles(block)){
        number_of_positive_particles_deleted+=Delete_Deep_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*positive_particles(block),escaped_positive_particles(block),1,radius_fraction,
            need_to_identify_escaped_particles);
        if(positive_particles(block)->array_collection->Size()==0) Free_Particle_And_Clear_Pointer(positive_particles(block));}}
    if(verbose) ss << "deleted " << number_of_positive_particles_deleted << " particles" << std::endl;
    if(verbose) ss << "negative particles - ";int number_of_negative_particles_deleted=0;
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(negative_particles(block)){
        number_of_negative_particles_deleted+=Delete_Deep_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*negative_particles(block),escaped_negative_particles(block),-1,radius_fraction,
            need_to_identify_escaped_particles);
        if(negative_particles(block)->array_collection->Size()==0) Free_Particle_And_Clear_Pointer(negative_particles(block));}}
    if(verbose) ss << "deleted " << number_of_negative_particles_deleted << " particles" << std::endl;
    LOG::filecout(ss.str());
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();assert(deletion_list(block).m==0);}
}
//#####################################################################
// Function Delete_Deep_Escaped_Particles
//#####################################################################
// delete those farther than radius_fraction too deep
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Deep_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign,const T radius_fraction,
    const bool need_to_identify_escaped_particles)
{
    PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=&particles;
    int deleted=0;int minus_sign=-sign;
    while(cell_particles){
        if(need_to_identify_escaped_particles) Identify_Escaped_Particles(block,*cell_particles,escaped,sign);
        for(int k=cell_particles->array_collection->Size();k>=1;k--) if(escaped(k) && minus_sign*levelset.Phi(particles.X(k))>radius_fraction*particles.radius(k)){Add_Particle_To_Deletion_List(this->deletion_list(block.block_index),*cell_particles,k);deleted++;}}
    Delete_Particles_From_Deletion_List(this->deletion_list(block.block_index),particles); //TODO: Make these deletion list's per domain
    return deleted;
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Outside_Grid()
{
    RANGE<TV_INT> local_domain(levelset.grid.Domain_Indices());local_domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(local_domain,positive_particles);
    Consistency_Check(local_domain,negative_particles);

    RANGE<TV_INT> real_domain(levelset.grid.Domain_Indices());real_domain.max_corner+=TV_INT::All_Ones_Vector();
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices(3));domain.max_corner+=TV_INT::All_Ones_Vector();
    for(int axis=1;axis<=TV::dimension;axis++) for(int side=1;side<=2;side++){
        RANGE<TV_INT> ghost_domain(domain);
        if(side==1) ghost_domain.max_corner(axis)=real_domain.min_corner(axis)-1;
        else ghost_domain.min_corner(axis)=real_domain.max_corner(axis)+1;
        DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(ghost_domain,thread_queue).Run(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Particles_Far_Outside_Grid);}
    for(int axis=1;axis<=TV::dimension;axis++) for(int side=1;side<=2;side++){
        RANGE<TV_INT> boundary_domain(real_domain);
        if(side==1) boundary_domain.max_corner(axis)=real_domain.min_corner(axis);
        else boundary_domain.min_corner(axis)=real_domain.max_corner(axis);
        DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(boundary_domain,thread_queue).template Run<int,int>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Particles_Near_Outside_Grid,axis,side);}
    
    Consistency_Check(local_domain,positive_particles);
    Consistency_Check(local_domain,negative_particles);
}
//#####################################################################
// Function Delete_Particles_Far_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Far_Outside_Grid(RANGE<TV_INT>& domain)
{
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()) Delete_All_Particles_In_Cell(iterator.Node_Index());
}
//#####################################################################
// Function Delete_Particles_Near_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Near_Outside_Grid(RANGE<TV_INT>& domain,int axis,int side)
{
    // delete particles barely outside grid
    RANGE<TV> local_domain=levelset.grid.domain;TV domain_boundaries[2]={local_domain.Minimum_Corner(),local_domain.Maximum_Corner()};
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(negative_particles(block)){
            Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*negative_particles(block),2*side-3);
            if(!negative_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(negative_particles(block));}
        if(positive_particles(block)){
            Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*positive_particles(block),2*side-3);
            if(!positive_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(positive_particles(block));}
        if(use_removed_negative_particles)if(removed_negative_particles(block)){
            Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*removed_negative_particles(block),2*side-3);
            if(!removed_negative_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(removed_negative_particles(block));}
        if(use_removed_positive_particles)if(removed_positive_particles(block)){
            Delete_Particles_Outside_Grid(domain_boundaries[side-1][axis],axis,*removed_positive_particles(block),2*side-3);
            if(!removed_positive_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(removed_positive_particles(block));}}
}
//#####################################################################
// Function Delete_Particles_In_Local_Maximum_Phi_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_In_Local_Maximum_Phi_Cells(const int sign)
{
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    T tolerance=levelset.small_number*levelset.grid.Minimum_Edge_Length();
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(levelset.grid.Domain_Indices(1),thread_queue).template Run<int,T>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Particles_In_Local_Maximum_Phi_Cells_Threaded,sign,tolerance);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV> threaded_iterator(domain,thread_queue);
    threaded_iterator.template Run<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Marked_Particles,positive_particles);
    threaded_iterator.template Run<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Marked_Particles,negative_particles);
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
}
//#####################################################################
// Function Delete_Particles_In_Local_Maximum_Phi_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_In_Local_Maximum_Phi_Cells_Threaded(RANGE<TV_INT>& domain,const int sign,const T tolerance)
{
    for(CELL_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        bool local_minima=true;
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(sign*levelset.phi(cell)<sign*levelset.phi(cell+axis_vector)-tolerance||sign*levelset.phi(cell)<sign*levelset.phi(cell-axis_vector)-tolerance){local_minima=false;break;}}
        if(!local_minima)continue;
        TV_INT blocks[T_GRID::number_of_nodes_per_cell];levelset.grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell,blocks);
        for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++){TV location=iterator.Location();
            if(sign==-1&&negative_particles(blocks[i])){PARTICLE_LEVELSET_PARTICLES<TV>* particles=negative_particles(blocks[i]);
                while(particles){for(int k=1;k<=particles->array_collection->Size();k++) if(sqr(particles->radius(k))>=(particles->X(k)-location).Magnitude_Squared()){particles->radius(k)=-1;}particles=particles->next;}}
            else if(sign==1&&positive_particles(blocks[i])){PARTICLE_LEVELSET_PARTICLES<TV>* particles=positive_particles(blocks[i]);
                while(particles){for(int k=1;k<=particles->array_collection->Size();k++) if(sqr(particles->radius(k))>=(particles->X(k)-location).Magnitude_Squared()){particles->radius(k)=-1;}particles=particles->next;}}}}
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int side)
{
    PARTICLE_LEVELSET_PARTICLES<TV>* local_particles=&particles;
    while(local_particles){for(int k=1;k<=local_particles->array_collection->Size();k++) if(side*local_particles->X(k)[axis]>side*domain_boundary){
        Delete_Particle_And_Clean_Memory(&particles,*local_particles,k);k--;}
    local_particles=local_particles->next;} // side should be -1 or 1
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int side)
{
    for(int k=particles.array_collection->Size();k>=1;k--) if(side*particles.X(k)[axis]>side*domain_boundary) Delete_Particle_And_Clean_Memory(&particles,particles,k);
}
//#####################################################################
// Function Identify_And_Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Identify_And_Remove_Escaped_Particles(const T_FACE_ARRAYS_SCALAR& V,const T radius_fraction,const T time,const bool verbose)
{
    escaped_positive_particles.Resize(levelset.grid.Domain_Indices(3));
    escaped_negative_particles.Resize(levelset.grid.Domain_Indices(3));
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(domain,thread_queue).template Run<const T_FACE_ARRAYS_SCALAR&,T,T,bool>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Identify_And_Remove_Escaped_Particles_Threaded,V,radius_fraction,time,verbose);
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
}
//#####################################################################
// Function Identify_And_Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Identify_And_Remove_Escaped_Particles_Threaded(RANGE<TV_INT>& domain,const T_FACE_ARRAYS_SCALAR& V,const T radius_fraction,const T time,const bool verbose)
{
    T_FACE_LOOKUP V_lookup(V);
    T_FACE_LOOKUP_COLLIDABLE V_lookup_collidable(V_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP V_lookup_collidable_lookup(V_lookup_collidable,V_lookup);

    Identify_Escaped_Particles(domain);int p=0,n=0;
    if(use_removed_positive_particles){
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(positive_particles(block)){
            p+=Remove_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*positive_particles(block),escaped_positive_particles(block),1,removed_positive_particles(block),
                V_lookup_collidable_lookup,radius_fraction,time);
            if(!positive_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(positive_particles(block));}}}
    else{
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(positive_particles(block)){
            p+=Remove_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*positive_particles(block),escaped_positive_particles(block),1,radius_fraction);
            if(!positive_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(positive_particles(block));}}}
    if(use_removed_negative_particles){
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(negative_particles(block)){
            n+=Remove_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*negative_particles(block),escaped_negative_particles(block),-1,removed_negative_particles(block),
                V_lookup_collidable_lookup,radius_fraction,time);
            if(!negative_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(negative_particles(block));}}}
    else{
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(negative_particles(block)){
            n+=Remove_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*negative_particles(block),escaped_negative_particles(block),-1,radius_fraction);
            if(!negative_particles(block)->array_collection->Size()) Free_Particle_And_Clear_Pointer(negative_particles(block));}}}
    //LOG::cout << "new removed positive particles = " << p << std::endl;LOG::cout << "new removed negative particles = " << n << std::endl;
}
//#####################################################################
// Function Identify_And_Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Identify_And_Remove_Escaped_Positive_Particles(const T_FACE_ARRAYS_SCALAR& V,const T radius_fraction,const T time,const bool verbose)
{
    T_FACE_LOOKUP V_lookup(V);
    T_FACE_LOOKUP_COLLIDABLE V_lookup_collidable(V_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP V_lookup_collidable_lookup(V_lookup_collidable,V_lookup);

    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    escaped_positive_particles.Resize(levelset.grid.Domain_Indices(3));
    escaped_negative_particles.Resize(levelset.grid.Domain_Indices(3));
    Identify_Escaped_Particles(domain);int p=0;
    if(use_removed_positive_particles){
        for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(positive_particles(block))
            p+=Remove_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*positive_particles(block),escaped_positive_particles(block),1,removed_positive_particles(block),
                V_lookup_collidable_lookup,radius_fraction,time);}}
    else{
        for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(positive_particles(block))
            p+=Remove_Escaped_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),*positive_particles(block),escaped_positive_particles(block),1,radius_fraction);}}
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> template<class T_FACE_LOOKUP_LOOKUP> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Remove_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles,const T_FACE_LOOKUP_LOOKUP& V,const T radius_fraction,const T time)
{
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block):false;if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    int removed=0;T one_over_radius_multiplier=-sign/radius_fraction;
    // TODO: limit the amount of mass removed - don't let the particles just set their own radii
    ARRAY_VIEW<int>* id=save_removed_particle_times && store_unique_particle_id?particles.array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID):0;
    ARRAY<PAIR<int,T> > removed_particle_times_local;
    ARRAY<bool> local_escaped(escaped);
    PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=&particles;
    while(cell_particles){for(int k=cell_particles->array_collection->Size();k>=1;k--) if(local_escaped(k) && one_over_radius_multiplier*levelset.Phi(cell_particles->X(k))>cell_particles->radius(k)){
        if(!removed_particles) removed_particles=Allocate_Particles(template_removed_particles);
        removed++;
        if(id) removed_particle_times_local.Append(PAIR<int,T>((*id)(k),time));
        Copy_Particle(*cell_particles,*removed_particles,k);
        Delete_Particle_And_Clean_Memory(&particles,*cell_particles,k);
        removed_particles->V.Last()=linear_interpolation_collidable.Clamped_To_Array_Face(levelset.grid,V,removed_particles->X.Last());}
        cell_particles=cell_particles->next;
        if(cell_particles) Identify_Escaped_Particles(block,*cell_particles,local_escaped,sign);}
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
#ifdef USE_PTHREADS
    if(thread_queue) pthread_mutex_lock(&cell_lock);
#endif
    removed_particle_times.Append_Elements(removed_particle_times_local);
#ifdef USE_PTHREADS
    if(thread_queue) pthread_mutex_unlock(&cell_lock);
#endif
    return removed;
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> bool PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Fix_Momentum_With_Escaped_Particles(const T_FACE_ARRAYS_SCALAR& V,const T_ARRAYS_SCALAR& momentum_lost,const T radius_fraction,const T mass_scaling,const T time,const bool force)
{
    T_FACE_LOOKUP V_lookup(V);
    T_FACE_LOOKUP_COLLIDABLE V_lookup_collidable(V_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP V_lookup_collidable_lookup(V_lookup_collidable,V_lookup);
    
    bool done=true;
    T_LINEAR_INTERPOLATION_SCALAR linear_interpolation;
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(positive_particles(block)){
        T local_momentum_lost=linear_interpolation.Clamped_To_Array(levelset.grid,momentum_lost,iterator.Location());
        if(local_momentum_lost==0) continue; //Nothing to add
        BLOCK_UNIFORM<T_GRID> block_uniform(levelset.grid,block);
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles=removed_negative_particles(block);
        if(!removed_particles){
            if(!force){done=false;continue;}
            else removed_particles=template_removed_particles.Clone();}
        T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
        bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block_uniform):false;if(near_objects) levelset.Enable_Collision_Aware_Interpolation(-1);
        for(int k=removed_particles->array_collection->Size();k>=1;k--){T fraction=local_momentum_lost/removed_particles->array_collection->Size();
            removed_particles->V(k)+=fraction/mass_scaling;}
        if(!removed_particles->array_collection->Size() && force) Add_Negative_Particle(iterator.Location(),local_momentum_lost/mass_scaling*linear_interpolation_collidable.Clamped_To_Array_Face(levelset.grid,V_lookup_collidable_lookup,iterator.Location()).Normalized(),(unsigned short)(random.Get_Number()*USHRT_MAX));
        if(!removed_particles->array_collection->Size()) done=false;
        if(near_objects) levelset.Disable_Collision_Aware_Interpolation();}}
    PHYSBAM_ASSERT(!force||done);
    return done;
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> bool PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Fix_Momentum_With_Escaped_Particles(const TV& location,const T_FACE_ARRAYS_SCALAR& V,const T momentum_lost,const T radius_fraction,const T mass_scaling,const T time,const bool force)
{
    T_FACE_LOOKUP V_lookup(V);
    T_FACE_LOOKUP_COLLIDABLE V_lookup_collidable(V_lookup,*levelset.collision_body_list,levelset.face_velocities_valid_mask_current);
    typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP V_lookup_collidable_lookup(V_lookup_collidable,V_lookup);
    T_LINEAR_INTERPOLATION_SCALAR linear_interpolation;
    T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR linear_interpolation_collidable;
    TV_INT block=levelset.grid.Closest_Node(location);
    
    BLOCK_UNIFORM<T_GRID> block_uniform(levelset.grid,block);
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles=removed_negative_particles(block);
    if(!removed_particles){
        if(!force) return false;
        else removed_particles=template_removed_particles.Clone();}
    bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block_uniform):false;if(near_objects) levelset.Enable_Collision_Aware_Interpolation(-1);
    for(int k=removed_particles->array_collection->Size();k>=1;k--){T fraction=momentum_lost/removed_particles->array_collection->Size();
        removed_particles->V(k)+=fraction/mass_scaling;}
    if(!removed_particles->array_collection->Size() && force) Add_Negative_Particle(location,momentum_lost/mass_scaling*linear_interpolation_collidable.Clamped_To_Array_Face(levelset.grid,V_lookup_collidable_lookup,location).Normalized(),(unsigned short)(random.Get_Number()*USHRT_MAX));
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
    PHYSBAM_ASSERT(!force||removed_particles->array_collection->Size());
    return (removed_particles->array_collection->Size())==0?false:true;
}
//#####################################################################
// Function Remove_Escaped_Particles
//#####################################################################
template<class T_GRID> int PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Remove_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,const T radius_fraction)
{
    bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block):false;if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    int deleted=0;T one_over_radius_multiplier=-sign/radius_fraction;
    PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=&particles;
    while(cell_particles){for(int k=cell_particles->array_collection->Size();k>=1;k--) if(escaped(k) && one_over_radius_multiplier*levelset.Phi(cell_particles->X(k))>cell_particles->radius(k)){
        Delete_Particle_And_Clean_Memory(&particles,*cell_particles,k);deleted++;k--;}
        cell_particles=cell_particles->next;}
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
    return deleted;
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reincorporate_Removed_Particles(const T radius_fraction,const T mass_scaling,T_FACE_ARRAYS_SCALAR* V,const bool conserve_momentum_for_removed_negative_particles)
{
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(domain,thread_queue).template Run<T,T,T_FACE_ARRAYS_SCALAR*>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Reincorporate_Removed_Particles_Threaded,radius_fraction,mass_scaling,conserve_momentum_for_removed_negative_particles?V:0);
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reincorporate_Removed_Particles_Threaded(RANGE<TV_INT>& domain,const T radius_fraction,const T mass_scaling,T_FACE_ARRAYS_SCALAR* V)
{
    if(use_removed_positive_particles){
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(removed_positive_particles(block))
            Reincorporate_Removed_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),positive_particles(block),1,*removed_positive_particles(block),radius_fraction,mass_scaling,0);}}
    if(use_removed_negative_particles){
        for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();if(removed_negative_particles(block))
            Reincorporate_Removed_Particles(BLOCK_UNIFORM<T_GRID>(levelset.grid,block),negative_particles(block),-1,*removed_negative_particles(block),radius_fraction,mass_scaling,V);}}
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Reincorporate_Removed_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>*& particles,const int sign,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& removed_particles,
    const T radius_fraction,const T mass_scaling,T_FACE_ARRAYS_SCALAR* V)
{
    bool near_objects=levelset.collision_body_list?levelset.collision_body_list->Occupied_Block(block):false;if(near_objects) levelset.Enable_Collision_Aware_Interpolation(sign);
    T one_over_radius_multiplier=-sign/radius_fraction;
    T half_unit_sphere_size_over_cell_size=(T).5*(T)unit_sphere_size<TV::dimension>::value/levelset.grid.Cell_Size();
    ARRAY_VIEW<T>* material_volume=vof_advection?removed_particles.array_collection->template Get_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME):0;
    for(int k=removed_particles.array_collection->Size();k>=1;k--) if(reincorporate_removed_particles_everywhere || (one_over_radius_multiplier*levelset.Phi(removed_particles.X(k))<removed_particles.radius(k))){
        if(!particles) particles=Allocate_Particles(template_particles);
        // rasterize the velocity onto the velocity field
        if(V){
            TV_INT cell_index=levelset.grid.Cell(removed_particles.X(k),0);
            T r=removed_particles.radius(k);
            TV half_impulse;
            if(material_volume) half_impulse=removed_particles.V(k)*mass_scaling*(T).5*(*material_volume)(k);
            else half_impulse = removed_particles.V(k)*mass_scaling*pow<TV::dimension>(r)*half_unit_sphere_size_over_cell_size;
            for(int i=1;i<=T_GRID::dimension;++i){
                TV_INT face_index(cell_index);
                typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE &face_velocity=V->Component(i);
#ifdef USE_PTHREADS
                if(thread_queue) pthread_mutex_lock(&cell_lock); //TODO (mlentine) Remove this with thread_safe code
#endif                
                if(face_velocity.Valid_Index(face_index)) face_velocity(face_index)+=half_impulse[i];
                face_index[i]+=1;
                if(face_velocity.Valid_Index(face_index)) face_velocity(face_index)+=half_impulse[i];
#ifdef USE_PTHREADS
                if(thread_queue) pthread_mutex_unlock(&cell_lock);
#endif
                }}
        removed_particles.radius(k)=clamp(removed_particles.radius(k),minimum_particle_radius,maximum_particle_radius);
        Move_Particle(removed_particles,*particles,k);}// maybe recalculate radius to get rid of raining effect?
    if(near_objects) levelset.Disable_Collision_Aware_Interpolation();
}
//#####################################################################
// Function Delete_Particles_Far_From_Interface
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Far_From_Interface(const int discrete_band)
{
    RANGE<TV_INT> domain(levelset.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    assert(discrete_band<8);
    T_ARRAYS_CHAR near_interface(levelset.grid.Cell_Indices(3));
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(levelset.grid.Domain_Indices(3),thread_queue).template Run<T_ARRAYS_CHAR&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Particles_Far_From_Interface_Part_One,near_interface);
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV> part_two_iterator(levelset.grid.Domain_Indices(3),thread_queue);
    for(int distance=1;distance<=discrete_band;distance++){
        part_two_iterator.template Run<T_ARRAYS_CHAR&,int>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Particles_Far_From_Interface_Part_Two,near_interface,distance);}
    DOMAIN_ITERATOR_THREADED_ALPHA<PARTICLE_LEVELSET_UNIFORM<T_GRID>,TV>(domain,thread_queue).template Run<T_ARRAYS_CHAR&>(*this,&PARTICLE_LEVELSET_UNIFORM<T_GRID>::Delete_Particles_Far_From_Interface_Part_Three,near_interface);
    Consistency_Check(domain,positive_particles);
    Consistency_Check(domain,negative_particles);
}
//#####################################################################
// Function Delete_Particles_Far_From_Interface
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Far_From_Interface_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface)
{
    const T_ARRAYS_BOOL_DIMENSION& cell_neighbors_visible=levelset.collision_body_list->cell_neighbors_visible;
    for(CELL_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT neighbor=cell+TV_INT::Axis_Vector(axis);
            if(near_interface.Valid_Index(neighbor) && cell_neighbors_visible(cell)(axis) && LEVELSET_UTILITIES<T>::Interface(levelset.phi(cell),levelset.phi(neighbor))){
            near_interface(cell)=near_interface(neighbor)=1;}}}
}
//#####################################################################
// Function Delete_Particles_Far_From_Interface
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Far_From_Interface_Part_Two(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface,const int distance)
{
    const T_ARRAYS_BOOL_DIMENSION& cell_neighbors_visible=levelset.collision_body_list->cell_neighbors_visible;
    RANGE<TV_INT> ghost_domain(domain);if(ghost_domain.min_corner.x>0) ghost_domain.min_corner.x-=1;
    char new_mask=1<<distance,old_mask=new_mask-1;
    for(CELL_ITERATOR iterator(levelset.grid,ghost_domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT neighbor=cell+TV_INT::Axis_Vector(axis);
            if(near_interface.Valid_Index(neighbor) && cell_neighbors_visible(cell)(axis) && (near_interface(cell)|near_interface(neighbor))&old_mask){
                if(domain.Lazy_Inside(cell)) near_interface(cell)|=new_mask;if(domain.Lazy_Inside(neighbor)) near_interface(neighbor)|=new_mask;}}}
}
//#####################################################################
// Function Delete_Particles_Far_From_Interface
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Delete_Particles_Far_From_Interface_Part_Three(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface)
{
    for(NODE_ITERATOR iterator(levelset.grid,domain);iterator.Valid();iterator.Next()){TV_INT b=iterator.Node_Index();if(negative_particles(b) || positive_particles(b)){
        T_BLOCK block(levelset.grid,b);
        for(int i=0;i<T_GRID::number_of_cells_per_block;i++)if(near_interface.Valid_Index(block.Cell(i))&&near_interface(block.Cell(i))) goto NEAR_INTERFACE;
        // if not near interface, delete particles
        Free_Particle_And_Clear_Pointer(negative_particles(b));
        Free_Particle_And_Clear_Pointer(positive_particles(b));
        NEAR_INTERFACE:;}}
}
//#####################################################################
// Function Add_Negative_Particle
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Add_Negative_Particle(const TV& location,const TV& particle_velocity,const unsigned short quantized_collision_distance)
{
    //PHYSBAM_NOT_IMPLEMENTED(); // TODO: figure out whether/how this should be collision aware
    // TODO: still need to figure out whether this should be collision aware, but for now I need it

    TV_INT block=levelset.grid.Block_Index(location,3);
    if(!negative_particles.Valid_Index(block)) return;
    if(levelset.Phi(location)<minimum_particle_radius){
        if(!negative_particles(block)) negative_particles(block)=Allocate_Particles(template_particles);
        PARTICLE_LEVELSET_PARTICLES<TV>* cell_particles=negative_particles(block);
        int index=Add_Particle(cell_particles);
        cell_particles->X(index)=location;negative_particles(block)->radius(index)=minimum_particle_radius;
        cell_particles->quantized_collision_distance(index)=quantized_collision_distance;
        if(store_unique_particle_id){
            ARRAY_VIEW<int>* id=cell_particles->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);
            PHYSBAM_ASSERT(id);
#ifndef WIN32
            (*id)(index)=__exchange_and_add(&last_unique_particle_id,1);}
#else
            PHYSBAM_ASSERT(!thread_queue);
            (*id)(index)=++last_unique_particle_id;}
#endif
        if(vof_advection){
            ARRAY_VIEW<T>* material_volume=cell_particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME);
            PHYSBAM_ASSERT(material_volume);
            (*material_volume)(index)=0;}}
    else{
        if(!removed_negative_particles(block)) removed_negative_particles(block)=Allocate_Particles(template_removed_particles);
        int index=Add_Particle(removed_negative_particles(block));
        removed_negative_particles(block)->X(index)=location;removed_negative_particles(block)->radius(index)=minimum_particle_radius;
        removed_negative_particles(block)->V(index)=particle_velocity;removed_negative_particles(block)->quantized_collision_distance(index)=quantized_collision_distance;
        if(store_unique_particle_id){
            ARRAY_VIEW<int>* id=removed_negative_particles(block)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);
            PHYSBAM_ASSERT(id);
#ifndef WIN32
            (*id)(index)=__exchange_and_add(&last_unique_particle_id,1);}}
#else
            PHYSBAM_ASSERT(!thread_queue);
            (*id)(index)=last_unique_particle_id++;}}
#endif
}
//#####################################################################
// Function Exchange_Overlap_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Exchange_Overlap_Particles()
{
    if(mpi_grid){
        Exchange_Overlapping_Block_Particles(*mpi_grid,template_particles,negative_particles,(int)(cfl_number+(T)3.5),*this);
        Exchange_Overlapping_Block_Particles(*mpi_grid,template_particles,positive_particles,(int)(cfl_number+(T)3.5),*this);}
}
//#####################################################################
// Function Copy_From_Move_List_Threaded
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_UNIFORM<T_GRID>::
Copy_From_Move_List_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,ARRAY<ARRAY<TRIPLE<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT,int> >,TV_INT>& move_particles,const ARRAY<int,TV_INT>* domain_index)
{
    PHYSBAM_NOT_IMPLEMENTED();
}

//#####################################################################
#define INSTANTIATION_HELPER(T_GRID) \
    template class PARTICLE_LEVELSET_UNIFORM<T_GRID >;                  \
    template void PARTICLE_LEVELSET_UNIFORM<T_GRID >::Euler_Step_Particles_Wrapper(const GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& V,PARTICLE_LEVELSET_UNIFORM<T_GRID >::T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt, \
        const T time,const bool update_particle_cells_after_euler_step,const bool assume_particles_in_correct_blocks,const bool enforce_domain_boundaries); \
    template void PARTICLE_LEVELSET_UNIFORM<T_GRID >::Euler_Step_Particles_Wrapper(const GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& V,PARTICLE_LEVELSET_UNIFORM<T_GRID >::T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt, \
        const T time,const bool update_particle_cells_after_euler_step,const bool assume_particles_in_correct_blocks,const bool enforce_domain_boundaries);

#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >));
#endif
