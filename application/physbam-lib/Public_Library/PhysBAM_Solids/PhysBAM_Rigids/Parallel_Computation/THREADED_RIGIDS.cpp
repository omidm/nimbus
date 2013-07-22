//#####################################################################
// Copyright 2011, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/THREADED_RIGIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>

using namespace PhysBAM;

//#####################################################################
// Function THREADED_RIGIDS
//#####################################################################
template<class TV> THREADED_RIGIDS<TV>::
THREADED_RIGIDS(ARRAY<THREAD_PACKAGE>& buffers_input,const int tid_input,const int number_of_threads)
    :buffers(buffers_input),tid(tid_input),number_of_processors(number_of_threads)
{
    LOG::SCOPE scope("THREADED INITIALIZE","Initializing Threading",tid);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    std::stringstream ss;ss<<"number of processes = "<<number_of_processors<<std::endl;LOG::filecout(ss.str(),tid);
#endif

#ifdef USE_PTHREADS
    // initialize pthread mutex and barrier
    if(tid==1){
        lock=new pthread_mutex_t;
        pthread_mutex_init(lock,NULL);
        THREAD_PACKAGE pack(sizeof(pthread_mutex_t*));
        *(pthread_mutex_t**)(&pack.buffer(1))=lock;
        buffers.Append(pack);}
    else{
        while(buffers.m==0) PROCESS_UTILITIES::PB_Sleep(1);
        lock=*(pthread_mutex_t**)(&buffers(1).buffer(1));}
    if(tid==1){
        barr=new pthread_barrier_t;
        pthread_barrier_init(barr,NULL,number_of_processors);
        THREAD_PACKAGE pack(sizeof(pthread_barrier_t*));
        *(pthread_barrier_t**)(&pack.buffer(1))=barr;
        buffers.Append(pack);}
    else{
        while(buffers.m<=1) PROCESS_UTILITIES::PB_Sleep(1);
        barr=*(pthread_barrier_t**)(&buffers(2).buffer(1));}
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}    
//#####################################################################
// Function ~THREADED_RIGIDS
//#####################################################################
template<class TV> THREADED_RIGIDS<TV>::
~THREADED_RIGIDS()
{
#ifdef USE_PTHREADS
    // initialize pthread mutex and barrier
    if(tid==1){
        delete lock;
        lock=NULL;
        THREAD_PACKAGE pack(sizeof(pthread_mutex_t*));
        *(pthread_mutex_t**)(&pack.buffer(1))=lock;
        buffers.Append(pack);}
    else{
        while(buffers.m==0) PROCESS_UTILITIES::PB_Sleep(1);
        lock=*(pthread_mutex_t**)(&buffers(1).buffer(1));}
    if(tid==1){
        delete barr;
        barr=NULL;
        THREAD_PACKAGE pack(sizeof(pthread_barrier_t*));
        *(pthread_barrier_t**)(&pack.buffer(1))=barr;
        buffers.Append(pack);}
    else{
        while(buffers.m<=1) PROCESS_UTILITIES::PB_Sleep(1);
        barr=*(pthread_barrier_t**)(&buffers(2).buffer(1));}
    if(tid==1) buffers.m=0;
#endif
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Initialize(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,TV_INT processes_per_dimension_input,
    COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition)
{
    processes_per_dimension=processes_per_dimension_input;
        
    if(!rigid_body_collection_input.rigid_body_particle.store_id){
        if(ATTRIBUTE_INDEX index=rigid_body_collection_input.rigid_body_particle.array_collection->Get_Attribute_Index(ATTRIBUTE_ID_ID)){
            rigid_body_collection_input.rigid_body_particle.store_id=true;
            ARRAY_VIEW<int>* id=rigid_body_collection_input.rigid_body_particle.array_collection->template Get_Array_From_Index<int>(index);
            rigid_body_collection_input.rigid_body_particle.id.base_pointer=id->base_pointer;
            rigid_body_collection_input.rigid_body_particle.id.m=id->m;}
        else{
            rigid_body_collection_input.rigid_body_particle.Store_Id();
            for(int i=1;i<=rigid_body_collection_input.rigid_body_particle.array_collection->Size();i++){
                rigid_body_collection_input.rigid_body_particle.id(i)=i;}}}
    
    global_domain=spatial_partition.Scene_Bounding_Box();
    local_domain=Split_Range(processes_per_dimension,global_domain);
    Simple_Partition(rigid_body_collection_input,spatial_partition);
}
//#####################################################################
// Function Broadcast_Positions
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Broadcast_Positions(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Update_Partitions
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Update_Partitions(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
    COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition)
{
#ifdef USE_PTHREADS
    RANGE<TV> new_global_domain=spatial_partition.Scene_Bounding_Box();
    if(new_global_domain!=global_domain){
        global_domain=new_global_domain;
        local_domain=Split_Range(processes_per_dimension,global_domain);}
    
    ARRAY<ARRAY<int> > bodies_to_send(number_of_processors);
    TV domain_size=local_domain.Edge_Lengths();
    for(int body=1;body<=particles_of_partition(PARTITION_ID(tid)).m;body++){
        int body_id=particles_of_partition(PARTITION_ID(tid))(body);
        if(!local_domain.Lazy_Inside(rigid_body_collection_input.rigid_body_particle.X(body_id))){
            TV_INT partition_coordinate;
            for(int d=1;d<=TV::dimension;d++)
                partition_coordinate(d)=(int)ceil((rigid_body_collection_input.rigid_body_particle.X(body_id)(d)-global_domain.min_corner(d))/domain_size(d));
            bodies_to_send(partition_coordinates_to_rank(partition_coordinate)).Append(body_id);}}
    
    // update particles of partition
    int buffer_size=1; // 1 for the size of the array
    for(int i=1;i<=bodies_to_send.m;i++) buffer_size+=(bodies_to_send(i).m+1);
    buffer_size+=1; // extra padding
    THREAD_PACKAGE pack(buffer_size*sizeof(int));pack.send_tid=tid;pack.recv_tid=0;int position=0;
    *(int*)(&pack.buffer(position+1))=bodies_to_send.m;position+=sizeof(int); // pack size
    for(int i=1;i<=bodies_to_send.m;i++){
        *(int*)(&pack.buffer(position+1))=bodies_to_send(i).m;position+=sizeof(int); // pack size
        for(int j=1;j<=bodies_to_send(i).m;j++){
            *(int*)(&pack.buffer(position+1))=bodies_to_send(i)(j);
            position+=sizeof(int);}}

    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);

    pthread_barrier_wait(barr);    

    for(int i=1;i<=buffers.m;i++){
        if(buffers(i).send_tid!=tid){
            int position=0;
            int proc_count=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
            for(int new_proc=1;new_proc<=proc_count;new_proc++){
                int body_count=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
                for(int body=1;body<=body_count;body++){
                    int body_id=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
                    partition_id_from_particle_index(body_id)=PARTITION_ID(new_proc);
                    for(int b=1;b<=particles_of_partition(PARTITION_ID(buffers(i).send_tid)).m;b++){ // TODO faster way to do this?
                        if(particles_of_partition(PARTITION_ID(buffers(i).send_tid))(b)==body_id){
                            particles_of_partition(PARTITION_ID(buffers(i).send_tid)).Remove_Index_Lazy(b);
                            break;}}
                    particles_of_partition(PARTITION_ID(new_proc)).Append(body_id);}}}
        else{
            for(int new_proc=1;new_proc<=bodies_to_send.m;new_proc++){
                for(int body=1;body<=bodies_to_send(new_proc).m;body++){
                    partition_id_from_particle_index(bodies_to_send(new_proc)(body))=PARTITION_ID(new_proc);
                    for(int b=1;b<=particles_of_partition(PARTITION_ID(tid)).m;b++){ // TODO faster way to do this?
                        if(particles_of_partition(PARTITION_ID(tid))(b)==bodies_to_send(new_proc)(body)){
                            particles_of_partition(PARTITION_ID(tid)).Remove_Index_Lazy(b);
                            break;}}
                    particles_of_partition(PARTITION_ID(new_proc)).Append(bodies_to_send(new_proc)(body));}}}}

    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Split_Range
//#####################################################################
template<class TV> RANGE<TV> THREADED_RIGIDS<TV>::
Split_Range(TV_INT& processes_per_dimension,const RANGE<TV>& global_range)
{
    TV range=global_range.max_corner-global_range.min_corner;
    processes_per_dimension=Split_Range(range,processes_per_dimension);
    GRID<TV> process_grid(processes_per_dimension,RANGE<TV>::Centered_Box());
    ARRAY<TV_INT> all_coordinates;all_coordinates.Resize(number_of_processors);
    partition_coordinates_to_rank.Resize(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),processes_per_dimension));
    int count=0;
    for(typename GRID<TV>::NODE_ITERATOR iterator(process_grid);iterator.Valid();iterator.Next()){
        all_coordinates(++count)=iterator.Node_Index();
        partition_coordinates_to_rank(iterator.Node_Index())=count;}
    TV_INT coordinates=all_coordinates(tid);
    TV start,end;
    for(int axis=1;axis<=TV::dimension;axis++){
        start[axis]=boundaries(axis)(coordinates[axis]);
        end[axis]=boundaries(axis)(coordinates[axis]+1);}
    return RANGE<TV>(start+global_range.min_corner,end+global_range.min_corner);
}
template<class T>
static bool Minimize_2D_Surface_Area(const int number_of_processes,const T x,const T y,int& count_x)
{
    count_x=max(1,(int)sqrt((T)number_of_processes*x/y));
    if(number_of_processes%count_x==0){
        if(number_of_processes%(count_x+1)==0 && y*count_x+x*number_of_processes/count_x >= y*(count_x+1)+x*number_of_processes/(count_x+1)) count_x++;}
    else if(number_of_processes%(count_x+1)==0) count_x++;
    else return false;
    return true;
}
template<class TV> VECTOR<int,1> THREADED_RIGIDS<TV>::
Split_Range(const VECTOR<T,1>& global_range,const VECTOR<int,1>& processes_per_dimension)
{
    PHYSBAM_ASSERT(!processes_per_dimension.x || processes_per_dimension.x==number_of_processors);
    boundaries.Resize(1);
    Split_Dimension(global_range.x,number_of_processors,boundaries(1));
    return VECTOR<int,1>(number_of_processors);
}
template<class TV> VECTOR<int,2> THREADED_RIGIDS<TV>::
Split_Range(const VECTOR<T,2>& global_range,const VECTOR<int,2>& processes_per_dimension)
{
    T x=global_range.x,y=global_range.y;VECTOR<int,2> count;
    if(processes_per_dimension!=VECTOR<int,2>()) count=processes_per_dimension;
    else{ // try to figure out counts by minimizing surface area between processes
        if(!Minimize_2D_Surface_Area<T>(number_of_processors,x,y,count.x)){
            LOG::cerr<<"Don't know how to divide domain in both directions."<<std::endl;PHYSBAM_NOT_IMPLEMENTED();}
        count.y=number_of_processors/count.x;}
    PHYSBAM_ASSERT(count.x*count.y==number_of_processors);
    std::stringstream ss;ss<<"dividing domain into "<<count.x<<" by "<<count.y<<" processor grid"<<std::endl;LOG::filecout(ss.str(),tid);
    boundaries.Resize(2);
    Split_Dimension(x,count.x,boundaries(1));
    Split_Dimension(y,count.y,boundaries(2));
    return count;
}
template<class TV> VECTOR<int,3> THREADED_RIGIDS<TV>::
Split_Range(const VECTOR<T,3>& global_range,const VECTOR<int,3>& processes_per_dimension)
{
    T x=global_range.x,y=global_range.y,z=global_range.z;VECTOR<int,3> count;
    if(processes_per_dimension!=VECTOR<int,3>()) count=processes_per_dimension;
    else{ // try to figure out counts by minimizing surface area between processes
        T minimum_surface_area=FLT_MAX;VECTOR<int,3> test_count;
        for(test_count.z=1;test_count.z<=number_of_processors;test_count.z++) if(number_of_processors%test_count.z==0){
                if(Minimize_2D_Surface_Area<T>(number_of_processors/test_count.z,x,y,test_count.x)){
                    test_count.y=number_of_processors/(test_count.x*test_count.z);
                    T surface_area=test_count.x*(y*z)+test_count.y*(x*z)+test_count.z*(x*y);
                    if(surface_area<minimum_surface_area){count=test_count;minimum_surface_area=surface_area;}}
                if(Minimize_2D_Surface_Area<T>(number_of_processors/test_count.z,x,y,test_count.y)){
                    test_count.x=number_of_processors/(test_count.y*test_count.z);
                    T surface_area=test_count.x*(y*z)+test_count.y*(x*z)+test_count.z*(x*y);
                    if(surface_area<minimum_surface_area){count=test_count;minimum_surface_area=surface_area;}}}
        if(minimum_surface_area==INT_MAX){LOG::cerr<<"Don't know how to divide domain in all directions."<<std::endl;PHYSBAM_NOT_IMPLEMENTED();}}
    PHYSBAM_ASSERT(count.x*count.y*count.z==number_of_processors);
    std::stringstream ss;ss<<"dividing domain into "<<count.x<<" by "<<count.y<<" by "<<count.z<<" processor grid"<<std::endl;LOG::filecout(ss.str(),tid);
    boundaries.Resize(3);
    Split_Dimension(x,count.x,boundaries(1));
    Split_Dimension(y,count.y,boundaries(2));
    Split_Dimension(z,count.z,boundaries(3));
    return count;
}
//#####################################################################
// Function Split_Dimension
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Split_Dimension(const T x,const int processes,ARRAY<T>& boundaries)
{
    T range_over_processes=x/processes;
    boundaries.Resize(processes+1);boundaries(1)=0;
    for(int p=1;p<=processes;p++)boundaries(p+1)=boundaries(p)+range_over_processes;
}
//#####################################################################
// Function Simple_Partition
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Simple_Partition(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
    COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition)
{
    particles_of_partition.Resize(PARTITION_ID(number_of_processors+1)); // Last index is for global particles
    for(PARTITION_ID id(1);id<=particles_of_partition.Size();id++) particles_of_partition(id).Resize(0);
    partition_id_from_particle_index.Remove_All();
    partition_id_from_particle_index.Resize(rigid_body_collection_input.rigid_body_particle.array_collection->Size());
    
    OPERATION_HASH<COLLISION_GEOMETRY_ID> already_added;
    ARRAY<COLLISION_GEOMETRY_ID> object_indices;spatial_partition.Get_Potential_Collisions(local_domain,object_indices,already_added);
    // Determine which bodies belong on which processors
    PARTITION_ID id(tid);
    for(int i=1;i<=partition_id_from_particle_index.m;i++) partition_id_from_particle_index(i)=PARTITION_ID(-1);
    for(int i=1;i<=rigid_body_collection_input.rigid_body_particle.array_collection->Size();i++){
        int p=rigid_body_collection_input.rigid_body_particle.id(i);
        if(local_domain.Lazy_Inside(rigid_body_collection_input.rigid_body_particle.X(p)))
            partition_id_from_particle_index(p)=id;}
    
    //Add bodies not in partition to final extra partition
    for(int i=1;i<=spatial_partition.bodies_not_in_partition.m;i++){
        int p=rigid_body_collection_input.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(spatial_partition.bodies_not_in_partition(i));
        partition_id_from_particle_index(p)=particles_of_partition.Size();}
    
    for(int i=1;i<=partition_id_from_particle_index.m;i++){
        partition_id_from_particle_index(i)=PARTITION_ID(Reduce_Max(Value(partition_id_from_particle_index(i))));
        particles_of_partition(partition_id_from_particle_index(i)).Append(i);}
}
//#####################################################################
// Function Clear_Impulse_Accumulators
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Clear_Impulse_Accumulators(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
{
    for(int j=1;j<=rigid_body_collection_input.rigid_body_particle.array_collection->Size();j++) rigid_body_collection_input.Rigid_Body(j).impulse_accumulator->Reset();
}
//#####################################################################
// Function Exchange_All_Impulses
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Exchange_All_Impulses(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TWIST<TV> >& mpi_rigid_velocity_save,ARRAY<T_SPIN>& mpi_rigid_angular_momentum_save,
    RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const bool euler_step_with_new_velocity,const T dt,const T time)
{
    ARRAY<bool> need_to_reevolve(rigid_body_collection.rigid_body_particle.array_collection->Size());
    Prune_And_Exchange_Impulses(rigid_body_collection,need_to_reevolve);
    
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        if(need_to_reevolve(i)){
            int p=rigid_body_collection.rigid_body_particle.id(i);
            RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
            // Restore Velocities
            rigid_body.V()=mpi_rigid_velocity_save(p).linear;
            rigid_body.Angular_Momentum()=mpi_rigid_angular_momentum_save(p);
            // Update velocities
            RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1> *impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1>*>(rigid_body.impulse_accumulator);
            rigid_body.V()+=impulse_accumulator->accumulated_impulse.linear/rigid_body.Mass();
            rigid_body.Angular_Momentum()+=impulse_accumulator->accumulated_impulse.angular;
            rigid_body.Update_Angular_Velocity();
            // Euler step
            if(euler_step_with_new_velocity){
                rigid_body_collisions.collision_callbacks.Euler_Step_Position_With_New_Velocity(p,dt,time);
                rigid_body_collisions.skip_collision_check.Set_Last_Moved(p);}
            else{
                rigid_body_collisions.collision_callbacks.Restore_Position(p);
                rigid_body_collisions.collision_callbacks.Euler_Step_Position(p,dt,time);
                rigid_body_collisions.rigid_body_collection.Rigid_Body(p).Update_Bounding_Box();
                rigid_body_collisions.skip_collision_check.Set_Last_Moved(p);}}}
}
//#####################################################################
// Function Exchange_All_Pushes
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Exchange_All_Pushes(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TV>& mpi_rigid_X_save,ARRAY<ROTATION<TV> >& mpi_rigid_rotation_save,RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions)
{
    ARRAY<bool> need_to_reevolve(rigid_body_collection.rigid_body_particle.array_collection->Size());
    Prune_And_Exchange_Impulses(rigid_body_collection,need_to_reevolve);
        
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        if(need_to_reevolve(i)){
            int p=rigid_body_collection.rigid_body_particle.id(i);
            RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
            // Restore positions and rotations
            rigid_body_collection.rigid_body_particle.X(p)=mpi_rigid_X_save(p);
            rigid_body_collection.rigid_body_particle.rotation(p)=mpi_rigid_rotation_save(p);
            // Update positions and rotations
            RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1> *impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1>*>(rigid_body.impulse_accumulator);
            rigid_body.X()+=impulse_accumulator->accumulated_impulse.linear;
            rigid_body.Rotation()=ROTATION<TV>::From_Rotation_Vector(impulse_accumulator->accumulated_impulse.angular)*rigid_body.Rotation();rigid_body.Rotation().Normalize();
            rigid_body.Update_Angular_Velocity();
            rigid_body.Update_Bounding_Box();
            rigid_body_collisions.skip_collision_check.Set_Last_Moved(p);}}
}
//#####################################################################
// Function Prune_And_Exchange_Impulses
//#####################################################################
namespace PhysBAM{
template<> void THREADED_RIGIDS<VECTOR<float,1> >::Prune_And_Exchange_Impulses(RIGID_BODY_COLLECTION<VECTOR<float,1> >& rigid_body_collection,ARRAY<bool>& need_to_reevolve)
{PHYSBAM_NOT_IMPLEMENTED();}
template<> void THREADED_RIGIDS<VECTOR<double,1> >::Prune_And_Exchange_Impulses(RIGID_BODY_COLLECTION<VECTOR<double,1> >& rigid_body_collection,ARRAY<bool>& need_to_reevolve)
{PHYSBAM_NOT_IMPLEMENTED();}
}
template<class TV> void THREADED_RIGIDS<TV>::
Prune_And_Exchange_Impulses(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<bool>& need_to_reevolve)
{
#ifdef USE_PTHREADS
    // Figure out what impulse_accumulators this node needs to send
    ARRAYS_COMPUTATIONS::Fill(need_to_reevolve,false);
    ARRAY<int> accumulators_to_send;
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        int p=rigid_body_collection.rigid_body_particle.id(i);
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
        RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1> *impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1>*>(rigid_body.impulse_accumulator);
        if((impulse_accumulator->accumulated_impulse.linear!=TV()) || (impulse_accumulator->accumulated_impulse.angular!=T_SPIN())){
            need_to_reevolve(p)=true;
            accumulators_to_send.Append(p);}}
    
    // Send the impulses
    int buffer_size=sizeof(int)*(accumulators_to_send.m+1)+sizeof(TV)*accumulators_to_send.m+sizeof(T_SPIN)*accumulators_to_send.m+1;
    THREAD_PACKAGE pack(buffer_size);pack.send_tid=tid;pack.recv_tid=0;int position=0;
    *(int*)(&pack.buffer(position+1))=accumulators_to_send.m;position+=sizeof(int); // pack size
    for(int b=1;b<=accumulators_to_send.m;b++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(accumulators_to_send(b));
        RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1> *impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1>*>(rigid_body.impulse_accumulator);
        accumulators_to_send.Pack(pack.buffer,position,b);
        *(TV*)(&pack.buffer(position+1))=impulse_accumulator->accumulated_impulse.linear;position+=sizeof(TV); // pack linear impulse
        *(T_SPIN*)(&pack.buffer(position+1))=impulse_accumulator->accumulated_impulse.angular;position+=sizeof(T_SPIN);} // pack angular impulse

    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);

    pthread_barrier_wait(barr);

    // Process the received impulses
    for(int i=1;i<=buffers.m;i++){
        if(buffers(i).send_tid!=tid){
            int position=0;
            int accumulators_to_recv=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
            for(int b=1;b<=accumulators_to_recv;b++){
                int body_id=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
                TV linear_impulse=*(const TV*)(&buffers(i).buffer(position+1));position+=sizeof(TV);
                T_SPIN angular_impulse=*(const T_SPIN*)(&buffers(i).buffer(position+1));position+=sizeof(T_SPIN);
                RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(body_id);
                RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1> *impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1>*>(rigid_body.impulse_accumulator);
                impulse_accumulator->accumulated_impulse+=TWIST<TV>(linear_impulse,angular_impulse);
                need_to_reevolve(rigid_body.particle_index)=true;}}}

    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Exchange_Bounding_Box_Collision_Pairs
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Exchange_Bounding_Box_Collision_Pairs(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<VECTOR<int,2> >& pairs,const bool contact_pairs)
{
#ifdef USE_PTHREADS
    // Figure out what pairs to send to which nodes
    ARRAY<ARRAY<VECTOR<int,2> > > pairs_per_node(number_of_processors);
    if(!contact_pairs)
        for(int i=1;i<=pairs.m;i++){
            int j=pairs(i)(2);
            if(!Is_Real_Body(j) && !rigid_body_collection.Rigid_Body(j).Has_Infinite_Inertia()){
                int ghost_partition=Value(partition_id_from_particle_index(j));
                pairs_per_node(ghost_partition).Append(VECTOR<int,2>(j,pairs(i)(1)));}}
    else
        for(int i=1;i<=pairs.m;i++){
            int j=pairs(i)(1);
            if(!Is_Real_Body(j) && !rigid_body_collection.Rigid_Body(j).Has_Infinite_Inertia()){
                int ghost_partition=Value(partition_id_from_particle_index(j));
                pairs_per_node(ghost_partition).Append(VECTOR<int,2>(j,pairs(i)(2)));}}
    
    // Send the pairs
    int num_pairs_to_send=0;for(int i=1;i<=pairs_per_node.m;i++) num_pairs_to_send+=pairs_per_node(i).m;
    int buffer_size=sizeof(int)*(1+pairs_per_node.m+2*num_pairs_to_send)+1; //pairs_per_node.m,length of each sub array,2*total pairs sent,+1 for padding
    THREAD_PACKAGE pack(buffer_size);pack.send_tid=tid;pack.recv_tid=0;int position=0;
    *(int*)(&pack.buffer(position+1))=pairs_per_node.m;position+=sizeof(int); // pack size
    for(int ghost_partition=1;ghost_partition<=pairs_per_node.m;ghost_partition++){
        *(int*)(&pack.buffer(position+1))=pairs_per_node(ghost_partition).m;position+=sizeof(int); // pack sub array size
        for(int pair=1;pair<=pairs_per_node(ghost_partition).m;pair++){
            *(int*)(&pack.buffer(position+1))=pairs_per_node(ghost_partition)(pair)(1);position+=sizeof(int);
            *(int*)(&pack.buffer(position+1))=pairs_per_node(ghost_partition)(pair)(2);position+=sizeof(int);}}

    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);

    pthread_barrier_wait(barr);

    // Process the received pairs
    for(int i=1;i<=buffers.m;i++){
        if(buffers(i).send_tid!=tid){
            int position=0;
            int arrays_to_recv=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
            for(int array_index=1;array_index<=arrays_to_recv;array_index++){
                int pairs_to_recv=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
                if(array_index==tid){
                    for(int pair=1;pair<=pairs_to_recv;pair++){
                        int id_1=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
                        int id_2=*(const int*)(&buffers(i).buffer(position+1));position+=sizeof(int);
                        pairs.Append(VECTOR<int,2>(id_1,id_2));}}
                else {position+=pairs_to_recv*2*sizeof(int);}}}}

    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Synchronize_Dt
//#####################################################################
template<class TV> void THREADED_RIGIDS<TV>::
Synchronize_Dt(T& dt) const
{
#ifdef USE_PTHREADS
    THREAD_PACKAGE pack(sizeof(T));*(T*)(&pack.buffer(1))=dt;
    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers.m;i++) dt=min(*(T*)(&buffers(i).buffer(1)),dt);
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Reduce_Min
//#####################################################################
template<class TV> template<class T_DATA> T_DATA THREADED_RIGIDS<TV>::
Reduce_Min(const T_DATA local_value) const
{
#ifdef USE_PTHREADS
    T_DATA global_value=local_value;
    THREAD_PACKAGE pack(sizeof(T_DATA));
    *(T_DATA*)(&pack.buffer(1))=local_value;
    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers.m;i++) global_value=min(global_value,*(T_DATA*)(&buffers(i).buffer(1)));
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
    return global_value;
#else
    PHYSBAM_NOT_IMPLEMENTED();
#endif
}
//#####################################################################
// Function Reduce_Max
//#####################################################################
template<class TV> template<class T_DATA> T_DATA THREADED_RIGIDS<TV>::
Reduce_Max(const T_DATA local_value) const
{
#ifdef USE_PTHREADS
    T_DATA global_value=local_value;
    THREAD_PACKAGE pack(sizeof(T_DATA));
    *(T_DATA*)(&pack.buffer(1))=local_value;
    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers.m;i++) global_value=max(global_value,*(T_DATA*)(&buffers(i).buffer(1)));
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
    return global_value;
#else
    PHYSBAM_NOT_IMPLEMENTED();
#endif
}
//#####################################################################
// Function Reduce_Add
//#####################################################################
template<class TV> template<class T_DATA> T_DATA THREADED_RIGIDS<TV>::
Reduce_Add(const T_DATA local_value) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER(TV,T)       \
    template class THREADED_RIGIDS<TV>; \
    template T THREADED_RIGIDS<TV>::Reduce_Min(const T) const; \
    template T THREADED_RIGIDS<TV>::Reduce_Max(const T) const; \
    template int THREADED_RIGIDS<TV>::Reduce_Max(const int) const; \
    template T THREADED_RIGIDS<TV>::Reduce_Add(const T) const;
INSTANTIATION_HELPER(P(VECTOR<float,1>),float);
INSTANTIATION_HELPER(P(VECTOR<float,2>),float);
INSTANTIATION_HELPER(P(VECTOR<float,3>),float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(P(VECTOR<double,1>),double);
INSTANTIATION_HELPER(P(VECTOR<double,2>),double);
INSTANTIATION_HELPER(P(VECTOR<double,3>),double);
#endif
