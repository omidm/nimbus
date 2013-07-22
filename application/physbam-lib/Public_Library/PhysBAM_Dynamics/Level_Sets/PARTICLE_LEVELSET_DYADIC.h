#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_DYADIC
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_DYADIC__
#define __PARTICLE_LEVELSET_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
namespace PhysBAM{

template<class T> class COLLISION_BODY_LIST_3D;

template<class T_GRID>
class PARTICLE_LEVELSET_DYADIC:public PARTICLE_LEVELSET<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
    typedef typename T_GRID::UNIFORM_GRID UNIFORM_GRID;
    typedef typename UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;
    typedef typename UNIFORM_GRID::NODE_ITERATOR UNIFORM_NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<UNIFORM_GRID>::ARRAYS_SCALAR UNIFORM_ARRAYS;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    typedef PARTICLE_LEVELSET<T_GRID> BASE;
    using BASE::half_band_width;using BASE::use_removed_negative_particles;using BASE::use_removed_positive_particles;using BASE::store_unique_particle_id;using BASE::last_unique_particle_id;
    using BASE::number_particles_per_cell;using BASE::maximum_particle_radius;using BASE::minimum_particle_radius;using BASE::random;using BASE::outside_particle_distance_multiplier;
    using BASE::maximum_iterations_for_attraction;using BASE::bias_towards_negative_particles;using BASE::Set_Number_Particles_Per_Cell;using BASE::particle_pool;
    using BASE::min_collision_distance_factor;using BASE::max_minus_min_collision_distance_factor_over_max_short;
    using BASE::velocity_interpolation_collidable_contour_value;using BASE::reincorporate_removed_particles_everywhere;using BASE::Set_Band_Width;using BASE::Adjust_Particle_For_Objects;
    using BASE::levelset;using BASE::negative_particles;using BASE::positive_particles;using BASE::removed_negative_particles;
    using BASE::removed_positive_particles;using BASE::escaped_negative_particles;using BASE::escaped_positive_particles;

    PARTICLE_LEVELSET_DYADIC(T_GRID& grid_input,ARRAY<T>& phi_input);
    ~PARTICLE_LEVELSET_DYADIC();

private:
    template<class T_PARTICLES> T_PARTICLES* Allocate_Particle()
    {return PARTICLE_LEVELSET<T_GRID>::template Allocate_Particle<T_PARTICLES>();}
public:

    void Delete_All_Particles_In_Cell(const int i)
    {Free_Particle_And_Clear_Pointer(positive_particles(i));
    Free_Particle_And_Clear_Pointer(negative_particles(i));
    if(use_removed_positive_particles) Free_Particle_And_Clear_Pointer(removed_positive_particles(i));
    if(use_removed_negative_particles) Free_Particle_And_Clear_Pointer(removed_negative_particles(i));}

    void Update_Particle_Cells(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles)
    {ARRAY<int> number_of_particles_per_cell(particles.m);
    for(int i=1;i<=particles.m;i++)if(particles(i)) number_of_particles_per_cell(i)=particles(i)->Number();
    ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*> particles_links;particles_links.Preallocate(10);
    for(int i=1;i<=particles.m;i++){
        CELL* current_cell=cell_pointer_from_index(i);
        particles_links.Remove_All();PARTICLE_LEVELSET_PARTICLES<TV>* particles_link=particles(i);
        while(particles_link){particles_links.Append(particles_link);particles_link=particles_link->next;}
        for(int link=particles_links.m;link>=1;link--) for(int k=particles_links(link)->array_collection->Size();k>=1;k--){
            CELL* final_base_cell=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,particles_links(link)->X(k));
            if(!final_base_cell){Delete_Particle(*particles_links(link),k);continue;}
            if(final_base_cell!=current_cell){
                if(!particles(final_base_cell->Cell())) particles(final_base_cell->Cell())=Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >();
                int absolute_index=(link-1)*particle_pool.number_particles_per_cell+k;
                Move_Particle(*particles_links(link),*particles(final_base_cell->Cell()),absolute_index);}}}
    for(int i=1;i<=particles.m;i++)if(particles(i) && !particles(i)->array_collection->Size()) Free_Particle_And_Clear_Pointer(particles(i));}

    void Update_Particle_Cells(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles)
    {ARRAY<int> number_of_particles_per_cell(particles.m);
    for(int i=1;i<=particles.m;i++)if(particles(i)) number_of_particles_per_cell(i)=particles(i)->array_collection->Size();
    ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    for(int i=1;i<=particles.m;i++)for(int k=number_of_particles_per_cell(i);k>=1;k--){
        CELL* current_cell=cell_pointer_from_index(i);
        CELL* final_base_cell=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,particles(i)->X(k));
        if(!final_base_cell){Delete_Particle(*particles(i),k);continue;}
        if(final_base_cell!=current_cell){
            if(!particles(final_base_cell->Cell())) particles(final_base_cell->Cell())=Allocate_Particle<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >();
            Move_Particle(*particles(i),*particles(final_base_cell->Cell()),k);}}
    for(int i=1;i<=particles.m;i++)if(particles(i) && !particles(i)->array_collection->Size()){delete particles(i);particles(i)=0;}}
      
    void Compact_Particles_Into_Single_Particle_Bin()
    {ARRAY<CELL*>& cell_pointer_from_index=levelset.grid.Cell_Pointer_From_Index();
    for(int i=1;i<=negative_particles.m;i++) if(negative_particles(i))
        Compact_Particles_Into_Single_Particle_Bin(*negative_particles(i),BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(i)),-1);
    for(int i=1;i<=positive_particles.m;i++) if(positive_particles(i))
        Compact_Particles_Into_Single_Particle_Bin(*positive_particles(i),BLOCK_DYADIC<T_GRID>(levelset.grid,cell_pointer_from_index(i)),1);}
    
    void Compact_Particles_Into_Single_Particle_Bin(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const BLOCK_DYADIC<T_GRID>& block,const int sign)
    {if(!particles.next) return;
    PARTICLE_LEVELSET_PARTICLES<TV>* particles_link=&particles;
    int heap_size=0;ARRAY<int> heap_particle_indices(particle_pool.number_particles_per_cell);ARRAY<T> heap_phi_minus_radius(particle_pool.number_particles_per_cell);
    for(int index=particles_link->array_collection->Size();index>=1;index--){
        T phi_minus_radius=sign*levelset.Phi(block,particles_link->X(index))-particles_link->radius(index);
        heap_size++;heap_particle_indices(heap_size)=index;heap_phi_minus_radius(heap_size)=phi_minus_radius;}
    ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices); // when heap is full, order values with largest on top 
    while(particles_link->next){
        particles_link=particles_link->next;assert(!particles_link->next||particles_link->array_collection->Size()==particle_pool.number_particles_per_cell);
        for(int index=particles_link->array_collection->Size();index>=1;index--){
            T phi_minus_radius=sign*levelset.Phi(block,particles_link->X(index))-particles_link->radius(index);
            if(phi_minus_radius < heap_phi_minus_radius(1)){ // delete particle on top of heap & add new particle
                particles.array_collection->Copy_Element(*particles_link->array_collection,index,heap_particle_indices(1));
                heap_phi_minus_radius(1)=phi_minus_radius;
                ARRAYS_COMPUTATIONS::Heapify(heap_phi_minus_radius,heap_particle_indices,1,heap_phi_minus_radius.m);}}}
    Free_Particle_And_Clear_Pointer(particles.next);}

//#####################################################################       
    void Seed_Particles(const T time,const bool verbose=true);
    void Adjust_Particle_Radii();
    void Modify_Levelset_Using_Escaped_Particles();
    void Euler_Step_Particles(const ARRAY<T>& face_velocities,const T dt,const T time,const bool use_second_order_for_nonremoved_particles=false,
        const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    void Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    void Euler_Step_Removed_Particles(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    void Euler_Step_Particles(const ARRAY<T>& face_velocities,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
        const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool assume_particles_in_correct_cells=true);
    void Second_Order_Runge_Kutta_Step_Particles(const ARRAY<T>& face_velocities,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    int Reseed_Particles(const T time,ARRAY<bool>* cell_centered_mask=0);
    int Reseed_Delete_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const int sign); 
    void Identify_Escaped_Particles(const int sign=0);
    void Identify_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign);    
    void Delete_Deep_Escaped_Particles(const T radius_fraction=1.5,const bool need_to_identify_escaped_particles=true,const bool verbose=false);
    void Delete_Particles_Outside_Grid();
    void Identify_And_Remove_Escaped_Particles(const ARRAY<T>& face_velocities,const T radius_fraction=1.5,const bool verbose=true);
    void Reincorporate_Removed_Particles(const T radius_fraction);
    void Delete_Particles_Far_From_Interface(const int discrete_band=4);
    void Add_Negative_Particle(const TV& particle_location,const T collision_distance_percentage);
private: 
    bool Attract_Individual_Particle_To_Interface_And_Adjust_Radius(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const T phi_min,const T phi_max,const BLOCK_DYADIC<T_GRID>& block,
        const int index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const bool delete_particles_that_leave_original_cell=false);
    void Adjust_Particle_Radii(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign);
    template<class T_PARTICLES> void Modify_Levelset_Using_Escaped_Particles(ARRAY<T>& phi,ARRAY<T_PARTICLES*>& particles,const int sign);
    int Reseed_Add_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& other_particles,const int sign,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,ARRAY<bool>* cell_centered_mask); 
    int Delete_Deep_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign,const T radius_fraction,
        const bool need_to_identify_escaped_particles);    
    void Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int side);
    template<class T_FACE_LOOKUP_LOOKUP> int Remove_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles,const T_FACE_LOOKUP_LOOKUP& face_velocities,const T radius_fraction);
    int Remove_Escaped_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,const T radius_fraction);
    void Reincorporate_Removed_Particles(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>*& particles,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& removed_particles,const T radius_fraction);
//#####################################################################
};   
}
#endif
#endif
