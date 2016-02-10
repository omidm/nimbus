//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_UNIFORM
//#####################################################################
#ifndef __PARTICLE_LEVELSET_UNIFORM__
#define __PARTICLE_LEVELSET_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <PhysBAM_Tools/Parallel_Computation/PTHREAD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET.h>
namespace PhysBAM{

template<class T_GRID> class FLUID_MASS_CONSERVATION_MODIFY_PHI;

template<class T_GRID>
class PARTICLE_LEVELSET_UNIFORM:public PARTICLE_LEVELSET<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<char>::TYPE T_ARRAYS_CHAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<bool> >::TYPE T_ARRAYS_ARRAY_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<TV>*>::TYPE T_ARRAYS_ARRAY_TV;
    typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_LINEAR_INTERPOLATION_SCALAR::template REBIND<TV>::TYPE T_LINEAR_INTERPOLATION_VECTOR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    typedef PARTICLE_LEVELSET<T_GRID> BASE;
    using BASE::half_band_width;using BASE::use_removed_negative_particles;using BASE::use_removed_positive_particles;using BASE::store_unique_particle_id;using BASE::last_unique_particle_id;
    using BASE::number_particles_per_cell;using BASE::maximum_particle_radius;using BASE::minimum_particle_radius;using BASE::random;using BASE::outside_particle_distance_multiplier;
    using BASE::maximum_iterations_for_attraction;using BASE::bias_towards_negative_particles;using BASE::Set_Number_Particles_Per_Cell;
    using BASE::min_collision_distance_factor;using BASE::max_minus_min_collision_distance_factor_over_max_short;
    using BASE::reincorporate_removed_particles_everywhere;using BASE::save_removed_particle_times;using BASE::removed_particle_times;using BASE::only_use_negative_particles;
    using BASE::velocity_interpolation_collidable_contour_value;using BASE::template_particles;using BASE::template_removed_particles;
    using BASE::levelset;using BASE::negative_particles;using BASE::positive_particles;using BASE::removed_negative_particles;using BASE::removed_positive_particles;
    using BASE::escaped_negative_particles;using BASE::escaped_positive_particles;using BASE::Set_Band_Width;using BASE::cfl_number;using BASE::number_of_ghost_cells;using BASE::particle_pool;

    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;

    VOF_ADVECTION<TV>* vof_advection;
    T_ARRAYS_ARRAY_TV positive_particle_positions,negative_particle_positions;

    THREAD_QUEUE* thread_queue;
#ifdef USE_PTHREADS
    pthread_mutex_t cell_lock;
    pthread_barrier_t cell_barr;
#endif

    PARTICLE_LEVELSET_UNIFORM(T_GRID& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells_input);
    ~PARTICLE_LEVELSET_UNIFORM();

    void Set_Thread_Queue(THREAD_QUEUE* thread_queue_input)
    {thread_queue=thread_queue_input;}

//#####################################################################
    void Seed_Particles(const T time,const bool verbose=true);
    void Adjust_Particle_Radii();
    void Adjust_Particle_Radii_Threaded(RANGE<TV_INT>& domain);
    void Modify_Levelset_Using_Escaped_Particles(T_FACE_ARRAYS_SCALAR* V,ARRAY<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES*>* other_positive_particles=0);
    void Update_Particles_To_Reflect_Mass_Conservation(T_ARRAYS_SCALAR& phi_old,const bool update_particle_cells=true,const bool verbose=false);
    void Update_Particles_To_Reflect_Mass_Conservation(const T_FAST_LEVELSET& levelset_old,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const bool update_particle_cells,const bool verbose);
    void Euler_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,const T dt,const T time,const bool use_second_order_for_nonremoved_particles=false,
        const bool update_particle_cells_after_euler_step=true,const bool verbose=true,const bool analytic_test=false,const bool default_to_forward_euler_outside=false);
    void Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    void Euler_Step_Removed_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    template<class T_ARRAYS_PARTICLES> void Euler_Step_Particles_Wrapper(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool update_particle_cells_after_euler_step=true,const bool assume_particles_in_correct_blocks=true,const bool enforce_domain_boundaries=true);
    template<class T_ARRAYS_PARTICLES> void Euler_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool update_particle_cells_after_euler_step=true,const bool assume_particles_in_correct_blocks=true,const bool enforce_domain_boundaries=true);
    template<class T_ARRAYS_PARTICLES> void Euler_Step_Particles_Threaded(RANGE<TV_INT>& domain,const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool assume_particles_in_correct_blocks=true,const bool enforce_domain_boundaries=true);
    template<class T_ARRAYS_PARTICLES> void Second_Order_Runge_Kutta_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true,const bool default_to_forward_euler_outside=false);
    template<class T_ARRAYS_PARTICLES> void Second_Order_Runge_Kutta_Step_Particles_Threaded(RANGE<TV_INT>& domain,const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool verbose=true,const bool default_to_forward_euler_outside=false);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells(T_ARRAYS_PARTICLES& particles);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_One_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_Two_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_Three_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Find_List(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process);
    template<class T_ARRAYS_PARTICLES> void Delete_Marked_Particles(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles);
    template<class T_ARRAYS_PARTICLES> void Consistency_Check(RANGE<TV_INT> domain,T_ARRAYS_PARTICLES& particles);
    void Advect_Particles_In_Ghost_Region(T_GRID& grid,T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T dt,int number_of_ghost_cells);
    int Reseed_Particles_In_Ghost_Region(const T time,int number_of_ghost_cells);
    int Reseed_Delete_Particles_In_Whole_Region(int number_of_ghost_cells,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign,T_ARRAYS_BOOL* cell_centered_mask=0);
    int Reseed_Add_Particles_In_Whole_Region(int number_of_ghost_cells,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,T_ARRAYS_BOOL* cell_centered_mask=0);
    int Reseed_Particles(const T time,T_ARRAYS_BOOL* cell_centered_mask=0);
    int Reseed_Delete_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign);
    void Identify_Escaped_Particles(RANGE<TV_INT>& domain,const int sign=0);
    void Delete_All_Particles_In_Cell(const TV_INT& block);
    void Delete_Positive_Particles_In_Cell(const TV_INT& block);
    void Delete_Deep_Escaped_Particles(const T radius_fraction=1.5,const bool need_to_identify_escaped_particles=true,const bool verbose=false);
    void Delete_Particles_Outside_Grid();
    void Delete_Particles_Far_Outside_Grid(RANGE<TV_INT>& domain);
    void Delete_Particles_Near_Outside_Grid(RANGE<TV_INT>& domain,int axis,int side);
    void Delete_Particles_In_Local_Maximum_Phi_Cells(const int sign);
    void Delete_Particles_In_Local_Maximum_Phi_Cells_Threaded(RANGE<TV_INT>& domain,const int sign,const T tolerance);
    void Identify_And_Remove_Escaped_Particles(const T_FACE_ARRAYS_SCALAR& V,const T radius_fraction=1.5,const T time=0,const bool verbose=true);
    void Identify_And_Remove_Escaped_Particles_Threaded(RANGE<TV_INT>& domain,const T_FACE_ARRAYS_SCALAR& V,const T radius_fraction=1.5,const T time=0,const bool verbose=true);
    void Identify_And_Remove_Escaped_Positive_Particles(const T_FACE_ARRAYS_SCALAR& V,const T radius_fraction=1.5,const T time=0,const bool verbose=true);
    void Reincorporate_Removed_Particles(const T radius_fraction,const T mass_scaling,T_FACE_ARRAYS_SCALAR* V,const bool conserve_momentum_for_removed_negative_particles=true);
    void Reincorporate_Removed_Particles_Threaded(RANGE<TV_INT>& domain,const T radius_fraction,const T mass_scaling,T_FACE_ARRAYS_SCALAR* V);
    bool Fix_Momentum_With_Escaped_Particles(const T_FACE_ARRAYS_SCALAR& V,const T_ARRAYS_SCALAR& momentum_lost,const T radius_fraction,const T mass_scaling,const T time,const bool force=true);
    bool Fix_Momentum_With_Escaped_Particles(const TV& location,const T_FACE_ARRAYS_SCALAR& V,const T momentum_lost,const T radius_fraction,const T mass_scaling,const T time,const bool force=true);

    void Delete_Particles_Far_From_Interface(const int discrete_band=4);
    void Delete_Particles_Far_From_Interface_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface);
    void Delete_Particles_Far_From_Interface_Part_Two(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface,const int discrete_band=4);
    void Delete_Particles_Far_From_Interface_Part_Three(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface);
    void Add_Negative_Particle(const TV& particle_location,const TV& particle_velocity,const unsigned short quantized_collision_distance);
    void Exchange_Overlap_Particles();
    template<class T_ARRAYS_PARTICLES> void Modify_Levelset_Using_Escaped_Particles_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi,T_ARRAYS_PARTICLES& particles,T_FACE_ARRAYS_SCALAR* V,const int sign,const T one_over_radius_multiplier);
protected:
    bool Attract_Individual_Particle_To_Interface_And_Adjust_Radius(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles,
        const T phi_min,const T phi_max,const BLOCK_UNIFORM<T_GRID>& block,const int index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const bool delete_particles_that_leave_original_cell,
        RANDOM_NUMBERS<T>& local_random,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index);
    void Adjust_Particle_Radii(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign);
    template<class T_ARRAYS_PARTICLES> void Modify_Levelset_Using_Escaped_Particles(T_ARRAYS_SCALAR& phi,T_ARRAYS_PARTICLES& particles,T_FACE_ARRAYS_SCALAR* V,const int sign);
    int Reseed_Add_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,T_ARRAYS_BOOL* cell_centered_mask);
    void Reseed_Add_Particles_Threaded_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,T_ARRAYS_BOOL* cell_centered_mask,ARRAY<int,TV_INT>& number_of_particles_to_add);
    void Reseed_Add_Particles_Threaded_Part_Two(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const ARRAY<int,TV_INT>& number_of_particles_to_add,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,const ARRAY<int,TV_INT>* domain_index);
    void Copy_From_Move_List_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,ARRAY<ARRAY<TRIPLE<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT,int> >,TV_INT>& move_particles,const ARRAY<int,TV_INT>* domain_index);
    void Identify_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign);
    int Delete_Deep_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign,const T radius_fraction,
        const bool need_to_identify_escaped_particles);
    void Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int side);
    void Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int side);
    template<class T_FACE_LOOKUP_LOOKUP> int Remove_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles,const T_FACE_LOOKUP_LOOKUP& V,const T radius_fraction,const T time);
    int Remove_Escaped_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,const T radius_fraction);
    void Reincorporate_Removed_Particles(const BLOCK_UNIFORM<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>*& particles,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& removed_particles,const T radius_fraction,const T mass_scaling,T_FACE_ARRAYS_SCALAR* V);
//#####################################################################
};
}
#endif
