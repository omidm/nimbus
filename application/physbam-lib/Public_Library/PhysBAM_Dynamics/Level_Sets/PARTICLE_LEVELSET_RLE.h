#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_RLE
//#####################################################################
#ifndef __PARTICLE_LEVELSET_RLE__
#define __PARTICLE_LEVELSET_RLE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
namespace PhysBAM{

template<class T_GRID>
class PARTICLE_LEVELSET_RLE:public PARTICLE_LEVELSET<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::BLOCK_ITERATOR BLOCK_ITERATOR;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename T_GRID::BOX_HORIZONTAL T_BOX_HORIZONTAL;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR T_LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
public:
    typedef PARTICLE_LEVELSET<T_GRID> BASE;
    using BASE::half_band_width;using BASE::use_removed_negative_particles;using BASE::use_removed_positive_particles;using BASE::store_unique_particle_id;using BASE::last_unique_particle_id;
    using BASE::number_particles_per_cell;using BASE::maximum_particle_radius;using BASE::minimum_particle_radius;using BASE::random;using BASE::outside_particle_distance_multiplier;
    using BASE::maximum_iterations_for_attraction;using BASE::bias_towards_negative_particles;using BASE::Set_Number_Particles_Per_Cell;using BASE::particle_pool;
    using BASE::min_collision_distance_factor;using BASE::max_minus_min_collision_distance_factor_over_max_short;using BASE::Set_Band_Width;using BASE::positive_particles;
    using BASE::negative_particles;using BASE::levelset;using BASE::removed_negative_particles;using BASE::removed_positive_particles;using BASE::escaped_positive_particles;
    using BASE::escaped_negative_particles;using BASE::cfl_number;using BASE::template_particles;using BASE::template_removed_particles;

    const T_GRID& grid;
    ARRAY<T> phi;
    bool use_removed_negative_particles_in_long_cells,use_removed_positive_particles_in_long_cells;
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> *removed_negative_particles_in_long_cells,*removed_positive_particles_in_long_cells;
    MPI_RLE_GRID<T_GRID>* mpi_grid;

    PARTICLE_LEVELSET_RLE(T_GRID& grid_input);
    virtual ~PARTICLE_LEVELSET_RLE();

    template<class T_PARTICLES> T_PARTICLES* Allocate_Particle()
    {return PARTICLE_LEVELSET<T_GRID>::template Allocate_Particle<T_PARTICLES>();}

//#####################################################################
    void Initialize_Particle_Levelset_Grid_Values();
    void Seed_Particles(const T time,const bool verbose=true);
    void Adjust_Particle_Radii();
    void Modify_Levelset_Using_Escaped_Particles();
    void Euler_Step_Particles(const ARRAY<T>& V,const T dt,const T time,const bool use_second_order_for_nonremoved_particles=false,const bool update_particle_cells_after_euler_step=true,
        const bool verbose=true);
    void Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    int Reseed_Particles(const T time,ARRAY<bool>* cell_centered_mask=0);
    void Delete_Particles_Outside_Grid();
    void Identify_And_Remove_Escaped_Particles(const ARRAY<T>& V,const T radius_fraction=1.5,const bool verbose=true);
    void Reincorporate_Removed_Particles(const T radius_fraction);
    void Compact_Particles_Into_Single_Particle_Bin();
    void Delete_Particles_Far_From_Interface(const int discrete_band=4);
    void Delete_Particles_In_Local_Maximum_Phi_Cells();
    void Transfer_Particles(const T_GRID& new_grid);
protected:
    void Update_Particle_Cells(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles);
    void Update_Particle_Cells(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles_in_long_cells,const T dt);
    void Compact_Particles_Into_Single_Particle_Bin(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign);
    void Euler_Step_Removed_Particles(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles_in_long_cells,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    void Euler_Step_Particles(const ARRAY<T>& V,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
        const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool assume_particles_in_correct_cells=true);
    void Second_Order_Runge_Kutta_Step_Particles(const ARRAY<T>& V,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
        const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    bool Attract_Individual_Particle_To_Interface_And_Adjust_Radius(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const T phi_min,const T phi_max,
        const BLOCK_ITERATOR& block,const int index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    void Adjust_Particle_Radii(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign);
    void Modify_Levelset_Using_Escaped_Particles(ARRAY<T>& phi,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const int sign);
    int Reseed_Delete_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const int sign);
    int Reseed_Add_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& other_particles,
        const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,ARRAY<bool>* cell_centered_mask);
    void Identify_Escaped_Particles(const int sign=0);
    void Identify_Escaped_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign);
    template<class T_FACE_LOOKUP_LOOKUP> int Remove_Escaped_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles,const T_FACE_LOOKUP_LOOKUP& V,const T radius_fraction);
    int Remove_Escaped_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,const T radius_fraction);
    void Reincorporate_Removed_Particles(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>*& particles,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& removed_particles,const T radius_fraction);
    void Delete_Particles_Outside_Grid(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_in_long_cells);
    template<class T_PARTICLE_CLASS> void Delete_Particles_Outside_Grid(ARRAY<T_PARTICLE_CLASS*>& particles,const RLE_GRID_2D<T>&);
    template<class T_PARTICLE_CLASS> void Delete_Particles_Outside_Grid(ARRAY<T_PARTICLE_CLASS*>& particles,const RLE_GRID_3D<T>&);
    void Transfer_Particles_In_Long_Cells(const T_GRID& new_grid,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& new_particles,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_in_long_cells);
//#####################################################################
};
}
#endif
#endif
