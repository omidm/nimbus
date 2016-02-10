#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VORTEX_PARTICLE_EVOLUTION_RLE_3D__
#define __VORTEX_PARTICLE_EVOLUTION_RLE_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
namespace PhysBAM{

template<class T>
class VORTEX_PARTICLE_EVOLUTION_RLE_3D:public NONCOPYABLE
{
    typedef RLE_GRID_3D<T> T_GRID;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename TV::SPIN T_SPIN;
    typedef typename T_GRID::BLOCK T_BLOCK;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
public:
    VORTICITY_PARTICLES<TV> vorticity_particles;
    T_GRID& grid;
    ARRAY<TV> grid_vorticity;
    T grid_confinement_parameter,particle_confinement_parameter,force_scaling;
    bool renormalize_vorticity_after_stretching_tilting;
    bool remove_grid_vorticity_from_particle_vorticity;
    T radius;
private:
    T radius_squared,one_over_radius_squared,one_over_radius_cubed;
    const NONLINEAR_FUNCTION<T(T)>* force_scaling_from_age;
public:

    VORTEX_PARTICLE_EVOLUTION_RLE_3D(T_GRID& grid_input)
        :grid(grid_input),grid_confinement_parameter(0),particle_confinement_parameter(0),force_scaling((T).03),renormalize_vorticity_after_stretching_tilting(false),
        remove_grid_vorticity_from_particle_vorticity(false),force_scaling_from_age(0)
    {
        Set_Radius();
    }

    void Initialize_Grid()
    {grid_vorticity.Resize(grid.number_of_cells,false,false);}

    void Set_Force_Scaling_From_Age(const NONLINEAR_FUNCTION<T(T)>& force_scaling_from_age_input)
    {force_scaling_from_age=&force_scaling_from_age_input;}

//#####################################################################
private:
    T Gaussian_Kernel(T distance_squared);
public:
    void Set_Radius(const T radius_input=(T).01);
    void Compute_Body_Force(const ARRAY<T>& V_ghost,ARRAY<T>& force,const T dt,const T time);
    void Euler_Step(const ARRAY<T>& V_ghost,const T dt,const T time);
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& input_directory,const int frame);
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const;
//#####################################################################
};
}
#endif
#endif
