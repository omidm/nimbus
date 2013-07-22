#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VORTEX_PARTICLE_EVOLUTION_OCTREE__
#define __VORTEX_PARTICLE_EVOLUTION_OCTREE__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
namespace PhysBAM{

template<class T>
class VORTEX_PARTICLE_EVOLUTION_OCTREE:public NONCOPYABLE
{
public:
    VORTICITY_PARTICLES<VECTOR<T,3> > vorticity_particles;
    OCTREE_GRID<T>& grid;
    RANDOM_NUMBERS<T> random;
    ARRAY<VECTOR<T,3> > grid_vorticity;
    T grid_confinement_parameter,particle_confinement_parameter,force_scaling;
    bool renormalize_vorticity_after_stretching_tilting;
    bool remove_grid_vorticity_from_particle_vorticity;
    bool apply_individual_particle_forces;
    T radius,radius_squared,one_over_radius_squared,one_over_radius_cubed;
public:

    VORTEX_PARTICLE_EVOLUTION_OCTREE(OCTREE_GRID<T>& grid_input)
        :grid(grid_input),grid_confinement_parameter(0),particle_confinement_parameter(0),force_scaling((T).03),renormalize_vorticity_after_stretching_tilting(false),remove_grid_vorticity_from_particle_vorticity(false),apply_individual_particle_forces(true)
    {
        Set_Radius();
    }

    void Initialize()
    {grid_vorticity.Resize(grid.number_of_nodes,false,false);}

//#####################################################################
    T Gaussian_Kernel(T distance_squared);
    void Set_Radius(const T radius_input=(T).01);
    void Compute_Body_Force(const ARRAY<VECTOR<T,3> >& V_ghost,ARRAY<VECTOR<T,3> >& force,const T dt,const T time);
    void Euler_Step(const ARRAY<VECTOR<T,3> >& V_ghost,const T dt,const T time);
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const;
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& input_directory,const int frame);
//#####################################################################
};
}
#endif
#endif
