//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VORTEX_PARTICLE_EVOLUTION_3D__
#define __VORTEX_PARTICLE_EVOLUTION_3D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <PhysBAM_Tools/Particles_Interpolation/SCATTERED_INTERPOLATION.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
namespace PhysBAM{

template<class T>
class VORTEX_PARTICLE_EVOLUTION_3D:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename GRID<TV>::FACE_ITERATOR T_FACE_ITERATOR;
    typedef typename MPI_GRID_POLICY<GRID<TV> >::MPI_GRID T_MPI_GRID;
public:
    VORTICITY_PARTICLES<VECTOR<T,3> > vorticity_particles;
    GRID<TV> grid;
    T_MPI_GRID* mpi_grid;
    RANDOM_NUMBERS<T> random;
    SCATTERED_INTERPOLATION<GRID<TV> > scattered_interpolation;
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > grid_vorticity,grid_vorticity_particles;
    T grid_confinement_parameter,particle_confinement_parameter,force_scaling;
    bool renormalize_vorticity_after_stretching_tilting;
    bool remove_grid_vorticity_from_particle_vorticity;
    bool apply_individual_particle_forces;
private:
//    T radius,radius_squared,one_over_radius_squared,one_over_radius_cubed;
public:

    VORTEX_PARTICLE_EVOLUTION_3D()
        :mpi_grid(0),grid_confinement_parameter(0),particle_confinement_parameter(0),force_scaling((T).03),renormalize_vorticity_after_stretching_tilting(false),
        remove_grid_vorticity_from_particle_vorticity(false),apply_individual_particle_forces(true)
    {
        scattered_interpolation.Use_Tent_Weights();
//        Set_Radius();
    }

    void Initialize(const GRID<TV>& grid_input)
    {grid=grid_input.Get_MAC_Grid();grid_vorticity.Resize(grid.Domain_Indices(2),false,false);grid_vorticity_particles.Resize(grid.Domain_Indices(2),false,false);
    /*scattered_interpolation.Set_Radius_Of_Influence(radius);*/}

//#####################################################################
private:
    T Gaussian_Kernel(T distance_squared);
public:
//    void Set_Radius(const T radius_input=(T).01);
    void Compute_Body_Force(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& force,const T dt,const T time);
    void Compute_Body_Force(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& force,const T dt,const T time);
    void Euler_Step(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T dt,const T time);
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const;
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& input_directory,const int frame);
//#####################################################################
};
}
#endif
