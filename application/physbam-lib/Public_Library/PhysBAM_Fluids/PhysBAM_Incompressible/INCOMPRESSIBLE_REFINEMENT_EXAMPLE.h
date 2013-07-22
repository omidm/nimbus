//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_REFINEMENT_EXAMPLE__
#define __INCOMPRESSIBLE_REFINEMENT_EXAMPLE__
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_2D.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_3D.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_POLICY.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class INCOMPRESSIBLE_REFINEMENT_EXAMPLE:public RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename TURBULENCE_POLICY<TV>::TURBULENCE T_TURBULENCE;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    bool write_debug_data;
    std::string frame_title;
    int write_substeps_level;
    std::string output_directory;
    std::string split_dir;
    bool use_coarse_forces,use_interpolated_vorticity;
    int sub_scale;
    int restart;
    T kolmogorov;

    int number_of_ghost_cells;
    T cfl;

    GRID<TV> coarse_mac_grid,fine_mac_grid;
    MPI_UNIFORM_GRID<GRID<TV> > *coarse_mpi_grid,*fine_mpi_grid;
    ARRAY<T,FACE_INDEX<TV::dimension> > coarse_face_velocities;
    ARRAY<T,FACE_INDEX<TV::dimension> > fine_face_velocities;
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > projection;
    INCOMPRESSIBLE_UNIFORM<GRID<TV> > incompressible;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T>* advection_scalar;
    THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T>* threaded_advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary,*boundary_coarse;
    ARRAY<T,TV_INT> density;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    RIGID_GEOMETRY_COLLECTION<TV> rigid_geometry_collection;
    T_TURBULENCE turbulence;
    THREAD_QUEUE* thread_queue;

    INCOMPRESSIBLE_REFINEMENT_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~INCOMPRESSIBLE_REFINEMENT_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize_Fields()=0;
    virtual void Initialize_Confinement()=0;
    virtual void Get_Scalar_Field_Sources(const T time)=0;
    virtual void Set_Boundary_Conditions(const T time)=0;
    virtual void Initialize_Bodies() {}
    virtual void Preprocess_Projection(const T dt,const T time) {}
    virtual void Postprocess_Projection(const T dt,const T time) {}
    virtual void Preprocess_Frame(const int frame) {}
    virtual void Postprocess_Frame(const int frame) {}
    virtual void Initialize_MPI()=0;
    virtual void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time){}

//#####################################################################
};
}
#endif
