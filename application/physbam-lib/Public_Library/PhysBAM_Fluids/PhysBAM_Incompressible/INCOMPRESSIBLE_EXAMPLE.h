//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_EXAMPLE__
#define __INCOMPRESSIBLE_EXAMPLE__
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
//#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Fluid_Evolution/INCOMPRESSIBLE_FLUID_EVOLUTION.h>
namespace PhysBAM{

template<class TV>
class INCOMPRESSIBLE_EXAMPLE:public RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    int restart;
    std::string frame_title;
    int write_substeps_level;
    bool write_debug_data;
    bool write_output_files;
    bool analytic_test;
    int order;
    std::string output_directory;

    int number_of_ghost_cells;
    T cfl;
    bool use_viscosity;

    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<GRID<TV> > *mpi_grid;
    THREAD_QUEUE* thread_queue;
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >* projection;
    INCOMPRESSIBLE_UNIFORM<GRID<TV> > incompressible;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities,face_velocities_save;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>,T, AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> > >,LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > > advection_scalar;
    //ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T> advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary;
    ARRAY<T,TV_INT> density;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;    
    RIGID_GEOMETRY_COLLECTION<TV> rigid_geometry_collection;

    INCOMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input,const int number_of_threads=1,const int scaling_factor=1);
    virtual ~INCOMPRESSIBLE_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize_Fields()=0;
    virtual void Initialize_Confinement() {}
    virtual void Get_Scalar_Field_Sources(const T time)=0;
    virtual void Set_Boundary_Conditions(const T time)=0;
    virtual void Initialize_Bodies() {}
    virtual void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time){}

//#####################################################################
};
}
#endif
