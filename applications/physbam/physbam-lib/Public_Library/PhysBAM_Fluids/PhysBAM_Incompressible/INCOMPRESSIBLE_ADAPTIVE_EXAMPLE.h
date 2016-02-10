//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_ADAPTIVE_EXAMPLE__
#define __INCOMPRESSIBLE_ADAPTIVE_EXAMPLE__
#include <PhysBAM_Tools/Grids_Uniform/GRID_ADAPTIVE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_ADAPTIVE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class INCOMPRESSIBLE_ADAPTIVE_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    std::string frame_title;
    int write_substeps_level;
    std::string output_directory;

    int number_of_ghost_cells;
    T cfl;

    GRID_ADAPTIVE<TV> mac_grid;
    FACE_ARRAY_ADAPTIVE<T,TV::dimension> face_velocities;
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > projection;
    INCOMPRESSIBLE_UNIFORM<GRID<TV> > incompressible;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T> advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    ARRAY<ARRAY<T,TV_INT>,TV_INT> density,temperature;

    INCOMPRESSIBLE_ADAPTIVE_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~INCOMPRESSIBLE_ADAPTIVE_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    virtual void Write_Output_Files(const int frame);
    virtual void Initialize_Fields()=0;
    virtual void Get_Scalar_Field_Sources(const T time,const GRID<TV>& mac_grid,ARRAY<T,TV_INT>& density)=0;
    virtual void Set_Boundary_Conditions(const T time,const GRID<TV>& mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >* projection)=0;
    virtual void Preprocess_Projection(const T dt,const T time) {}
    virtual void Postprocess_Projection(const T dt,const T time) {}

//#####################################################################
};
}
#endif
