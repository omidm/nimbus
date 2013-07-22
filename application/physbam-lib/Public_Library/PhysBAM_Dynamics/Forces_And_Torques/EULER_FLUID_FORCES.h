//#####################################################################
// Copyright 2007, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_FLUID_FORCES  
//#####################################################################
#ifndef __EULER_FLUID_FORCES__
#define __EULER_FLUID_FORCES__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/SOLIDS_FORCES.h>
#include <cfloat>
namespace PhysBAM{

template<class T_GRID> class FLUID_COLLISION_BODY_INACCURATE_UNION;
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class T_GRID>
class EULER_FLUID_FORCES:public SOLIDS_FORCES<typename T_GRID::VECTOR_T>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef SOLIDS_FORCES<typename T_GRID::VECTOR_T> BASE;
    typedef typename BASE::RIGID_FREQUENCY_DATA RIGID_FREQUENCY_DATA;
    typedef typename BASE::DEFORMABLE_FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;

    using BASE::particles;

    T_GRID grid;
    const T_FACE_ARRAYS_SCALAR& pressure_at_faces;
    const T_FACE_ARRAYS_BOOL& solid_fluid_face;
    const T_ARRAYS_BOOL& cells_inside_fluid;
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_bodies_affecting_fluid;

public:
    EULER_FLUID_FORCES(const T_GRID& grid_input,const T_FACE_ARRAYS_SCALAR& pressure_at_faces_input,
        const T_FACE_ARRAYS_BOOL& solid_fluid_face_input,const T_ARRAYS_BOOL& cells_inside_fluid_input,
        const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_bodies_affecting_fluid_input,PARTICLES<TV>& particles_input,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);

    virtual ~EULER_FLUID_FORCES();

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE
    {}

    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE
    {return FLT_MAX;} // TODO

    void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<RIGID_FREQUENCY_DATA> rigid_frequency) PHYSBAM_OVERRIDE
    {}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
