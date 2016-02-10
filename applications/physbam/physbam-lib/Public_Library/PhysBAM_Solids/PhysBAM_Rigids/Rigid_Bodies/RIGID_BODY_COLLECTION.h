//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_COLLECTION
//#####################################################################
#ifndef __RIGID_BODY_COLLECTION__
#define __RIGID_BODY_COLLECTION__

#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>

namespace PhysBAM{

template<class TV> class RIGID_GEOMETRY;
template<class TV> class RIGID_BODY_CLUSTER_BINDINGS;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class RIGIDS_FORCES;
template<class TV> class RIGID_BODY_EVOLUTION_PARAMETERS;

template<class TV>
class RIGID_BODY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    RIGID_BODY_PARTICLES<TV> rigid_body_particle;
    RIGID_GEOMETRY_COLLECTION<TV> rigid_geometry_collection;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings;
    RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>* rigids_example_forces_and_velocities;
    ARRAY<int> simulated_rigid_body_particles;
    ARRAY<int> dynamic_rigid_body_particles;
    ARRAY<int> &static_rigid_bodies,&kinematic_rigid_bodies,static_and_kinematic_rigid_bodies;
    ARRAY<RIGIDS_FORCES<TV>*> rigids_forces;
    
    bool print_diagnostics;
    bool print_residuals;
    bool print_energy;
    int iterations_used_diagnostic;

    RIGID_BODY_COLLECTION(RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>* rigids_example_forces_and_velocities_input,COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input,
        ARRAY_COLLECTION* array_collection=0);
    virtual ~RIGID_BODY_COLLECTION();

    RIGID_BODY_STATE<TV> State(const int particle) const
    {return RIGID_BODY_STATE<TV>(rigid_geometry_collection.Rigid_Geometry(particle).Frame(),rigid_geometry_collection.Rigid_Geometry(particle).Twist());}

    void Set_State(const int particle,const RIGID_BODY_STATE<TV>& state)
    {rigid_body_particle.X(particle)=state.frame.t;rigid_body_particle.rotation(particle)=state.frame.r;rigid_body_particle.V(particle)=state.twist.linear;rigid_body_particle.angular_velocity(particle)=state.twist.angular;}

    bool Exists(const int particle) const
    {return particle>0 && particle<=rigid_body_particle.array_collection->Size() && rigid_body_particle.rigid_geometry(particle);}

    bool Is_Active(const int particle) const
    {return Exists(particle) && Rigid_Body(particle).particle_index>0;}

    void Reset_Impulse_Accumulators();
    RIGID_BODY<TV>& Rigid_Body(const int particle_index);
    const RIGID_BODY<TV>& Rigid_Body(const int particle_index) const;
    int Add_Rigid_Body(RIGID_BODY<TV>* rigid_body,const int simplicial_boundary_id,const int implicit_object_id,const int simplicial_interior_id);
    int Add_Rigid_Body_And_Geometry(RIGID_BODY<TV>* rigid_body);
    int Add_Rigid_Body(const STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor=1,const bool read_simplicial_boundary=true,const bool read_implicit_object=true,
        const bool read_simplicial_interior=false,const bool read_rgd_file=true);
    int Add_Rigid_Body(const STREAM_TYPE stream_type,const bool thin_shell,const std::string& basename,const T scaling_factor=1,const bool read_simplicial_boundary=true,
        const bool read_implicit_object=true,const bool read_simplicial_interior=false,const bool read_rgd_file=true);
    void Update_Angular_Velocity();
    void Update_Angular_Momentum();
    void Update_Angular_Velocity(const ARRAY<int>& rigid_body_particles);
    void Update_Angular_Momentum(const ARRAY<int>& rigid_body_particles);
    void Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame,ARRAY<int>* needs_init=0,ARRAY<int>* needs_destroy=0);
    void Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const; // TODO: optionally skip certain kinds of structures in output
    void Update_Simulated_Particles();

    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
    void Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T scale,const T time) const;

    void Update_Position_Based_State(const T time);
    void Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const;
    void Print_Energy(const T time,const int step) const;
    T CFL_Rigid(const RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters,const bool verbose_dt);
    int Add_Force(RIGIDS_FORCES<TV>* force);
    void Save_Potential_Energy(const T time);
    void Compute_Energy_Error(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt);
    void Add_Energy_Correction_Force(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt);
    void Compute_Previously_Applied_Forces();
    void Store_Velocities();
};
}
#endif
