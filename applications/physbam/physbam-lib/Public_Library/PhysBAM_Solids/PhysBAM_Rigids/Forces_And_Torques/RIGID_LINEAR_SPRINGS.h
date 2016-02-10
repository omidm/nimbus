//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_LINEAR_SPRINGS
//#####################################################################
#ifndef __RIGID_LINEAR_SPRINGS__
#define __RIGID_LINEAR_SPRINGS__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGID_BODY;
template<class TV>
class RIGID_LINEAR_SPRINGS:public RIGIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef RIGIDS_FORCES<TV> BASE;
    using BASE::Invalidate_CFL;using BASE::cfl_number;using BASE::rigid_body_collection;
    using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::use_implicit_velocity_independent_forces;
    typedef typename FORCE_ELEMENTS::ITERATOR SEGMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

public:
    SEGMENT_MESH segment_mesh;
    ARRAY<T> youngs_modulus; // units of force (i.e. force per unit strain)
    ARRAY<T> restlength,visual_restlength,current_lengths; // visual restlength corresponds to length between particles; restlength may be larger than this to avoid zero/small restlength
    ARRAY<T> damping; // units of force*time (i.e. force per unit strain rate)
    ARRAY<T> extra_energy;
    ARRAY<T> potential_energy_save;
    ARRAY<T> force_correction;
    ARRAY<TV> previously_applied_forces;
    mutable ARRAY<VECTOR<T,2> > strains_of_segment; // VECTOR<T,2>(strain_rate, strain)
    ARRAY<VECTOR<TV,2> > attachment_radius;
    bool use_kinetic_energy_fix;
    bool relaxation_fraction;
    ARRAY<T> energy_correction_forces;
    bool use_gauss_seidel_in_energy_correction;
    bool allow_kd_direction_flip;

    struct STATE{
        STATE()
            :coefficient(0)
        {}

        VECTOR<int,2> nodes; // copy of nodes to reduce associativity needed in cache
        T coefficient;
        TV direction;
        VECTOR<TV,2> r;
    };
    ARRAY<STATE> states;
public:
    FORCE_ELEMENTS force_segments;
    RIGID_LINEAR_SPRINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    virtual ~RIGID_LINEAR_SPRINGS();

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {} // Add_Force_Differential always enforces definiteness

//#####################################################################
    void Add_Spring(int body1,int body2,const TV& r1,const TV& r2);
    TV Attachment_Location(int s,int b) const;
    TV Endpoint_Velocity(int s,int b) const;
    T Spring_Length(int s) const;
    void Set_Restlengths();
    void Set_Stiffness(int b,T stiffness);
    void Set_Damping(int b,T damp);
    void Set_Restlength(int b,T length=-1,T visual=-1);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE;
    void Set_Overdamping_Fraction(int b,const T overdamping_fraction=1); // 1 is critically damped
    void Update_Position_Based_State(const T time) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force(ARRAY_VIEW<TWIST<TV> > rigid_F,const STATE& state,const TV& force) const;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    T Average_Restlength() const;
    void Print_Restlength_Statistics() const;
    void Print_Deformation_Statistics() const;
    T Maximum_Compression_Or_Expansion_Fraction(int* index=0) const;
    void Print_All_Energy_Debts(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt) const;
    T Potential_Energy(int s,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Save_Potential_Energy(const T time);
    T Compute_Total_Energy(const T time) const;
    void Compute_Previously_Applied_Forces();
    T Compute_Spring_Work(int s,ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt) const;
    void Compute_Energy_Error(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt);
    void Apply_Energy_Correction_Impulse(const int s,const T force,const T dt);
    void Add_Energy_Correction_Force(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt);
    T Effective_Impulse_Factor(int s,int b) const;
    T Effective_Impulse_Factor(int s) const;
    const RIGID_BODY<TV>& Body(int s,int b) const;
    RIGID_BODY<TV>& Body(int s,int b);
//#####################################################################
};
}
#endif
