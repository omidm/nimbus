//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_BENDING_ELEMENTS
//#####################################################################
#ifndef __SEGMENT_BENDING_ELEMENTS__
#define __SEGMENT_BENDING_ELEMENTS__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class T_input>
class SEGMENT_BENDING_ELEMENTS:public DEFORMABLES_FORCES<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;
    using BASE::use_rest_state_for_strain_rate;using BASE::Limit_Time_Step_By_Strain_Rate;
    typedef typename FORCE_ELEMENTS::ITERATOR TRIPLE_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    ARRAY<VECTOR<int,3> > bending_triples; // for each triple, the bending point is vertex 2.
    ARRAY<T> length_scale;
    ARRAY<T> stiffness;
    ARRAY<T> sine_half_rest_angle;
    ARRAY<T> damping;
private:
    ARRAY<T> elastic_s,damping_coefficient;
    ARRAY<VECTOR<TV,2> > force_directions; // directions for nodes 1 and 3 (2 given by subtraction)
    FORCE_ELEMENTS force_triples;
public:

    SEGMENT_BENDING_ELEMENTS(PARTICLES<TV>& particles)
        :DEFORMABLES_FORCES<TV>(particles)
    {}

private:
    T Sine_Half_Angle_Between(const TV& n1,const TV& n2) const
    {T sine_half_psi=(T).5*(n1-n2).Magnitude();
    if(TV::Cross_Product(n1,n2).x<0) sine_half_psi=-sine_half_psi;
    return sine_half_psi;}
public:

    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE
    {} // TODO: Implement me?

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Set_Triples_From_Segment_Mesh(SEGMENT_MESH& mesh);
    void Set_Constants_From_Particles(const T material_stiffness,const T material_damping);
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T> SEGMENT_BENDING_ELEMENTS<T>*
Create_Bending_Elements(PARTICLES<VECTOR<T,2> >& particles,SEGMENT_MESH& mesh,const T stiffness=1e-3,const T damping=1,
    const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool verbose=true)
{
    SEGMENT_BENDING_ELEMENTS<T>* bend=new SEGMENT_BENDING_ELEMENTS<T>(particles);
    bend->Set_Triples_From_Segment_Mesh(mesh);
    bend->Set_Constants_From_Particles(stiffness,damping);
    bend->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    return bend;
}

template<class T> SEGMENT_BENDING_ELEMENTS<T>*
Create_Bending_Elements(SEGMENTED_CURVE<VECTOR<T,2> >& segmented_curve,const T stiffness=1e-3,
    const T damping=1,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool verbose=true)
{
    return Create_Bending_Elements(dynamic_cast<PARTICLES<VECTOR<T,2> >&>(segmented_curve.particles),segmented_curve.mesh,stiffness,damping,limit_time_step_by_strain_rate,
        max_strain_per_time_step,verbose);
}

}
#endif
