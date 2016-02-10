//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_BENDING_FORCES
//#####################################################################
#ifndef __BW_BENDING_FORCES__
#define __BW_BENDING_FORCES__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_FORCES.h>
namespace PhysBAM{

class TRIANGLE_MESH;

template<class TV>
class BW_BENDING_FORCES:public BW_FORCES<TV,1,4>
{
    typedef typename TV::SCALAR T;
public:
    typedef BW_FORCES<TV,1,4> BASE;
    using BASE::particles;using BASE::states;using BASE::force_simplices;using BASE::triangle_mesh;
    typedef typename FORCE_ELEMENTS::ITERATOR CONSTRAINT_ITERATOR;

protected:
    struct BENDING_STATE{
        BENDING_STATE() {}

        // These are set each time UPBS is called
        TV n_a;
        TV n_b;
        T n_a_magnitude;
        T n_b_magnitude;
        TV e;
        T e_magnitude;
    };
public:
    ARRAY<BENDING_STATE> bending_states;
    ARRAY<VECTOR<int,4> > constraint_particles; // spring is shortest line between segment with particles (1,2) and segment with particles (3,4)
    MATRIX<VECTOR<VECTOR<T,3>,3>,4> dq_a_i_j_t;
    MATRIX<VECTOR<VECTOR<T,3>,3>,4> dq_b_i_j_t;
    bool assume_constant_normal_length;

    BW_BENDING_FORCES(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input);

    virtual ~BW_BENDING_FORCES();

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    T Potential_Energy(int s,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class TV> BW_BENDING_FORCES<TV>*
Create_BW_Bending_Force(PARTICLES<TV>& particles,TRIANGLE_MESH& segment_mesh,const typename TV::SCALAR stiffness_coefficient_input,const typename TV::SCALAR damping_coefficient_input);
}
#endif
