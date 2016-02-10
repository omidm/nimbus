//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_MATERIAL_SPACE_FORCES
//#####################################################################
#ifndef __BW_MATERIAL_SPACE_FORCES__
#define __BW_MATERIAL_SPACE_FORCES__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_FORCES.h>
namespace PhysBAM{

class TRIANGLE_MESH;

template<class TV,int d>
class BW_MATERIAL_SPACE_FORCES:public BW_FORCES<TV,d,3>
{
    typedef typename TV::SCALAR T;
public:
    typedef BW_FORCES<TV,d,3> BASE;
    using BASE::particles;using BASE::states;using BASE::force_simplices;using BASE::triangle_mesh;
    typedef typename FORCE_ELEMENTS::ITERATOR TRIANGLE_ITERATOR;

protected:
    struct MATERIAL_FORCE_STATE{
        MATERIAL_FORCE_STATE()
            :b_u(1),b_v(1)
        {}

        // These don't change across the entire sim
        VECTOR<T,2> delta_u;
        VECTOR<T,2> delta_v;
        T rest_state_triangle_area;
        T b_u;
        T b_v;

        // These are set each time UPBS is called
        T denom;
        TV w_u;
        TV w_v;
        T w_u_magnitude;
        T w_v_magnitude;

        TV dwu_dx;
        TV dwv_dx;
    };
    ARRAY<MATERIAL_FORCE_STATE> material_force_states;
public:
    BW_MATERIAL_SPACE_FORCES(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input);

    virtual ~BW_MATERIAL_SPACE_FORCES();

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Compute_UV_Deformation(const int c);
//#####################################################################
};
}
#endif
