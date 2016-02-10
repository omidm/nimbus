//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_SHEAR_FORCES
//#####################################################################
#ifndef __BW_SHEAR_FORCES__
#define __BW_SHEAR_FORCES__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_MATERIAL_SPACE_FORCES.h>
namespace PhysBAM{

class TRIANGLE_MESH;

template<class TV>
class BW_SHEAR_FORCES:public BW_MATERIAL_SPACE_FORCES<TV,1>
{
    typedef typename TV::SCALAR T;
public:
    typedef BW_MATERIAL_SPACE_FORCES<TV,1> BASE;
    using BASE::particles;using BASE::force_simplices;using BASE::states;using BASE::material_force_states;
    using BASE::Potential_Energy;using BASE::Compute_UV_Deformation;
    typedef typename FORCE_ELEMENTS::ITERATOR TRIANGLE_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    BW_SHEAR_FORCES(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input);

    virtual ~BW_SHEAR_FORCES();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    T Potential_Energy(int s,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class TV> BW_SHEAR_FORCES<TV>*
Create_BW_Shear_Force(PARTICLES<TV>& particles,TRIANGLE_MESH& segment_mesh,const typename TV::SCALAR stiffness_coefficient_input,const typename TV::SCALAR damping_coefficient_input);
}
#endif
