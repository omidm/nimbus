//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_POINTWISE_FORCE
//#####################################################################
#ifndef __RIGID_POINTWISE_FORCE__
#define __RIGID_POINTWISE_FORCE__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
namespace PhysBAM{

template<class TV>
class RIGID_POINTWISE_FORCE:public RIGIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef RIGIDS_FORCES<TV> BASE;
    using BASE::rigid_body_collection;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;
protected:
    ARRAY<int> *influenced_rigid_body_particles;
    bool need_destroy_influenced_rigid_body_particles;
    bool influence_all_rigid_body_particles;
public:
    FORCE_ELEMENTS force_rigid_body_particles;

    RIGID_POINTWISE_FORCE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_rigid_body_particles_input)
        :RIGIDS_FORCES<TV>(rigid_body_collection_input),influenced_rigid_body_particles(influenced_rigid_body_particles_input),
        need_destroy_influenced_rigid_body_particles(false),influence_all_rigid_body_particles(false)
    {}

    RIGID_POINTWISE_FORCE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_rigid_body_particles_input)
        :RIGIDS_FORCES<TV>(rigid_body_collection_input),influenced_rigid_body_particles(0),need_destroy_influenced_rigid_body_particles(true),
        influence_all_rigid_body_particles(influence_all_rigid_body_particles_input)
    {}

    virtual ~RIGID_POINTWISE_FORCE()
    {if(need_destroy_influenced_rigid_body_particles) delete influenced_rigid_body_particles;}

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE
    {}

    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE
    {return FLT_MAX;}

    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> rigid_frequency) PHYSBAM_OVERRIDE
    {}

protected:
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE;

    template<class T_ARRAY>
    ARRAY<int> Get_Rigid_Body_Particle_List(const T_ARRAY& array);
//#####################################################################
};
}
#endif
