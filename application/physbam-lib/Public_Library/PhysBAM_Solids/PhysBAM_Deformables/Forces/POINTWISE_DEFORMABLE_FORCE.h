//#####################################################################
// Copyright 2007-2008, Michael Lentine, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINTWISE_DEFORMABLE_FORCE
//#####################################################################
#ifndef __POINTWISE_DEFORMABLE_FORCE__
#define __POINTWISE_DEFORMABLE_FORCE__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class POINTWISE_DEFORMABLE_FORCE:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;
protected:
    ARRAY<int> *influenced_particles;
    bool need_destroy_influenced_particles;
    bool influence_all_particles;
    MPI_SOLIDS<TV>* mpi_solids;
    ARRAY<bool> is_simulated;

public:
    FORCE_ELEMENTS force_particles;

    POINTWISE_DEFORMABLE_FORCE(PARTICLES<TV>& particles_input,ARRAY<int>* influenced_particles_input)
        :DEFORMABLES_FORCES<TV>(particles_input),influenced_particles(influenced_particles_input),
        need_destroy_influenced_particles(false),influence_all_particles(false),mpi_solids(0)
    {}

    POINTWISE_DEFORMABLE_FORCE(PARTICLES<TV>& particles_input,const bool influence_all_particles_input)
        :DEFORMABLES_FORCES<TV>(particles_input),influenced_particles(0),need_destroy_influenced_particles(true),
        influence_all_particles(influence_all_particles_input),mpi_solids(0)
    {}

    template<class T_MESH>
    POINTWISE_DEFORMABLE_FORCE(PARTICLES<TV>& particles_input,const T_MESH& mesh)
        :DEFORMABLES_FORCES<TV>(particles_input),influenced_particles(new ARRAY<int>),
        need_destroy_influenced_particles(true),influence_all_particles(false)
    {
        mesh.elements.Flattened().Get_Unique(*influenced_particles);
    }

    virtual ~POINTWISE_DEFORMABLE_FORCE()
    {if(need_destroy_influenced_particles) delete influenced_particles;}

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE
    {}

    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;

    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE
    {return FLT_MAX;}

    void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE
    {}

protected:
    template<class T_ARRAY>
    T_ARRAY Get_Particle_List(const T_ARRAY& array)
    {return array;}
//#####################################################################
};
}
#endif
