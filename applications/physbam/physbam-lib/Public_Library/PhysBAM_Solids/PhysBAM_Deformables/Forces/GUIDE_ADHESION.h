//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GUIDE_ADHESION
//#####################################################################
#ifndef __GUIDE_ADHESION__
#define __GUIDE_ADHESION__

#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/HAIR_ID.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{
class SEGMENT_MESH;

template<class TV>
class GUIDE_ADHESION_SPRING_STATE
{
    typedef typename TV::SCALAR T;
public:
    VECTOR<int,4> nodes;
    VECTOR<T,2> weights;
    TV normal;
    T distance;
    T damping;
    T restlength;

    GUIDE_ADHESION_SPRING_STATE()
        :distance(0),damping(0)
    {}
};

template<class TV>
class GUIDE_ADHESION:public DEFORMABLES_FORCES<TV>,public SPRINGS_TAG
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT T_SEGMENTED_CURVE;

public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;
    typedef GUIDE_ADHESION_SPRING_STATE<TV> SPRING_STATE;

    SEGMENT_MESH& mesh;
    SEGMENT_MESH& guide_mesh;
    T_SEGMENTED_CURVE curve;
    T_SEGMENTED_CURVE guide_curve;
    T overdamping_fraction;
    T youngs_modulus;
    T thickness;
    int max_connections;
    ARRAY<HAIR_ID>& particle_to_spring_id;
    ARRAY<int,HAIR_ID>& roots;
    const ARRAY_VIEW<T>& one_over_effective_mass;

    typedef HASHTABLE<VECTOR<int,2>,SPRING_STATE> T_SPRING_HASH;
    T_SPRING_HASH* springs;
    ARRAY<int> segments_with_springs;

    GUIDE_ADHESION(PARTICLES<TV>& particles,SEGMENT_MESH& mesh_input,SEGMENT_MESH& guide_mesh_input,
        ARRAY<HAIR_ID>& particle_to_spring_id_input,ARRAY<int,HAIR_ID>& roots_input);

    virtual ~GUIDE_ADHESION();

//#####################################################################
    void Update_Springs(const bool search_hierarchy);
    void Set_Parameters(const T youngs_modulus_input,const T overdamping_fraction_input,const T thickness_input,const int max_connections_input);
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    virtual void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    virtual void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    virtual void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    virtual void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    virtual void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    virtual void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    virtual T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    virtual void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    void Write_State(STREAM_TYPE type,const std::string& filename);
//#####################################################################
};
}
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Forces_And_Torques/READ_WRITE_GUIDE_ADHESION_SPRING_STATE.h>
#endif
