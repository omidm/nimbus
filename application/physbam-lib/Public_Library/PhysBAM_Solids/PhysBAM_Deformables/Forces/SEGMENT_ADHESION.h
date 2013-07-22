//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_ADHESION
//#####################################################################
#ifndef __SEGMENT_ADHESION__
#define __SEGMENT_ADHESION__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/HAIR_ID.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{
class SEGMENT_MESH;

template<class TV>
class SEGMENT_ADHESION_SPRING_STATE
{
    typedef typename TV::SCALAR T;
public:
    VECTOR<int,4> nodes;
    VECTOR<T,2> weights;
    TV normal;
    T distance;
    T damping;
    bool external;    

    SEGMENT_ADHESION_SPRING_STATE()
        :distance(0),damping(0),external(false)
    {}
};

template<class TV>
class SEGMENT_ADHESION:public DEFORMABLES_FORCES<TV>,public SPRINGS_TAG
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT T_SEGMENTED_CURVE;
    enum WORKAROUND {d=TV::m};
    MPI_SOLIDS<TV>* mpi_solids;
public:
    typedef SEGMENT_ADHESION_SPRING_STATE<TV> SPRING_STATE;
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    SEGMENT_MESH& mesh;
    T_SEGMENTED_CURVE curve;
    T on_distance;
    T off_distance;
    T overdamping_fraction;
    T youngs_modulus;
    T restlength;
    int max_connections;
    HASHTABLE<VECTOR<int,2> > existing_pairs;
    ARRAY<HAIR_ID>& particle_to_spring_id;
    
    ARRAY<int> internal_segment_indices;
    ARRAY<int> external_segment_indices;
    SEGMENT_MESH internal_mesh,external_mesh;
    T_SEGMENTED_CURVE internal_curve,external_curve;

    typedef HASHTABLE<VECTOR<int,2>,SPRING_STATE> T_SPRING_HASH;
    T_SPRING_HASH *springs; // Internal springs and external springs owned by this processor
    ARRAY<VECTOR<int,2> > external_spring_segments;
    ARRAY<SPRING_STATE> internal_springs; // cache for Add_Velocity_Dependent_Forces
    ARRAY<SPRING_STATE> external_springs;
    
    ARRAY<PAIR<ARRAY<PAIR<T,int> >,bool> > segments_with_springs;
    HASHTABLE<VECTOR<int,4> >& intersecting_edge_edge_pairs;
    HASHTABLE<VECTOR<int,4> > default_intersecting_edge_edge_pairs;

    SEGMENT_ADHESION(PARTICLES<TV>& particles,SEGMENT_MESH& mesh,ARRAY<HAIR_ID>& particle_to_spring_id,HASHTABLE<VECTOR<int,4> >& intersecting_edge_edge_pairs);
    SEGMENT_ADHESION(PARTICLES<TV>& particles,SEGMENT_MESH& mesh,ARRAY<HAIR_ID>& particle_to_spring_id);

    virtual ~SEGMENT_ADHESION()
    {
        delete springs;
    }

//#####################################################################
    void Update_Hierarchy();
    void Update_Springs(const bool search_hierarchy);
    void Update_Collisions_List();
    void Update_Partitions(bool restart,MPI_SOLIDS<TV>* mpi_solids,const std::string output_directory);
    void Set_Parameters(const T youngs_modulus_input,const T overdamping_fraction_input,const T on_distance,const T off_distance, const int max_connections_input);
    void Set_Restlength(const T restlength);
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
    void Read_State(STREAM_TYPE type,const std::string& filename);
//#####################################################################
};
}
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Forces_And_Torques/READ_WRITE_SEGMENT_ADHESION_SPRING_STATE.h>
#endif
