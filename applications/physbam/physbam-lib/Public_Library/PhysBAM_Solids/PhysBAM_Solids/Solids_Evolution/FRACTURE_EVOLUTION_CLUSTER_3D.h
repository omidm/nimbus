//#####################################################################
// Copyright 2006, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_EVOLUTION_CLUSTER_3D
//#####################################################################
#ifndef __FRACTURE_EVOLUTION_CLUSTER_3D__
#define __FRACTURE_EVOLUTION_CLUSTER_3D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/PLASTICITY_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SOLIDS_FORCES_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION.h>
namespace PhysBAM{

template<class TV> class SOLIDS_PARAMETERS;
template<class TV>
class FRACTURE_EVOLUTION_CLUSTER_3D:public FRACTURE_EVOLUTION<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename EMBEDDING_POLICY<TV,3>::EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
    typedef typename EMBEDDING_POLICY<TV,3>::EMBEDDING T_EMBEDDING;
    typedef typename EMBEDDING_POLICY<TV,3>::EMBEDDED_MATERIAL_SURFACE T_EMBEDDED_MATERIAL_SURFACE;

    using FRACTURE_EVOLUTION<TV>::push_out;using FRACTURE_EVOLUTION<TV>::perturb_amount_for_collision_freeness;using FRACTURE_EVOLUTION<TV>::fractured_after_rebuild_topology;
    using FRACTURE_EVOLUTION<TV>::collision_iterations;

public:
    SOLIDS_PARAMETERS<VECTOR<T,3> >& solids_parameters;
//    ARRAY<RIGID_BODY_CLUSTER_3D<T>*> clusters_to_check;
    ARRAY<PAIR<int,int> > new_clusters; // where the pair is <new id,old id>
    T previous_time;

    FRACTURE_EVOLUTION_CLUSTER_3D(SOLIDS_PARAMETERS<VECTOR<T,3> >& solids_parameters_input)
        :FRACTURE_EVOLUTION<TV>(),solids_parameters(solids_parameters_input),previous_time(0)
    {}

    virtual ~FRACTURE_EVOLUTION_CLUSTER_3D()
    {}

    void Reset_New_Cluster_List()
    {new_clusters.Clean_Memory();}

    bool Add_Scripted_Cuts(const T time) PHYSBAM_OVERRIDE {return false;} // TEMPORARY

//#####################################################################
    void Initialize_Bodies() PHYSBAM_OVERRIDE {}
    void Initialize_Self_Collision(){}
    void Reinitialize_Bodies(){}
    void Rebuild_Topology(){}
    int Fracture_Where_High_Stress(const T small_number=1e-4){return 0;}
    TV Spatial_Fracture_Bias_Direction(const int t,const T small_number) const{return TV();}
    int Adjust_Nodes_For_Segment_Triangle_Intersections(T threshhold=1e-2){return 0;}
    void Process_Cluster_Fracture(const T dt,const T time,SOLIDS_EVOLUTION<TV>* rigid_body_evolution,SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callback,const bool force);
//#####################################################################
};
}
#endif
