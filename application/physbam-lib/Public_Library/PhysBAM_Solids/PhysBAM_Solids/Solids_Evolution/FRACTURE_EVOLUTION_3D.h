//#####################################################################
// Copyright 2006, Frank Losasso, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_EVOLUTION_3D
//#####################################################################
#ifndef __FRACTURE_EVOLUTION_3D__
#define __FRACTURE_EVOLUTION_3D__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/PLASTICITY_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SOLIDS_FORCES_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION.h>
namespace PhysBAM{

template<class TV> class RIGID_DEFORMABLE_EVOLUTION_OLD;
template<class TV> class FRACTURE_CALLBACKS;
template<class TV> class BINDING_LIST;
template<class TV> class SOLIDS_PARAMETERS;

template<class T_input>
class FRACTURE_EVOLUTION_3D:public FRACTURE_EVOLUTION<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef typename EMBEDDING_POLICY<TV,3>::EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
    typedef typename EMBEDDING_POLICY<TV,3>::EMBEDDING T_EMBEDDING;
    typedef typename EMBEDDING_POLICY<TV,3>::EMBEDDED_MATERIAL_SURFACE T_EMBEDDED_MATERIAL_SURFACE;
    typedef typename MATRIX_POLICY<TV>::DIAGONAL_MATRIX T_DIAGONAL_MATRIX;

public:
    using FRACTURE_EVOLUTION<TV>::push_out;using FRACTURE_EVOLUTION<TV>::perturb_amount_for_collision_freeness;using FRACTURE_EVOLUTION<TV>::fractured_after_rebuild_topology;

    SOLIDS_PARAMETERS<TV>& solids_parameters;
    FRACTURE_OBJECT<TV,3>* fracture_object;
    PLASTICITY_MODEL<T,3>* plasticity_model;
    PARTICLES<TV> rigid_body_deformable_body_particles;
    ARRAY<int> particle_to_rigid_body_id;
    ARRAY<int> deformable_to_rigid_particles;
    ARRAY<int> rigid_bodies_with_impulse;
    bool create_new_rigid_bodies;
    bool rigid_fracture;
    bool save_state;

public:
    FRACTURE_EVOLUTION_3D(SOLIDS_PARAMETERS<TV>& solids_parameters_input,const bool rigid_fracture_input=false)
        :FRACTURE_EVOLUTION<TV>(),solids_parameters(solids_parameters_input),fracture_object(0),plasticity_model(0),create_new_rigid_bodies(true),rigid_fracture(rigid_fracture_input),save_state(true)
    {}

    virtual ~FRACTURE_EVOLUTION_3D()
    {}

    bool Add_Scripted_Cuts(const T time) PHYSBAM_OVERRIDE {return false;} // TEMPORARY

//#####################################################################
    void Initialize_Bodies() PHYSBAM_OVERRIDE;
    void Initialize_Self_Collision();
    void Reinitialize_Bodies();
    void Rebuild_Topology();
    int Fracture_Where_High_Stress(const T small_number=1e-4);
    int Rigid_Fracture_Where_High_Stress(const T small_number=1e-4);
    TV Spatial_Fracture_Bias_Direction(const int t,const T small_number) const;
    int Adjust_Nodes_For_Segment_Triangle_Intersections(T threshhold=1e-2);
    void Create_Rigid_Body_Fracture_Object(const TV velocity,const TV angular_velocity,const ARRAY<TV>* seed_positions,const ARRAY<T>* seed_weakness_multipliers,const FRACTURE_CALLBACKS<TV>* fracture_callbacks);
    void Create_New_Rigid_Bodies_From_Fracture(ARRAY<int>& map_to_old_particles);
    // TODO: this is being pulled out of Solids Standard tests until after siggraph--don't want to mess up solids people
    void Substitute_Soft_Bindings_For_Embedded_Nodes(T_EMBEDDED_MATERIAL_SURFACE& object,const BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,HASHTABLE<int,int>* persistent_soft_bindings=0);
    void Process_Rigid_Fracture(const T dt,const T time,SOLIDS_EVOLUTION<TV>* rigid_deformable_evolution_old,SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks);
    bool Run_Quasistatics_And_Fracture(const T time,const int max_number_of_fracture_iterations);
//#####################################################################
};
}
#endif
