//#####################################################################
// Copyright 2006, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_EVOLUTION
//#####################################################################
#ifndef __FRACTURE_EVOLUTION__
#define __FRACTURE_EVOLUTION__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class TV> class SOLIDS_EVOLUTION;
template<class TV> class SOLIDS_EVOLUTION_CALLBACKS;

template<class TV>
class FRACTURE_EVOLUTION
{
    typedef typename TV::SCALAR T;
public:
    bool perform_fracture;
    T perturb_amount_for_collision_freeness;
    bool fractured_after_rebuild_topology;
    int substeps_before_rebuild;
    bool push_out,check_for_fracture_only_on_frame_boundaries;
    int collision_iterations;

    FRACTURE_EVOLUTION()
        :perform_fracture(true),perturb_amount_for_collision_freeness((T)1e-6),fractured_after_rebuild_topology(false),substeps_before_rebuild(1),push_out(false),
        check_for_fracture_only_on_frame_boundaries(false),collision_iterations(5)
    {}

    virtual ~FRACTURE_EVOLUTION()
    {}

    void Set_Collision_Iterations(const int collision_iterations_input=5)
    {collision_iterations=collision_iterations_input;}

    virtual bool Add_Scripted_Cuts(const T time){return false;} // TEMPORARY

//#####################################################################
    virtual void Initialize_Bodies(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();};
    virtual void Initialize_Self_Collision(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();};
    virtual void Reinitialize_Bodies(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();};
    virtual void Rebuild_Topology(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();};
    void Postprocess_Solids_Substep(const T time,const int substep);
    void Postprocess_Frame(const int frame);
    virtual int Fracture_Where_High_Stress(const T small_number=1e-4){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;};
    virtual TV Spatial_Fracture_Bias_Direction(const int t,const T small_number) const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return TV();};
    virtual int Adjust_Nodes_For_Segment_Triangle_Intersections(T threshhold=1e-2){return 0;};
    virtual void Preprocess_Fracture(const T dt,const T time,SOLIDS_EVOLUTION<TV>* rigid_body_evolution,SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callback)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();};
    virtual void Process_Rigid_Fracture(const T dt,const T time,SOLIDS_EVOLUTION<TV>* rigid_body_evolution,SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();};
//#####################################################################
};
}
#endif
