//#####################################################################
// Copyright 2006, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_EVOLUTION_S3D
//#####################################################################
#ifndef __FRACTURE_EVOLUTION_S3D__
#define __FRACTURE_EVOLUTION_S3D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/PLASTICITY_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SOLIDS_FORCES_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION.h>
namespace PhysBAM{

template<class TV> class SOLIDS_PARAMETERS;
template<class T_input>
class FRACTURE_EVOLUTION_S3D:public FRACTURE_EVOLUTION<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef typename EMBEDDING_POLICY<TV,2>::EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
    typedef typename EMBEDDING_POLICY<TV,2>::EMBEDDING T_EMBEDDING;
    typedef typename EMBEDDING_POLICY<TV,2>::EMBEDDED_MATERIAL_SURFACE T_EMBEDDED_MATERIAL_SURFACE;

    using FRACTURE_EVOLUTION<TV>::push_out;using FRACTURE_EVOLUTION<TV>::perturb_amount_for_collision_freeness;

public:
    SOLIDS_PARAMETERS<TV>& solids_parameters;
    FRACTURE_OBJECT<TV,2>* fracture_object;
    PLASTICITY_MODEL<T,2>* plasticity_model;
    ARRAY<VECTOR<T,3> > spatial_fracture_bias_direction;        // bias direction in spacial coordinates (**need to be converted to world space for diagonalized**)

    FRACTURE_EVOLUTION_S3D(SOLIDS_PARAMETERS<TV>& solids_parameters_input)
        :FRACTURE_EVOLUTION<TV>(),solids_parameters(solids_parameters_input),fracture_object(0),plasticity_model(0)
    {}

    virtual ~FRACTURE_EVOLUTION_S3D()
    {}

    bool Add_Scripted_Cuts(const T time) PHYSBAM_OVERRIDE {return false;} // TEMPORARY

//#####################################################################
    void Initialize_Bodies() PHYSBAM_OVERRIDE;
    void Initialize_Self_Collision();
    void Reinitialize_Bodies();
    void Rebuild_Topology();
    int Fracture_Where_High_Stress(const T small_number=1e-4);
    TV Spatial_Fracture_Bias_Direction(const int t,const T small_number) const;
//#####################################################################
};
}
#endif
