//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Unnur Gretarsdottir, Jon Gretarsson, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_BENDING_SPRINGS
//#####################################################################
#ifndef __SEGMENT_BENDING_SPRINGS__
#define __SEGMENT_BENDING_SPRINGS__

#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
namespace PhysBAM{

template<class TV>
class SEGMENT_BENDING_SPRINGS:public LINEAR_SPRINGS<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef LINEAR_SPRINGS<TV> BASE;
    using LINEAR_SPRINGS<TV>::Set_Stiffness;using LINEAR_SPRINGS<TV>::Set_Damping;using LINEAR_SPRINGS<TV>::segment_mesh;

    SEGMENT_MESH bending_segment_mesh;

    SEGMENT_BENDING_SPRINGS(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const bool implicit=false);

    virtual ~SEGMENT_BENDING_SPRINGS();

//#####################################################################
    void Initialize(SEGMENT_MESH& segment_mesh);
//#####################################################################
};

template<class TV> SEGMENT_BENDING_SPRINGS<TV>*
Create_Segment_Bending_Springs(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const typename TV::SCALAR stiffness=2/(1+sqrt((typename TV::SCALAR)2)),
    const typename TV::SCALAR overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const typename TV::SCALAR max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,
    const typename TV::SCALAR restlength_enlargement_fraction=0,const bool verbose=true,const bool implicit=false);

template<class TV> SEGMENT_BENDING_SPRINGS<TV>*
Create_Segment_Bending_Springs(SEGMENTED_CURVE<TV>& segmented_curve,const typename TV::SCALAR stiffness=2/(1+sqrt((typename TV::SCALAR)2)),
    const typename TV::SCALAR overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const typename TV::SCALAR max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const typename TV::SCALAR restlength_enlargement_fraction=0,const bool verbose=true,const bool implicit=false);

}
#endif
