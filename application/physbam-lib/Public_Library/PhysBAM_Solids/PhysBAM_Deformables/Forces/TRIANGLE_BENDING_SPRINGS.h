//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Unnur Gretarsdottir, Geoffrey Irving, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_BENDING_SPRINGS
//#####################################################################
#ifndef __TRIANGLE_BENDING_SPRINGS__
#define __TRIANGLE_BENDING_SPRINGS__

#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
namespace PhysBAM{

class TRIANGLE_MESH;

template<class T_input>
class TRIANGLE_BENDING_SPRINGS:public LINEAR_SPRINGS<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef LINEAR_SPRINGS<TV> BASE;
    using LINEAR_SPRINGS<TV>::Set_Stiffness;using LINEAR_SPRINGS<TV>::Set_Damping;using LINEAR_SPRINGS<TV>::segment_mesh;

    SEGMENT_MESH bending_segment_mesh;

    TRIANGLE_BENDING_SPRINGS(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh,const bool implicit=false);

    virtual ~TRIANGLE_BENDING_SPRINGS();

    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const PHYSBAM_OVERRIDE
    {if(force_name.empty()) BASE::Add_Force_Data(force_data_list,"TRIANGLE_BENDING_SPRINGS");
    else BASE::Add_Force_Data(force_data_list,force_name);}

//#####################################################################
    void Initialize(TRIANGLE_MESH& triangle_mesh);
//#####################################################################
};

template<class T> TRIANGLE_BENDING_SPRINGS<T>*
Create_Bending_Springs(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& triangle_mesh,const T stiffness=2/(1+sqrt((T)2)),
    const T overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,
    const T restlength_enlargement_fraction=0,const bool verbose=true,const bool implicit=false);

template<class T> TRIANGLE_BENDING_SPRINGS<T>*
Create_Bending_Springs(TRIANGULATED_SURFACE<T>& triangulated_surface,
    const T stiffness=2/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true,const bool implicit=false);

}
#endif
