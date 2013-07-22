//#####################################################################
// Copyright 2006-2008, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINDING_SPRINGS
//#####################################################################
#ifndef __BINDING_SPRINGS__
#define __BINDING_SPRINGS__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/IMPLICIT_ZERO_LENGTH_SPRINGS.h>
namespace PhysBAM{

template<class TV>
class BINDING_SPRINGS:public IMPLICIT_ZERO_LENGTH_SPRINGS<TV>
{
    typedef typename TV::SCALAR T;
public:
    BINDING_SPRINGS(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input)
        :IMPLICIT_ZERO_LENGTH_SPRINGS<TV>(particles,segment_mesh_input)
    {}

    ~BINDING_SPRINGS()
    {}

    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}
};

template<class TV> BINDING_SPRINGS<TV>*
Create_Edge_Binding_Springs(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const typename TV::SCALAR stiffness=2e3,
    const typename TV::SCALAR overdamping_fraction=1,const bool verbose=true)
{
    BINDING_SPRINGS<TV>* ls=new BINDING_SPRINGS<TV>(particles,segment_mesh);
    ls->Set_Stiffness(stiffness);
    ls->Set_Overdamping_Fraction(overdamping_fraction);
    return ls;
}

}
#endif
