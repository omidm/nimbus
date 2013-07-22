//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_GEOMETRY_STANDARD_TESTS
//#####################################################################
#ifndef __RIGIDS_GEOMETRY_STANDARD_TESTS__
#define __RIGIDS_GEOMETRY_STANDARD_TESTS__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
namespace PhysBAM{

template<class TV> class EXAMPLE;
template<class TV> class RIGID_GEOMETRY;
template<class TV> class RIGID_GEOMETRY_COLLECTION;

template<class TV>
class RIGIDS_GEOMETRY_STANDARD_TESTS
{
    typedef typename TV::SCALAR T;
public:
    EXAMPLE<TV>& example;
    RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection;

    RIGIDS_GEOMETRY_STANDARD_TESTS(EXAMPLE<TV>& example_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_body_collection_input);

//#####################################################################
    RIGID_GEOMETRY<TV>& Add_Rigid_Body(const std::string& name,const T scaling_factor,const T friction,const bool read_implicit=true,const bool always_read_object=false);
    RIGID_GEOMETRY<TV>& Add_Ground(const T friction=(T).3,const T height=0,const T scale=(T)1);
//#####################################################################
};
}
#endif
