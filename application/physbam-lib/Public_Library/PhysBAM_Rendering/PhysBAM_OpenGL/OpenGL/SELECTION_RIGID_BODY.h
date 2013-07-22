//#####################################################################
// Copyright 2011, Valeria Nikolaenko
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SELECTION_RIGID_BODY
//#####################################################################
#ifndef __SELECTION_RIGID_BODY__
#define __SELECTION_RIGID_BODY__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <climits>

namespace PhysBAM {

template<class TV>
class SELECTION_RIGID_BODY {
public:
    int body_id;
    virtual ~SELECTION_RIGID_BODY(){}
    virtual TV Get_Object_Location() = 0;
    virtual bool Was_Dragged() = 0;
    virtual int Get_Body_ID() = 0;
};

}

#endif


