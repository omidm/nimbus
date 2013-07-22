//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_COLLISION_MANAGER
//#####################################################################
#ifndef __RIGID_BODY_COLLISION_MANAGER__
#define __RIGID_BODY_COLLISION_MANAGER__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

class RIGID_BODY_COLLISION_MANAGER:public NONCOPYABLE
{
public:
    virtual ~RIGID_BODY_COLLISION_MANAGER()
    {}
    
    virtual bool Body_Collides_With_The_Other(int rigid_body_id_1,int rigid_body_id_2) const=0;
    virtual bool Either_Body_Collides_With_The_Other(int rigid_body_id_1,int rigid_body_id_2) const=0;

//#####################################################################
};
}
#endif
