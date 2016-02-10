//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_COLLISION_MANAGER_HASH
//#####################################################################
#ifndef __RIGID_BODY_COLLISION_MANAGER_HASH__
#define __RIGID_BODY_COLLISION_MANAGER_HASH__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
namespace PhysBAM{

class RIGID_BODY_COLLISION_MANAGER_HASH:public RIGID_BODY_COLLISION_MANAGER
{
public:
    bool default_collide;
    typedef PAIR<int,int> T_PAIR;
    HASHTABLE<T_PAIR> hash; // if default_collide=0 then holds colliders, else holds non-colliders

    RIGID_BODY_COLLISION_MANAGER_HASH()
        :default_collide(true)
    {}

    virtual ~RIGID_BODY_COLLISION_MANAGER_HASH()
    {}
   
    bool Body_Collides_With_The_Other(int rigid_body_id_1,int rigid_body_id_2) const PHYSBAM_OVERRIDE
    {return default_collide^hash.Contains(T_PAIR(rigid_body_id_1,rigid_body_id_2));}

    bool Either_Body_Collides_With_The_Other(int rigid_body_id_1,int rigid_body_id_2) const PHYSBAM_OVERRIDE
    {return Body_Collides_With_The_Other(rigid_body_id_1,rigid_body_id_2) || Body_Collides_With_The_Other(rigid_body_id_2,rigid_body_id_1);}

//#####################################################################
};
}
#endif
