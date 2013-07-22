//#####################################################################
// Copyright 2011 Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDER_WORLD_MOUSE_LISTENER
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDER_WORLD_MOUSE_LISTENER__
#define __OPTIX_RENDER_WORLD_MOUSE_LISTENER__

namespace PhysBAM{

class OPTIX_RENDER_WORLD_MOUSE_LISTENER {
public:
    virtual void Handle_Click(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed) = 0;
    virtual void Handle_Drag(int x,int y) = 0;
};
}
#endif
#endif

