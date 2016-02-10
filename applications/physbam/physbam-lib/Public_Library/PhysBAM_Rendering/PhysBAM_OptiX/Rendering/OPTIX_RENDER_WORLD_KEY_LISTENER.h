//#####################################################################
// Copyright 2011 Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDER_WORLD_KEY_LISTENER
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDER_WORLD_KEY_LISTENER__
#define __OPTIX_RENDER_WORLD_KEY_LISTENER__

namespace PhysBAM{

class OPTIX_RENDER_WORLD_KEY_LISTENER {
public:
    virtual void KeyPressed(unsigned char key, int x, int y) {}
};
}
#endif
#endif
