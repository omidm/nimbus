//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_MOUSE_HANDLER__
#define __OPENGL_MOUSE_HANDLER__
#include <iostream>
#include <string>
namespace PhysBAM{

class OPENGL_MOUSE_HANDLER
{
public:
    virtual ~OPENGL_MOUSE_HANDLER() {}
    virtual bool Handle_Click(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed)=0;
    virtual bool Handle_Drag(int x,int y){return false;}
//#####################################################################
};
}
#endif

