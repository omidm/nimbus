//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MOUSE_CAMERA_MANIPULATOR
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MOUSE_CAMERA_MANIPULATOR__
#define __OPTIX_MOUSE_CAMERA_MANIPULATOR__

#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_CAMERA.h>
#include <GL/glut.h>
#include <GL/glext.h>

namespace PhysBAM
{
// Note: after each invokation of a function, camera stays consistent
template<class T>
class OPTIX_MOUSE_CAMERA_MANIPULATOR:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
    OPTIX_CAMERA<T> *camera;
    int button,state,modifiers;
    int old_x,old_y;
public:
    OPTIX_MOUSE_CAMERA_MANIPULATOR(OPTIX_CAMERA<T> *_camera):camera(_camera),old_x(-1),old_y(-1){}

    void setState(int _button,int _state,int _modifiers) 
    {
        if ((_button==GLUT_RIGHT_BUTTON && _state==GLUT_UP)||(_button==GLUT_LEFT_BUTTON && _state==GLUT_UP)) old_x=old_y=-1;
        button=_button;state=_state;modifiers=_modifiers;
    }
    void setCamera(OPTIX_CAMERA<T>* _camera){camera=_camera;}

    void handleMouseMove(int x,int y) 
    {
        if (old_x==-1&&old_y==-1){old_x=x;old_y=y;return;}
        T relative_x=(x-old_x)/(T)(camera->screen_width);
        T relative_y=(y-old_y)/(T)(camera->screen_height);

        if(modifiers==0&&button==GLUT_LEFT_BUTTON){
            camera->rotateAroundRight(relative_y);
            camera->rotateY(relative_x);
        } 
        else if(modifiers==0&&button==GLUT_RIGHT_BUTTON){
            camera->moveStraight(relative_y * 10);
        }
        old_x=x;old_y=y;
    }
};
}
#endif
#endif
