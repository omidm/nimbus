//#####################################################################
// Class GLUT_RENDER_WORLD
//#####################################################################
#ifndef __GLUT_RENDER_WORLD__
#define __GLUT_RENDER_WORLD__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_KEY.h>
#include <ctime>

namespace PhysBAM{

class GLUT_RENDER_WORLD {
public:
    virtual ~GLUT_RENDER_WORLD(){}
    // Prompting
    virtual bool getPromptMode() { return false; }

    virtual void Handle_Idle() {}
    virtual void Handle_Timer() {}

    virtual void Render_World(bool selecting,bool swap_buffers=true) {}
    virtual void Handle_Reshape_Main() {}
    virtual void Handle_Keypress_Main(const OPENGL_KEY& key,int x,int y) {}
    virtual void Handle_Keypress_Prompt(unsigned char raw_key) {}
    virtual void Handle_Click_Main(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed) {}
    virtual void Handle_Drag_Main(int x,int y) {}
//#####################################################################
};
}
#endif


