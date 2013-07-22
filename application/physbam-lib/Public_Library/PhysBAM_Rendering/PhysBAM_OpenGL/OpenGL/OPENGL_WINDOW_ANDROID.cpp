//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_OPENGLES
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW_ANDROID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <android/input.h>
#include <android/log.h>
#include <jni.h>

namespace PhysBAM
{

OPENGL_WINDOW_ANDROID* OPENGL_WINDOW_ANDROID::single=0;
//#####################################################################
// OPENGL_WINDOW_ANDROID
//#####################################################################
OPENGL_WINDOW_ANDROID::
OPENGL_WINDOW_ANDROID(OPENGL_WORLD& opengl_world_input,const std::string& window_title_input,const int width_input,const int height_input)
    :OPENGL_WINDOW(opengl_world_input),width(width_input),height(height_input)
{
    if(single) PHYSBAM_FATAL_ERROR("Only one android context allowed");
    single=this;

    /*static int argc=1;static const char *(argv[1]);argv[0]="Visualization";
    glutInit(&argc,(char**)argv);
    glutInitDisplayMode(ANDROID_DOUBLE|ANDROID_RGBA|ANDROID_DEPTH|ANDROID_ALPHA);
    glutInitWindowSize(width,height);
    main_window=glutCreateWindow(window_title_input.c_str());
    glutDisplayFunc(Handle_Display_Glut);
    glutReshapeFunc(Handle_Reshape_Glut);
    glutKeyboardFunc(Handle_Keypress_Glut);
    glutSpecialFunc(Handle_Special_Keypress_Glut);
    glutMouseFunc(Handle_Click_Glut);
    glutMotionFunc(Handle_Drag_Glut);*/
}
//#####################################################################
// ~OPENGL_WINDOW_ANDROID
//#####################################################################
OPENGL_WINDOW_ANDROID::
~OPENGL_WINDOW_ANDROID()
{single=0;}

//#####################################################################
// Setup_Idle
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Setup_Idle(const bool use)
{
  //glutIdleFunc(use?Handle_Idle_Glut:0);
}
//#####################################################################
// Setup_Timer
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Setup_Timer(const float wait_milliseconds)
{
  //glutTimerFunc((int)(wait_milliseconds*1000)+1,Handle_Timer_Glut,0);
}
//#####################################################################
// Redisplay
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Redisplay()
{
  //glutPostRedisplay();
}
//#####################################################################
// Main_Loop
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Main_Loop()
{
  //glutMainLoop();
}
//#####################################################################
// Request_Resize
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Request_Resize(const int width,const int height)
{
  //glutReshapeWindow(width,height);
}
//#####################################################################
// Request_Move
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Request_Move(const int x,const int y)
{
  //glutPositionWindow(x,y);
}
//#####################################################################
// Width
//#####################################################################
int OPENGL_WINDOW_ANDROID::
Width() const
{
    return width;
}
//#####################################################################
// Height
//#####################################################################
int OPENGL_WINDOW_ANDROID::
Height() const
{
    return height;
}
//#####################################################################
// Handle_Display_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Display_Android()
{
    single->opengl_world.Render_World(false); // render, no selection
}
//#####################################################################
// Handle_Reshape_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Reshape_Android(int w,int h)
{
    single->width=w;single->height=h;
    single->opengl_world.Handle_Reshape_Main();
}
//#####################################################################
// Handle_Keypress_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Keypress_Android(unsigned char key,int ctrl,int shift)
{
    if(key>110) return;

    int x=0,y=0;
    OPENGL_KEY real_key(key,0);
    if(key>=AKEYCODE_A && key<=AKEYCODE_Z){
        if(shift) real_key.key=key-AKEYCODE_A+65;
        else real_key.key=key-AKEYCODE_A+97;}
    else if(key>=AKEYCODE_0 && key<=AKEYCODE_9) real_key.key=key-AKEYCODE_0+48;
    else if(key==AKEYCODE_BACKSLASH) real_key.key='\\';
    else if(key==AKEYCODE_LEFT_BRACKET) real_key.key='[';
    else if(key==AKEYCODE_RIGHT_BRACKET) real_key.key=']';
    else if(key==AKEYCODE_EQUALS) real_key.key='=';
    else if(key==AKEYCODE_MINUS) real_key.key='-';
  
    if(ctrl) real_key.modifiers|=OPENGL_KEY::CTRL;

    if(single->opengl_world.getPromptMode())
        single->opengl_world.Handle_Keypress_Prompt(real_key.key);
    else
        single->opengl_world.Handle_Keypress_Main(real_key,x,y);
}
//#####################################################################
// Handle_Click_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Click_Android(int button,int state,int x,int y)
{
    bool ctrl_pressed=false;
    bool shift_pressed=false;
    single->opengl_world.Handle_Click_Main(button,state,x,y,ctrl_pressed,shift_pressed);
}
//#####################################################################
// Handle_Pinch_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Pinch_Android(int button,int state,int x,int y)
{
    bool ctrl_pressed=false;
    bool shift_pressed=false;
    single->opengl_world.Handle_Click_Main(button,state,x,y,ctrl_pressed,shift_pressed);
}
//#####################################################################
// Handle_Drag_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Drag_Android(int x,int y)
{
    single->opengl_world.Handle_Drag_Main(x,y);
}
//#####################################################################
// Get_Camera
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Get_Camera(VECTOR<float,3> &camera,VECTOR<float,3> &target,VECTOR<float,3> &up)
{
    OPENGL_WORLD& real_opengl_world = dynamic_cast<OPENGL_WORLD&>(single->opengl_world);
    real_opengl_world.Get_Look_At(camera,target,up);
}
//#####################################################################
// Set_Camera
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Set_Camera(VECTOR<float,3> camera,VECTOR<float,3> target,VECTOR<float,3> up)
{
    OPENGL_WORLD& real_opengl_world = dynamic_cast<OPENGL_WORLD&>(single->opengl_world);
    real_opengl_world.Set_Look_At(camera,target,up);
}
//#####################################################################
// Handle_Idle_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Idle_Android()
{
    single->opengl_world.Handle_Idle();
}
//#####################################################################
// Handle_Timer_Android
//#####################################################################
void OPENGL_WINDOW_ANDROID::
Handle_Timer_Android(int value)
{
    single->opengl_world.Handle_Timer();
}
//#####################################################################
}
#endif
