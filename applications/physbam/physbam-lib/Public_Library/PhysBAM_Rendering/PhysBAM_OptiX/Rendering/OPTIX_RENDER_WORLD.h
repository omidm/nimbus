//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDER_WORLD
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDER_WORLD__
#define __OPTIX_RENDER_WORLD__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_CAMERA.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD_KEY_LISTENER.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD_MOUSE_LISTENER.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_MOUSE_CAMERA_MANIPULATOR.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Lights/OPTIX_RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Utilities/OPTIX_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_Optix/Rendering_Shaders/OPTIX_COMMONSTRUCTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_KEY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/GLUT_RENDER_WORLD.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>
#include <string>

namespace PhysBAM{
using namespace optix;
class OPENGL_WINDOW;
template<class T> class OPTIX_RENDERING_OBJECT;
template<class T> class OPTIX_RENDERING_SPOT_LIGHT;
template<class T> class OPTIX_RENDERING_FIRE_LIGHT;

template<class T> 
class OPTIX_RENDER_WORLD : public GLUT_RENDER_WORLD 
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_WINDOW* window;
    OPTIX_CAMERA<T>* camera;
    const std::string shader_prefix;
    const uint RAY_TYPE_RADIANCE;
    const uint RAY_TYPE_SHADOW;

    static OPTIX_RENDER_WORLD* Instance();
    ~OPTIX_RENDER_WORLD();

    void Clean_Scene(){spot_lights.Clean_Memory();fire_lights.Clean_Memory();optix_lights.Clean_Memory();objects.Clean_Memory();}
    void Add_Light(OPTIX_RENDERING_LIGHT<T>* light){light->Add_To_Optix_Render_World(this);}
    void Add_Spot_Light(OPTIX_RENDERING_SPOT_LIGHT<T>* spot_light);
    void Add_Fire_Light(OPTIX_RENDERING_FIRE_LIGHT<T>* fire_light);
    void Add_Object(OPTIX_RENDERING_OBJECT<T>* object){objects.Append(object);}
    void Add_Key_Listener(OPTIX_RENDER_WORLD_KEY_LISTENER *listener){key_listeners.Append(listener);}
    void Add_Mouse_Listener(OPTIX_RENDER_WORLD_MOUSE_LISTENER *listener){mouse_listeners.Append(listener);}
    
    void Set_Ray_Depth(int max_depth){rt_context["max_depth"]->setInt(max_depth);}
    void Set_Scene_Epsilon(T epsilon){rt_context["scene_epsilon"]->setFloat(1e-4f);}    
    void Set_Background_Color(const TV& color){rt_context["bg_color"]->setFloat(color.x,color.y,color.z);}
    void Set_Ambient_Light_Color(const TV& ambient){rt_context["ambient_light_color"]->setFloat(ambient.x,ambient.y,ambient.z);}
    void Set_Camera(OPTIX_CAMERA<T>* camera_input);

    const ARRAY<OPTIX_RENDERING_SPOT_LIGHT<T>*>& Spot_Lights(){return spot_lights;}
    const ARRAY<OPTIX_RENDERING_FIRE_LIGHT<T>*>& Fire_Lights(){return fire_lights;}
    const ARRAY<BasicLight>& Optix_Lights(){return optix_lights;}

    const ARRAY<OPTIX_RENDERING_OBJECT<T>*>& Objects(){return objects;}
    Context RTContext(){return rt_context;}
    int Screen_Width(){return output_width;}
    int Screen_Height(){return output_height;}
    TV getWorldDeltaFromScreenDelta(int dx,int dy,VECTOR<T,3> v){return camera->getWorldDeltaFromScreenDelta(dx,dy,v);}

    void Enable_Fps_Display(bool display_fps_input=true){display_fps=display_fps_input;}
    void Enable_Debug(bool use_debug_input=true){use_debug=use_debug_input;}
    void Set_Image_Output(bool write_images_input,std::string image_path_input){write_images=write_images_input;image_path=write_images?image_path_input:"";}
    void Pause(bool setPause=true){isPaused=setPause;}  ////toggle for internal paused state (used for toggling image writing)
    bool Check_If_Paused(){return isPaused;};
    bool Toggle_Image_Writing(){write_images=!write_images;return write_images;}    ////returns true if set on, false if set off

    void Initialize();
    void Render_Frame(){window->Main_Loop();}

private:
    static OPTIX_RENDER_WORLD* global_render_world;

    ARRAY<OPTIX_RENDER_WORLD_KEY_LISTENER*> key_listeners;
    ARRAY<OPTIX_RENDER_WORLD_MOUSE_LISTENER*> mouse_listeners;

    ARRAY<OPTIX_RENDERING_SPOT_LIGHT<T>*> spot_lights;    ////on cpu end
    ARRAY<OPTIX_RENDERING_FIRE_LIGHT<T>*> fire_lights;    ////on cpu end  
    ARRAY<BasicLight> optix_lights;             ////on gpu end
    ARRAY<OPTIX_RENDERING_OBJECT<T>*> objects;  
    OPTIX_MOUSE_CAMERA_MANIPULATOR<T> *mouse;   ////controlling the camera

    bool render_shadows;
    bool display_fps;
    bool write_images;
    bool isPaused;
    bool use_vbo;
    bool use_debug;
    std::string image_path;
    
    uint entry_point_index;
    int ray_depth_limit;
    T scene_epsilon;
    TV background_color;
    TV ambient_light_color;

    Context rt_context;
    Program rt_ray_gen_program;
    Program rt_miss_program;
    Program rt_exception_program;
    
    Acceleration top_level_acceleration;
    Acceleration top_level_opaque_acceleration;

    Buffer light_buffer;

    int output_width,output_height;
    RTformat output_buffer_format;
    GLuint output_vbo_index;
    GLuint output_texture_index;

    TIMER timer;
    int timer_index_fps,timer_index_rendering_fps,timer_index_handle_timer,timer_index_handle_idle;

    //////////////////////////////////////////////////////////////////////////
    OPTIX_RENDER_WORLD();
    
    Buffer Create_Output_Buffer_VBO(RTformat format=RT_FORMAT_UNSIGNED_BYTE4);
    Buffer Create_Output_Buffer_MEM(RTformat format=RT_FORMAT_UNSIGNED_BYTE4);
    Buffer Get_Output_Buffer(){return global_render_world->rt_context["output_buffer"]->getBuffer();}

    void Create_Lights_Buffer();
    void Create_Output_Buffer();    ////should be called only after OpenGL context is initialized
    void Create_Geometry();

    void Initialize_OpenGL_Perspective();
    void Initialize_Optix_Context();
    void Initialize_Optix_Programs();
    void Initialize_Output_Texture();
    
    void Update_Optix_Camera();
    void Update_Objects(T time);
    void Update_Lights(T time);

    void Draw_Output_Buffer_VBO();
    void Draw_Output_Buffer_MEM();
    
    void Launch();

    virtual bool getPromptMode(){return false;}
    virtual void Render_World(bool selecting,bool swap_buffers=true);
    virtual void Handle_Idle();
    virtual void Handle_Keypress_Main(const OPENGL_KEY& key,int x,int y);
    virtual void Handle_Click_Main(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed);
    virtual void Handle_Drag_Main(int x,int y);
    virtual void Handle_Timer();
    virtual void Handle_Reshape_Main(){}
    virtual void Handle_Keypress_Prompt(unsigned char raw_key){}
};
}
#endif
#endif
