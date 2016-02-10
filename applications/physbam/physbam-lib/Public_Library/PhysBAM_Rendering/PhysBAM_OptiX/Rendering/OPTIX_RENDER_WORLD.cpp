//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDER_WORLD
//#####################################################################
#ifdef USE_OPTIX
#include <PhysBAM_Tools/Images/JPG_FILE.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/GLUT_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW_GLUT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_GEOMETRY.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Optix/Rendering_Lights/OPTIX_RENDERING_SPOT_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Lights/OPTIX_RENDERING_FIRE_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Shaders/OPTIX_COMMONSTRUCTS.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include <string>
using namespace PhysBAM;
using namespace optix;

template <class T> OPTIX_RENDER_WORLD<T>* OPTIX_RENDER_WORLD<T>::global_render_world=NULL;

////initialize OpenGL extended APIs
PFNGLGENBUFFERSARBPROC glGenBuffers=NULL;
PFNGLBINDBUFFERPROC glBindBuffer=NULL;
PFNGLBUFFERDATAPROC glBufferData=NULL;
void Initialize_OpenGL_Extended_APIs()
{
    static bool init=false;
    if(!init){
        glGenBuffers=(PFNGLGENBUFFERSARBPROC)wglGetProcAddress("glGenBuffers");
        if(!glGenBuffers)PHYSBAM_FATAL_ERROR("glGenBuffers not supported");
        glBindBuffer=(PFNGLBINDBUFFERPROC)wglGetProcAddress("glBindBuffer");
        if(!glBindBuffer)PHYSBAM_FATAL_ERROR("glBindBuffer not supported");
        glBufferData=(PFNGLBUFFERDATAPROC)wglGetProcAddress("glBufferData");
        if(!glBufferData)PHYSBAM_FATAL_ERROR("glBufferData not supported");
        init=true;
    }
}
#ifdef WIN32
std::string optix_ptx_file_path="\\Public_Library\\PhysBAM_Rendering\\PhysBAM_OptiX\\Rendering_Shaders\\Compiled\\cuda_compile_ptx";
#else
std::string optix_ptx_file_path="/Public_Library/PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Shaders/Compiled/cuda_compile_ptx";
#endif
//#####################################################################
// Constructor
//#####################################################################
template<class T> 
OPTIX_RENDER_WORLD<T>::OPTIX_RENDER_WORLD()
:camera(new OPTIX_CAMERA<T>()),mouse(new OPTIX_MOUSE_CAMERA_MANIPULATOR<T>(camera)),render_shadows(false),shader_prefix(getenv("PHYSBAM")+optix_ptx_file_path),
entry_point_index(0u),RAY_TYPE_RADIANCE(0u),RAY_TYPE_SHADOW(1u),ray_depth_limit(6),scene_epsilon((T)1e-4),background_color(0,0,0),ambient_light_color((T)0.1,(T)0.1,(T)0.1),
use_vbo(true),use_debug(true),output_width(camera->screen_width),output_height(camera->screen_height),output_buffer_format(RT_FORMAT_UNSIGNED_BYTE4),output_vbo_index(-1),output_texture_index(-1),
timer_index_fps(1),timer_index_rendering_fps(2),timer_index_handle_timer(3),timer_index_handle_idle(4)
{
}

template<class T>
void PhysBAM::OPTIX_RENDER_WORLD<T>::Set_Camera( OPTIX_CAMERA<T>* camera_input )
{
    if(camera)delete camera;
    camera=camera_input;
    mouse->setCamera(camera);
    output_width=camera->screen_width;
    output_height=camera->screen_height;
}

template<class T>
void PhysBAM::OPTIX_RENDER_WORLD<T>::Initialize_Optix_Context()
{
    try{
        rt_context=Context::create();
        rt_context->setRayTypeCount(2);
        rt_context["RAY_TYPE_R"]->setUint(RAY_TYPE_RADIANCE);
        rt_context["RAY_TYPE_S"]->setUint(RAY_TYPE_SHADOW);

        // one entry point for phong attenutation and another for photon emission
        ////only one entry point used here...
        rt_context->setEntryPointCount(1);
        rt_context->setStackSize(1024);
    }
    catch(Exception e){PHYSBAM_FATAL_ERROR(e.getErrorString());}

    ////default parameter values
    Set_Ray_Depth(ray_depth_limit);
    Set_Scene_Epsilon(scene_epsilon);
    Set_Background_Color(background_color);
    Set_Ambient_Light_Color(ambient_light_color);
}

template<class T>
void PhysBAM::OPTIX_RENDER_WORLD<T>::Initialize_Optix_Programs()
{
    // add ray generation, exception, and miss programs
    std::string path_to_ptx=shader_prefix+"_generated_OPTIX_RAY_GEN.cu.ptx";
    try{
        rt_ray_gen_program=rt_context->createProgramFromPTXFile(path_to_ptx,"ray_gen");
        rt_context->setRayGenerationProgram(entry_point_index,rt_ray_gen_program);
        rt_exception_program=rt_context->createProgramFromPTXFile(path_to_ptx,"exception");
        rt_context->setExceptionProgram(entry_point_index,rt_exception_program);
        path_to_ptx=shader_prefix+"_generated_OPTIX_MISS.cu.ptx";
        rt_miss_program=rt_context->createProgramFromPTXFile(path_to_ptx,"miss");
        rt_context->setMissProgram(entry_point_index,rt_miss_program);
    }
    catch(Exception e){PHYSBAM_FATAL_ERROR(e.getErrorString());}
}

template<class T>
void PhysBAM::OPTIX_RENDER_WORLD<T>::Create_Output_Buffer()
{
    Buffer rt_result_color_buffer=use_vbo?Create_Output_Buffer_VBO():Create_Output_Buffer_MEM();
    rt_context["output_buffer"]->set(rt_result_color_buffer);
    if(use_vbo)Initialize_Output_Texture();
}

template<class T> 
void OPTIX_RENDER_WORLD<T>::Update_Optix_Camera() 
{
    // pass camera parameters to ray generation program
    // pass projection parameters
    rt_context["halfed_proj_metrics"]->setFloat(camera->proj_width*0.5f,camera->proj_height*0.5f);
    rt_context["screen_metrics"]->setUint(camera->screen_width,camera->screen_height);

    // pass camera orientation parameteres
    // TODO: expand this piece as another function as it will be invoked every frame
    TV camera_from_loc_to_focal=(camera->focal_point-camera->position).Normalized();

    rt_context["up"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(camera->vertical_vector));
    rt_context["right"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(camera->horizontal_vector));
    rt_context["dir"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(camera->look_vector));
    rt_context["loc"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(camera->position));
    rt_context["from_loc_to_focal"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(camera_from_loc_to_focal*camera->proj_distance));
}

template<class T>
Buffer PhysBAM::OPTIX_RENDER_WORLD<T>::Create_Output_Buffer_VBO( RTformat format )
{
    Initialize_OpenGL_Extended_APIs();
    output_buffer_format=format;
    glGenBuffers(1,&output_vbo_index);
    glBindBuffer(GL_ARRAY_BUFFER,output_vbo_index);
    size_t element_size;
    rtuGetSizeForRTformat(format,&element_size);
    glBufferData(GL_ARRAY_BUFFER,element_size*output_width*output_height,0,GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    Buffer buffer=rt_context->createBufferFromGLBO(RT_BUFFER_OUTPUT,output_vbo_index);
    buffer->setFormat(output_buffer_format);buffer->setSize(output_width,output_height);
    return buffer;
}

template<class T>
Buffer PhysBAM::OPTIX_RENDER_WORLD<T>::Create_Output_Buffer_MEM( RTformat format )
{
    output_buffer_format=format;
    return rt_context->createBuffer(RT_BUFFER_OUTPUT,output_buffer_format,output_width,output_height);
}

//#####################################################################
// Launch - Getting OPTIX_RENDER_WORLD launching context for a 0-th entry point
//#####################################################################
template<class T> void OPTIX_RENDER_WORLD<T>::Launch() 
{
    rt_context->launch(entry_point_index,output_width,output_height);
}

//#####################################################################
// Instance - Getting OPTIX_RENDER_WORLD instance (it's a singleton class)
//#####################################################################
template<class T> OPTIX_RENDER_WORLD<T>* OPTIX_RENDER_WORLD<T>::Instance()
{
  if(global_render_world==NULL)global_render_world=new OPTIX_RENDER_WORLD<T>();return global_render_world;
}

//#####################################################################
// Destructor
//#####################################################################
template<class T> OPTIX_RENDER_WORLD<T>::~OPTIX_RENDER_WORLD() 
{
    rt_context->destroy();
    Clean_Scene();
}

// creates rt light buffers from 'lights' variable
// TODO: separate initialization from lights update (or do we need moving or dissapearing lights at all?)
template<class T> void OPTIX_RENDER_WORLD<T>::Create_Lights_Buffer() 
{
    light_buffer=rt_context->createBuffer(RT_BUFFER_INPUT,RT_FORMAT_USER);
    light_buffer->setElementSize(sizeof(BasicLight));
    light_buffer->setSize(optix_lights.m);
    rt_context["lights"]->set(light_buffer);
}

template<class T> void PhysBAM::OPTIX_RENDER_WORLD<T>::Add_Spot_Light( OPTIX_RENDERING_SPOT_LIGHT<T>* spot_light )
{
    spot_light->light_index=spot_lights.m+1;
    spot_light->start_index=optix_lights.m+1;
    spot_lights.Append(spot_light);
    optix_lights.Append(*(spot_light->Get_Basic_Light()));
}

template<class T> void PhysBAM::OPTIX_RENDER_WORLD<T>::Add_Fire_Light( OPTIX_RENDERING_FIRE_LIGHT<T>* fire_light )
{
    fire_lights.Append(fire_light);
    fire_light->start_index=optix_lights.m+1;
    for(int i=1;i<=fire_light->spot_light_number;i++){
        optix_lights.Append(*(fire_light->Get_Basic_Light(i)));
    }
}

template<class T> void OPTIX_RENDER_WORLD<T>::Create_Geometry() 
{
    // place all geometry under top level group
    Group top_level_group = rt_context->createGroup();

    ARRAY<Transform> childTransforms;
    for (int i=1;i<=objects.m;i++){
        objects(i)->Initialize(this);
        childTransforms.Append_Elements(objects(i)->getTransforms());
    }
    top_level_group->setChildCount(childTransforms.m);

    for (int i=1;i<=childTransforms.m;i++){
        top_level_group->setChild(i-1,childTransforms(i));
    }
    rt_context["top_object"]->set(top_level_group);

    top_level_acceleration=rt_context->createAcceleration("NoAccel","NoAccel");
    top_level_acceleration->markDirty();
    top_level_group->setAcceleration(top_level_acceleration);

    Group top_level_opaque_group=rt_context->createGroup();
    ARRAY<Transform> childOpaqueTransforms;
    for(int i=1;i<=objects.m;i++){
        if(objects(i)->getMaterial()->opaque)childOpaqueTransforms.Append_Elements(objects(i)->getTransforms());
    }
    top_level_opaque_group->setChildCount(childOpaqueTransforms.m);

    for(int i=1;i<=childOpaqueTransforms.m;i++){
        top_level_opaque_group->setChild(i-1,childOpaqueTransforms(i));
    }
    rt_context["top_opaque_object"]->set(top_level_opaque_group);

    top_level_opaque_acceleration=rt_context->createAcceleration("NoAccel","NoAccel");
    top_level_opaque_group->setAcceleration(top_level_opaque_acceleration);
    top_level_opaque_acceleration->markDirty();
}

template<class T> void OPTIX_RENDER_WORLD<T>::Render_World(bool selecting, bool swap_buffers) 
{
    timer.Reset_Time(timer_index_rendering_fps);
    Update_Objects(0);
    Update_Lights(0);
    Update_Optix_Camera();
    for(int i=1;i<=objects.m;i++){
        objects(i)->PrepareForDraw();
    }
    Launch();
    if(use_vbo)Draw_Output_Buffer_VBO();
    else Draw_Output_Buffer_MEM();
    
    if(display_fps){
        double elapsed_t=timer.Peek_And_Reset_Time(timer_index_fps);
        double rendering_elapsed_t=timer.Peek_And_Reset_Time(timer_index_rendering_fps);
        char text[64];
        sprintf_s(text,"fps: %f, rendering fps: %f",(float)(1000.0/elapsed_t),(float)(1000.0/rendering_elapsed_t));
        OpenGL_String(VECTOR<T,2>(0,(T)0.98),text);
    }
}

template<class T>
void PhysBAM::OPTIX_RENDER_WORLD<T>::Draw_Output_Buffer_MEM()
{
    GLvoid* imageData;
    GLsizei width,height;
    RTsize buffer_width,buffer_height;
    RTformat buffer_format;

    Buffer buffer=Get_Output_Buffer();

    imageData=buffer->map();
    buffer_format=buffer->getFormat();
    PHYSBAM_ASSERT(buffer_format==RT_FORMAT_UNSIGNED_BYTE4);

    if (0==imageData){
        std::cerr<<"Data in rt_result_color_buffer is null\n";
        exit(2);
    }
    buffer->getSize(buffer_width,buffer_height);
    width=static_cast<GLsizei>(buffer_width);
    height=static_cast<GLsizei>(buffer_height);

    GLenum gl_data_type=GL_UNSIGNED_BYTE;
    GLenum gl_format=GL_BGRA;

    if(write_images&&!isPaused) {
        static unsigned fileNumber=0;
        std::string filename=image_path;
        char* fileNumStr=(char*)malloc(std::strlen("render-000.jpg")*sizeof(char));
        std::sprintf(fileNumStr,"/render-%03d.jpg",fileNumber);
        filename.append(fileNumStr);
        unsigned char *bufferPointer=(unsigned char*)imageData;
        ARRAY<VECTOR<T,4>,VECTOR<int,2> > outputImage(1,height,1,width);
        for (int i=0;i<height;i++) {
            for (int j=0;j<width;j++) {
                unsigned t=(i*width+j)*4;
                //reading BGRA, writing RGBA
                outputImage(j+1,i+1)=VECTOR<T,4>(IMAGE<T>::Byte_Color_To_Scalar_Color(bufferPointer[t+2]),
                    IMAGE<T>::Byte_Color_To_Scalar_Color(bufferPointer[t+1]),
                    IMAGE<T>::Byte_Color_To_Scalar_Color(bufferPointer[t+0]),
                    IMAGE<T>::Byte_Color_To_Scalar_Color(bufferPointer[t+3]));
            }
        }
        JPG_FILE<T>::Write(filename,outputImage);
        fileNumber++;
    }

    glDrawPixels(width,height,gl_format,gl_data_type,imageData);
    // Now unmap the buffer
    buffer->unmap();
}

template<class T>
void PhysBAM::OPTIX_RENDER_WORLD<T>::Initialize_Output_Texture()
{
    glGenTextures(1,&output_texture_index);
    glBindTexture(GL_TEXTURE_2D,output_texture_index);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D,0);
}

template<class T>
void PhysBAM::OPTIX_RENDER_WORLD<T>::Draw_Output_Buffer_VBO()
{
    glBindTexture(GL_TEXTURE_2D,output_texture_index);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER,output_vbo_index);
    RTsize elementSize=Get_Output_Buffer()->getElementSize();
    if((elementSize%8)==0)glPixelStorei(GL_UNPACK_ALIGNMENT,8);
    else if((elementSize%4)==0)glPixelStorei(GL_UNPACK_ALIGNMENT,4);
    else if((elementSize%2)==0)glPixelStorei(GL_UNPACK_ALIGNMENT,2);
    else glPixelStorei(GL_UNPACK_ALIGNMENT,1);

    if(output_buffer_format==RT_FORMAT_UNSIGNED_BYTE4)glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,output_width,output_height,0,GL_BGRA,GL_UNSIGNED_BYTE,0);
    else if(output_buffer_format==RT_FORMAT_FLOAT4)glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F_ARB,output_width,output_height,0,GL_RGBA,GL_FLOAT,0);
    else if(output_buffer_format==RT_FORMAT_FLOAT3)glTexImage2D(GL_TEXTURE_2D,0,GL_RGB32F_ARB,output_width,output_height,0,GL_RGB,GL_FLOAT,0);
    else if(output_buffer_format==RT_FORMAT_FLOAT)glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE32F_ARB,output_width,output_height,0,GL_LUMINANCE,GL_FLOAT,0);
    else PHYSBAM_FATAL_ERROR("Unknown buffer format");

    glBindBuffer(GL_PIXEL_UNPACK_BUFFER,0);
    float u=0.5f/output_width;  
    float v=0.5f/output_height; ////pixel offsets
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glBegin(GL_QUADS);
    glTexCoord2f(u,v);
    glVertex2f(0,0);
    glTexCoord2f(1.0f,v);
    glVertex2f(1.0f,0);
    glTexCoord2f(1.0f-u,1.0f-v);
    glVertex2f(1.0f,1.0f);
    glTexCoord2f(u,1.0f-v);
    glVertex2f(0,1.0f);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

template<class T> void OPTIX_RENDER_WORLD<T>::Handle_Keypress_Main(const OPENGL_KEY &key,int x,int y) 
{
    for(int i=1;i<=global_render_world->key_listeners.m;i++){
        global_render_world->key_listeners(i)->KeyPressed(key.key,x,y);
    }
}

template<class T> void OPTIX_RENDER_WORLD<T>::Handle_Drag_Main(int x,int y) 
{
    for (int i=1;i<=global_render_world->mouse_listeners.m;i++) {
        global_render_world->mouse_listeners(i)->Handle_Drag(x,y);
    }
    global_render_world->mouse->handleMouseMove(x,y);
    window->Redisplay();
}

template<class T> void OPTIX_RENDER_WORLD<T>::Handle_Click_Main(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed)
{
    for(int i=1;i<=global_render_world->mouse_listeners.m;i++){
        global_render_world->mouse_listeners(i)->Handle_Click(button,state,x,y,ctrl_pressed,shift_pressed);
    }
    switch(state){
    case GLUT_DOWN:
        global_render_world->mouse->setState(button,state,glutGetModifiers());
        global_render_world->mouse->handleMouseMove(x,y);
        window->Redisplay();
        break;
    case GLUT_UP:
        global_render_world->mouse->setState(button,state,glutGetModifiers());
        break;
    }
}

template<class T> void OPTIX_RENDER_WORLD<T>::Handle_Timer() 
{
    //// Updating
    double t=timer.Peek_And_Reset_Time(timer_index_handle_timer);

    // time at first frame is equal to zero, if decide to change check nobidy rely on that
    global_render_world->Update_Objects((T)t);
    window->Redisplay();
}

template<class T> void OPTIX_RENDER_WORLD<T>::Handle_Idle() 
{
    //// Updating
    double t=timer.Peek_And_Reset_Time(timer_index_handle_idle);

    //// time at first frame is equal to zero, if decide to change check nobidy rely on that
    global_render_world->Update_Objects((T)t);
    window->Redisplay();
}

template<class T> void OPTIX_RENDER_WORLD<T>::Update_Objects(T time) 
{
    for(int i=1;i<=objects.m;i++){
        objects(i)->Update(time);
    }
    // TODO!!!! Check that anything changed to geometry
    top_level_acceleration->markDirty();
    top_level_opaque_acceleration->markDirty();
}

template<class T> void OPTIX_RENDER_WORLD<T>::Update_Lights(T time) 
{
    ////update spot lights
    for(int i=1;i<=spot_lights.m;i++){
        int start=spot_lights(i)->start_index;
        memcpy(&optix_lights(start),spot_lights(i)->Get_Basic_Light(),sizeof(BasicLight));
    }
    ////update fire lights
    for(int i=1;i<=fire_lights.m;i++){
        int start=fire_lights(i)->start_index;
        memcpy(&optix_lights(start),fire_lights(i)->Get_Basic_Light(1),fire_lights(i)->spot_light_number*sizeof(BasicLight));
    }
    ////send updated lights to optix
    if(optix_lights.m>0)memcpy(light_buffer->map(),&optix_lights(1),sizeof(BasicLight)*optix_lights.m);light_buffer->unmap();
}

template<class T> void OPTIX_RENDER_WORLD<T>::Initialize_OpenGL_Perspective() 
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,1,0,1,-1,1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0,0,output_width,output_height);
}

////Initialize() should be called after all objects and materials are added
template<class T> void OPTIX_RENDER_WORLD<T>::Initialize() 
{
    window->Setup_Idle(true);    
    Initialize_OpenGL_Perspective();
    Initialize_Optix_Context();
    Initialize_Optix_Programs();
    Update_Optix_Camera();

    Create_Lights_Buffer();
    Create_Output_Buffer();
    Create_Geometry();

    for(int i=1;i<=objects.m;i++)objects(i)->BeforeLaunch();
    rt_context->setPrintEnabled(use_debug);

    try{rt_context->validate();rt_context->compile();}
    catch(Exception e){PHYSBAM_FATAL_ERROR(e.getErrorString());}
}

template class OPTIX_RENDER_WORLD<float>;
// template class OPTIX_RENDER_WORLD<double>;
#endif
