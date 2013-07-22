//#####################################################################
// Copyright 2006, Kevin Der.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_MUSCLE_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_MUSCLE_3D__
#define __OPENGL_COMPONENT_MUSCLE_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T> class OPENGL_TRIANGULATED_SURFACE;
template<class TV> class ARTICULATED_RIGID_BODY;

template<class T>
class OPENGL_COMPONENT_MUSCLE_3D:public OPENGL_COMPONENT
{
public:
    static const int length_resolution=60;
    static const int radial_resolution=20;
    typedef VECTOR<T,3> TV;
    bool draw_linear_muscles;
    bool draw_surface_muscles;
    bool draw_muscle_internal_particles;
    OPENGL_COLOR_MAP<T>* muscle_color_map;
    OPENGL_SELECTION* current_selection;
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body;

    ARRAY<OPENGL_TRIANGULATED_SURFACE<T> *> opengl_triangulated_surface;
    ARRAY<PAIR<int,int> > surface_muscle_indices;
    ARRAY<TV> muscle_internal_particles;

    OPENGL_COMPONENT_MUSCLE_3D(OPENGL_COLOR_MAP<T>* muscle_color_map_input);
    virtual ~OPENGL_COMPONENT_MUSCLE_3D();

    void Toggle_Linear_Muscles() {draw_linear_muscles=!draw_linear_muscles;}
    void Toggle_Surface_Muscles() {draw_surface_muscles=!draw_surface_muscles;}
    void Toggle_Muscle_Internal_Particles() {draw_muscle_internal_particles=!draw_muscle_internal_particles;}
    
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MUSCLE_3D, Toggle_Linear_Muscles, "Toggle linear muscles");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MUSCLE_3D, Toggle_Surface_Muscles, "Toggle surface muscles");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MUSCLE_3D, Toggle_Muscle_Internal_Particles, "Toggle internal muscle particles");

//#####################################################################
    virtual void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION* Get_Muscle_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE; // not called by OPENGL_COMPONENT_RIGID_BODIES_3D
    void Clear_Highlight() PHYSBAM_OVERRIDE; // not called by OPENGL_COMPONENT_RIGID_BODIES_3D
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    void Initialize(ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body_input,std::string basedir,int frame);
    void Read_Muscle_Internal_Particles(const std::string& filename);
//#####################################################################
};

template<class T>
class OPENGL_SELECTION_MUSCLE_3D:public OPENGL_SELECTION
{
public:
    int muscle_id,segment_id;
    OPENGL_SELECTION_MUSCLE_3D(OPENGL_OBJECT *object,const int muscle_id_input,const int segment_id_input)
        :OPENGL_SELECTION(OPENGL_SELECTION::MUSCLE_3D,object),muscle_id(muscle_id_input),segment_id(segment_id_input)
    {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_MUSCLE_SURFACE_3D:public OPENGL_SELECTION
{
public:
    int muscle_id,segment_id,surface_index;
    OPENGL_SELECTION_MUSCLE_SURFACE_3D(OPENGL_OBJECT *object,const int muscle_id_input,const int segment_id_input,const int surface_index_input)
        :OPENGL_SELECTION(OPENGL_SELECTION::MUSCLE_SURFACE_3D,object),muscle_id(muscle_id_input),segment_id(segment_id_input),surface_index(surface_index_input)
    {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}
#endif
