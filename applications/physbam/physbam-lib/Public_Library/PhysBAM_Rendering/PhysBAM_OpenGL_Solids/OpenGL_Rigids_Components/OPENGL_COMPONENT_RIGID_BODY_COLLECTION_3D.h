//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_IMPLICIT_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D.h>
namespace PhysBAM{

template<class T> class OPENGL_COMPONENT_MUSCLE_3D;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGID_BODY_COLLECTION;

template<class T,class RW=T>
class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D:public OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>
{
private:
    typedef OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW> BASE;
    typedef VECTOR<T,3> TV;

protected:
    using BASE::basedir;using BASE::use_display_lists;using BASE::frame_loaded;using BASE::valid;using BASE::draw_velocity_vectors;
    using BASE::draw_individual_axes;using BASE::draw_triangulated_surface;using BASE::draw_tetrahedralized_volume;
    using BASE::draw_implicit_surface;using BASE::needs_init;using BASE::needs_destroy;using BASE::has_init_destroy_information;
public:
    using BASE::frame;using BASE::draw;using BASE::is_animation;using BASE::slice;using BASE::is_interactive;
    using BASE::rigid_geometry_collection;using BASE::opengl_colors;using BASE::opengl_triangulated_surface;
protected:
    using BASE::opengl_tetrahedralized_volume;using BASE::opengl_levelset;using BASE::opengl_octree_levelset_surface;using BASE::opengl_axes;using BASE::draw_object;
    using BASE::velocity_field;using BASE::angular_velocity_field;using BASE::draw_simplicial_object_particles;

    bool draw_articulation_points;
    int draw_joint_frames;
    bool draw_forces_and_torques;
public:
    RIGID_BODY_COLLECTION<TV> &rigid_body_collection;
    RIGID_BODY_COLLECTION<TV>* rigid_body_collection_simulated;
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body;
    OPENGL_COMPONENT_MUSCLE_3D<T>* opengl_component_muscle_3d;
private:
    ARRAY<ARRAY<OPENGL_COMPONENT*>,int> extra_components;
    ARRAY<VECTOR<T,3> > articulation_points;
    ARRAY<FRAME<TV> > joint_frames;
    ARRAY<PAIR<VECTOR<T,3>,VECTOR<T,3> >,int> forces_and_torques;
    VECTOR<T,3> projected_COM;
    OPENGL_SELECTION* current_selection;
    bool need_destroy_rigid_body_collection;

public:
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(bool use_display_lists=true);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(const std::string& basedir, bool use_display_lists=true);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir, bool use_display_lists=true);
    virtual ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D();

    virtual void Display(const int in_color=1) const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;

    void Read_Articulated_Information(const std::string& filename);
    void Update_Articulation_Points();

    void Toggle_Articulation_Points();
    void Toggle_Joint_Frames();
    void Toggle_Forces_And_Torques();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Articulation_Points, "Toggle articulation points");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Joint_Frames, "Toggle joint frames");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Forces_And_Torques, "Toggle forces and torques");

    void Resize_Structures(const int size);
    virtual void Reinitialize(const bool force=false,const bool read_geometry=true) PHYSBAM_OVERRIDE; // Needs to be called after some state changes
    virtual void Reinitialize_Without_Files(const bool force=false);    // Needs to be called after some state changes
    void Update_Bodies(const bool update_arb_points=true);
    void Initialize_One_Body(const int body_id,const bool force=false);
private:
    void Initialize();
    void Create_Geometry(const int id);
    void Update_Geometry(const int id);
    void Destroy_Geometry(const int id);
};

template<class T>
class OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D : public OPENGL_SELECTION
{
public:
    int joint_id;

    OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D(OPENGL_OBJECT *object,const int joint_id)
        :OPENGL_SELECTION(OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_3D,object),joint_id(joint_id)
    {}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
