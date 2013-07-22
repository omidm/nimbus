//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Tamar Shinar, Rachel Weinstein, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D.h>
namespace PhysBAM{

template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGID_BODY_COLLECTION;
template<class T,class RW=T>
class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D : public OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>
{
protected:
    typedef VECTOR<T,2> TV;
    typedef OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW> BASE;

    using BASE::basedir;
    using BASE::frame_loaded;
    using BASE::valid;
    using BASE::show_object_names;
    using BASE::output_positions;
    using BASE::draw_velocity_vectors;
    using BASE::draw_individual_axes;
    using BASE::draw_node_velocity_vectors;
    using BASE::draw_segmented_curve;
    using BASE::draw_triangulated_area;
    using BASE::draw_implicit_curve;
    using BASE::needs_init;
    using BASE::needs_destroy;
    using BASE::has_init_destroy_information;

    bool draw_articulation_points;
    bool draw_forces_and_torques;
    bool draw_linear_muscles;
public:
    using BASE::frame;using BASE::draw;using BASE::is_animation;using BASE::slice;
    using BASE::rigid_geometry_collection;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body;
protected:
    using BASE::opengl_segmented_curve;
    using BASE::opengl_triangulated_area;
    using BASE::opengl_levelset;
    using BASE::extra_components;
    using BASE::opengl_axes;
    using BASE::draw_object;
    using BASE::use_object_bounding_box;
    using BASE::positions;
    using BASE::velocity_vectors;
    using BASE::node_positions;
    using BASE::node_velocity_vectors;
    using BASE::velocity_field;
    using BASE::node_velocity_field;
    using BASE::current_selection;

    bool need_destroy_rigid_body_collection;
    ARRAY<VECTOR<T,2> > articulation_points;
    ARRAY<PAIR<VECTOR<T,2>,T>,int> forces_and_torques;

public:
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(const std::string& basedir);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir);
    ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D();
    
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;

    void Read_Articulated_Information(const std::string& filename);

    void Toggle_Articulation_Points();
    void Toggle_Linear_Muscles();
    void Toggle_Forces_And_Torques();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D, Toggle_Articulation_Points, "Toggle articulation points");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D, Toggle_Linear_Muscles, "Toggle linear muscles");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D, Toggle_Forces_And_Torques, "Toggle forces and torques");

public:
    virtual void Reinitialize(const bool force=false,const bool read_geoemtry=true);    // Needs to be called after some state changes
protected:
    using BASE::Create_Geometry;
    using BASE::Update_Geometry;
    using BASE::Destroy_Geometry;
    virtual void Update_Object_Labels();
};

template<class T>
class OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D : public OPENGL_SELECTION
{
public:
    int joint_id;

    OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D(OPENGL_OBJECT *object,const int joint_id)
        :OPENGL_SELECTION(OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_2D,object),joint_id(joint_id)
    {}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D : public OPENGL_SELECTION
{
public:
    int muscle_id;
    OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D(OPENGL_OBJECT *object,const int muscle_id)
        :OPENGL_SELECTION(OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_MUSCLE_2D,object),muscle_id(muscle_id)
    {}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
