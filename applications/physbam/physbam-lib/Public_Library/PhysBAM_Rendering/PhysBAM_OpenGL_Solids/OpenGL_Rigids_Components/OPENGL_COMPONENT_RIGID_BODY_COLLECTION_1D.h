//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class T,class RW=T>
class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D : public OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>
{
protected:
    typedef VECTOR<T,1> TV;
    typedef OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW> BASE;

    using BASE::basedir;using BASE::frame_loaded;using BASE::valid;using BASE::show_object_names;using BASE::output_positions;
    using BASE::draw_velocity_vectors;using BASE::draw_node_velocity_vectors;using BASE::draw_point_simplices;
    using BASE::needs_init;using BASE::needs_destroy;using BASE::has_init_destroy_information;
public:
    using BASE::rigid_geometry_collection;using BASE::is_animation;using BASE::draw;using BASE::frame;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
protected:
    using BASE::opengl_point_simplices;using BASE::opengl_axes;using BASE::draw_object;using BASE::use_object_bounding_box;
    using BASE::positions;using BASE::node_positions;using BASE::current_selection;
private:
    bool need_destroy_rigid_body_collection;

public:
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(const std::string& basedir);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir);
    ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D();
    
//#####################################################################
    virtual void Reinitialize(const bool force=false,const bool read_geometry=true) PHYSBAM_OVERRIDE; // Needs to be called after some state changes
protected:
    virtual void Update_Object_Labels();
//#####################################################################
};

}

#endif
