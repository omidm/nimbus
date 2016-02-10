//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D__
#define __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D__

#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D.h>
namespace PhysBAM{

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T> class RED_GREEN_GRID_2D;
#endif

template<class T,class RW=T>
class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D:public OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D<T,RW>
{
    typedef VECTOR<T,2> TV;
    typedef OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D<T,RW> BASE;
    using BASE::draw;using BASE::is_animation;using BASE::frame;
    using BASE::frame_loaded;using BASE::prefix;using BASE::display_mode;using BASE::valid;
public:
    using BASE::deformable_geometry_collection;
    using BASE::triangulated_area_objects;
    using BASE::triangles_of_material_objects;
    using BASE::segmented_curve_objects;
    using BASE::velocity_scale;
    using BASE::draw_velocities;
    using BASE::color_map;
    using BASE::free_particles_objects;
    using BASE::free_particles_indirect_arrays;

    COLLISION_GEOMETRY_COLLECTION<TV> collision_body_list;
    DEFORMABLE_BODY_COLLECTION<TV> deformable_body_collection;
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*> embedded_curve_objects;
    bool has_embedded_objects;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    ARRAY<RED_GREEN_GRID_2D<T>*> grid_list;
#endif
    ARRAY<ARRAY<T>*> phi_list;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    ARRAY<OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<T> >*> phi_objects;
#endif

    OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D(const std::string& prefix,const int start_frame);
    ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D();
    
    virtual void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION* Get_Selection(GLuint* buffer, int buffer_size);
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    void Set_Vector_Size(const T vector_size);

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Draw_Velocities();

    // TODO: Is overiding callbacks problematic?
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D,Decrease_Vector_Size,"Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D,Toggle_Draw_Velocities,"Toggle draw velocities");

    virtual void Reinitialize(bool force=false) PHYSBAM_OVERRIDE;    // Needs to be called after some state changes
};
}
#endif
