//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D__
#define __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D__

#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLES_COLLISIONS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D:public OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>
{
    typedef VECTOR<T,3> TV;
    typedef OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW> BASE;
    using BASE::draw;using BASE::is_animation;using BASE::frame;
    using BASE::frame_loaded;using BASE::prefix;using BASE::smooth_shading;using BASE::display_mode;using BASE::valid;
protected:
    int display_soft_bound_surface_mode,display_hard_bound_surface_mode,display_forces_mode,interaction_pair_display_mode;
public:
    using BASE::slice;using BASE::is_interactive;
    using BASE::deformable_geometry;using BASE::active_list;
    using BASE::triangulated_surface_objects;using BASE::tetrahedralized_volume_objects;
    COLLISION_GEOMETRY_COLLECTION<TV> collision_body_list;
    DEFORMABLE_BODY_COLLECTION<TV> deformable_body_collection;
    DEFORMABLE_BODY_COLLECTION<TV>* deformable_body_collection_simulated;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> boundary_surface_objects;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> embedded_surface_objects;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> hard_bound_boundary_surface_objects;
    bool has_embedded_objects,has_soft_bindings;
    ARRAY<POINT_FACE_REPULSION_PAIR<TV> > point_triangle_interaction_pairs;
    ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> > edge_edge_interaction_pairs;
    ARRAY<FORCE_DATA<TV> > force_data_list;
    OPENGL_COLOR_RAMP<T>* color_map_forces;

    OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D();
    OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D(const std::string& prefix,const int start_frame);
    ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D();
    
    virtual void Set_Display_Modes_For_Geometry_Collection(bool& display_triangulated_surface_objects,
            bool& display_tetrahedralized_volume_objects,bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects) const PHYSBAM_OVERRIDE;
    void Set_Display_Modes(bool& display_triangulated_surface_objects,bool& display_tetrahedralized_volume_objects,
        bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects,bool& display_boundary_surface_objects,
        bool& display_hard_bound_boundary_surface_objects) const;
    virtual void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE
    {BASE::Turn_Smooth_Shading_On();
    for(int i=1;i<=boundary_surface_objects.m;i++)if(boundary_surface_objects(i))boundary_surface_objects(i)->Turn_Smooth_Shading_On();
    for(int i=1;i<=embedded_surface_objects.m;i++)if(embedded_surface_objects(i))embedded_surface_objects(i)->Turn_Smooth_Shading_On();}
    
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE
    {BASE::Turn_Smooth_Shading_Off();
    for(int i=1;i<=boundary_surface_objects.m;i++)if(boundary_surface_objects(i))boundary_surface_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=1;i<=embedded_surface_objects.m;i++)if(embedded_surface_objects(i))embedded_surface_objects(i)->Turn_Smooth_Shading_Off();}

    void Cycle_Forces_Mode();
    void Cycle_Hard_Bound_Surface_Display_Mode();
    void Cycle_Interaction_Pair_Display_Mode();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D,Cycle_Forces_Mode,"Cycle forces display mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D,Cycle_Hard_Bound_Surface_Display_Mode,"Cycle hard bound embedded surface display mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D,Cycle_Interaction_Pair_Display_Mode,"Cycle display of interaction pairs");
    
    TRIANGULATED_SURFACE<T>& Create_Hard_Bound_Boundary_Surface(TRIANGULATED_SURFACE<T>& boundary_surface);
    virtual void Reinitialize(bool force=false,bool read_geometry=true) PHYSBAM_OVERRIDE;    // Needs to be called after some state changes
};
}
#endif
