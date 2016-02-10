//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_BODY_MOTION_SEQUENCE
//#####################################################################
#ifndef __OPENGL_COMPONENT_BODY_MOTION_SEQUENCE__
#define __OPENGL_COMPONENT_BODY_MOTION_SEQUENCE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D.h>
#include <PhysBAM_Dynamics/Motion/BODY_MOTION_SEQUENCE.h>

namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_BODY_MOTION_SEQUENCE:public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_COMPONENT_BODY_MOTION_SEQUENCE(const std::string& filename_input,const bool frame_dependent_data_input,const std::string& rigid_body_input,const T scale_input=1);
    ~OPENGL_COMPONENT_BODY_MOTION_SEQUENCE();

    void Save_Frame();
    int Get_Frame();
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Prev_Frame() PHYSBAM_OVERRIDE;
    void Next_Frame() PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid && opengl_component_rigid_body_collection.Use_Bounding_Box(); }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Reinitialize(bool force=false);

public:
    std::string filename,rigid_filename;
    HASHTABLE<int,int> id_to_index;
    BODY_MOTION_SEQUENCE<T> body_motion;
    RIGID_BODY_COLLECTION<TV> rigid_body_collection;
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T> opengl_component_rigid_body_collection;
    FRAME<TV> rigid_body_base_transform;
    T scale,ui_scale;

protected:
    int frame_loaded;
    bool valid;
    bool frame_dependent_data;
    T frame_rate,one_over_frame_rate;
    T default_length;

//#####################################################################
};
}
#endif
