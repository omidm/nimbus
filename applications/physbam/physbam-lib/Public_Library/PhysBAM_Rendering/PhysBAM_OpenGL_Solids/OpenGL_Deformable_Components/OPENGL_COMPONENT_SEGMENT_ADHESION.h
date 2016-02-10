//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SEGMENT_ADHESION
//#####################################################################
#ifndef __OPENGL_COMPONENT_SEGMENT_ADHESION__
#define __OPENGL_COMPONENT_SEGMENT_ADHESION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM{
template<class T,class RW> class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D;

template<class T,class RW=T>
class OPENGL_COMPONENT_SEGMENT_ADHESION:public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_COMPONENT_SEGMENT_ADHESION(const std::string& filename,const OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>& opengl_component_deformable_body,
        const OPENGL_COLOR& material_input=OPENGL_COLOR::Yellow());
    ~OPENGL_COMPONENT_SEGMENT_ADHESION();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid && opengl_segmented_curve.Use_Bounding_Box(); }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
protected:
    virtual void Reinitialize(bool force=false);

public:
    std::string filename;
    const OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>& opengl_component_deformable_body;
    PARTICLES<TV> particles;
    SEGMENT_MESH segment_mesh;
    SEGMENTED_CURVE<TV> segmented_curve;
    OPENGL_SEGMENTED_CURVE_3D<T> opengl_segmented_curve;

protected:
    int frame_loaded;
    bool valid;

//#####################################################################
};
}
#endif
