//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS
//#####################################################################
#ifndef __OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS__
#define __OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_PARAMETERS.h>

namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS:public OPENGL_COMPONENT
{
public:
    std::string basedir;
    std::string filename;
    bool valid;
    int frame_loaded;

    //ARRAY<VECTOR<T,3> > X_dummy;
    //ARRAY<ARRAY<int> > attached_nodes_dummy;
    //ACTIVATION_CONTROL_SET<T> activation_control_set;
    //ATTACHMENT_FRAME_CONTROL_SET<T> attachment_frame_control_set;
    FACE_CONTROL_PARAMETERS<T> face_control_parameters;

    OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS(const std::string& basedir_input,const std::string& filename);
    ~OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS();

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return false;}

//#####################################################################
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Reinitialize(bool force=false);
//#####################################################################
};
}
#endif
