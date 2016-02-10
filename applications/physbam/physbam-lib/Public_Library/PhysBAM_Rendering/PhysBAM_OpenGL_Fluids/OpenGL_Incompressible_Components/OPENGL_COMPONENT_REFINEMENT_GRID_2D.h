//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_REFINEMENT_GRID_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_REFINEMENT_GRID_2D__
#define __OPENGL_COMPONENT_REFINEMENT_GRID_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{
template<class TV> class GRID;

template<class T,class RW=T>
class OPENGL_COMPONENT_REFINEMENT_GRID_2D:public OPENGL_COMPONENT
{
private:
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
public:
    ARRAY<OPENGL_GRID_2D<T>*,TV_INT> opengl_grids;
    ARRAY<GRID<TV>,TV_INT> sub_grids;
private:
    GRID<TV> grid;
    std::string filename;
    int frame_loaded;
    bool valid;
public:

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw && valid;}

//#####################################################################
    OPENGL_COMPONENT_REFINEMENT_GRID_2D(const GRID<TV> &grid,const std::string& filename_input);
    ~OPENGL_COMPONENT_REFINEMENT_GRID_2D();
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

private:
    void Reinitialize();
//#####################################################################
};
}
#endif
