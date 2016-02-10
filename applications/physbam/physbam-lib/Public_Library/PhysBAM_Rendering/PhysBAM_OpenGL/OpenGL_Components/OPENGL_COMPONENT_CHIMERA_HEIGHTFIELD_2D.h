//#####################################################################
// Copyright 2004-2011, Eran Guendelman, Linhai Qiu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D__
#define __OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
public:
    OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D(ARRAY<GRID<TV> > &chimera_grid_input, 
                                    const std::string& height_filename_input,
                                    bool height_is_levelset_input=false);
    ~OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    void Set_Scale(T scale_input);

    void Increase_Scale();
    void Decrease_Scale();
    void Increase_Displacement_Scale();
    void Decrease_Displacement_Scale();
    void Toggle_Subdivision();
    void Toggle_Draw_All_Grids();
    void Next_Grid();
    void Previous_Grid();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Increase_Scale, "Increase scale");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Decrease_Scale, "Decrease scale");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Increase_Displacement_Scale, "Increase displacement scale");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Decrease_Displacement_Scale, "Decrease displacement scale");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Toggle_Subdivision, "Toggle subdivision");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Toggle_Draw_All_Grids, "Toggle draw all grids");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Next_Grid, "Next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D, Previous_Grid, "Previous grid");

public:
    void Reinitialize(bool force=false);
private:
    void Update_Surface();

    int To_Linear_Index(int grid_index, int i, int j) const
    {return (i-domain_array(grid_index).min_corner.x) + (j-domain_array(grid_index).min_corner.y)*counts_array(grid_index).x + 1;}

    VECTOR<int,2> From_Linear_Index(int grid_index,int idx) const
    {return VECTOR<int,2>((idx-1) % counts_array(grid_index).x + domain_array(grid_index).min_corner.x, ((idx-1)/counts_array(grid_index).x)+domain_array(grid_index).min_corner.y);}

public:
    ARRAY<GRID<TV> >& chimera_grid;
    ARRAY<TRIANGULATED_SURFACE<T>*> triangulated_surface_array;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> opengl_triangulated_surface_array;
    T vertical_offset;
    bool allow_smooth_shading;
    bool subdivide_surface;
    ARRAY<ARRAY<T,VECTOR<int,2> >*> height_array;
    ARRAY<LEVELSET_2D<GRID<TV> >*> levelset_array;

private:
    std::string height_filename;
    T scale;
    T displacement_scale;
    int frame_loaded;
    bool valid;
    ARRAY<RANGE<VECTOR<int,2> > > domain_array;
    ARRAY<VECTOR<int,2> > counts_array;
    bool height_is_levelset;
    int current_grid;
    bool draw_all_grids;
};

template<class T>
class OPENGL_SELECTION_COMPONENT_CHIMERA_HEIGHTFIELD_2D : public OPENGL_SELECTION
{
public:
    VECTOR<int,2> index;

    OPENGL_SELECTION_COMPONENT_CHIMERA_HEIGHTFIELD_2D(OPENGL_OBJECT *object) : OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_HEIGHTFIELD_2D, object) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
