//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D__
#define __OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>

namespace PhysBAM
{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
public:
    OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D(GRID<TV> &grid_input,const std::string &sub_grids_filename_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input);
    OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D(GRID<TV> &grid_input,const std::string &sub_grids_filename_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,typename OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_MODE draw_mode_input);
    ~OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Set_Uniform_Contour_Values(const T2 min_avlue,const T2 max_value,const T2 increment);
    void Update();  // Call when values or other attributes have changed

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Toggle_Smooth();
    void Toggle_Draw_Mode();
    void Toggle_Color_Map();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D, Toggle_Smooth, "Toggle smooth");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D, Toggle_Draw_Mode, "Toggle draw mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D, Toggle_Color_Map, "Toggle color map");

private:
    void Reinitialize();

public:
    GRID<TV>& grid;
    ARRAY<GRID<TV>,TV_INT> sub_grids;
    ARRAY<OPENGL_SCALAR_FIELD_2D<T,T2>*,TV_INT>  opengl_scalar_field;
    ARRAY<ARRAY<T2,TV_INT>,TV_INT> scalar_fields;

private:
    std::string scalar_field_filename,sub_grids_filename;
    int frame_loaded;
    bool valid;

};

}

#endif
