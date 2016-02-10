//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SCALAR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_SCALAR_FIELD_2D__
#define __OPENGL_COMPONENT_SCALAR_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>

namespace PhysBAM
{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_SCALAR_FIELD_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    // Should be able to combine these two constructors into one (with a default arg) but for some reason I can't get it to compile in linux...
    OPENGL_COMPONENT_SCALAR_FIELD_2D(GRID<TV> &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool is_moving_grid_input=false,const std::string rigid_grid_frame_filename_input="");
    OPENGL_COMPONENT_SCALAR_FIELD_2D(GRID<TV> &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,typename OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_MODE draw_mode_input,bool is_moving_grid_input=false,const std::string rigid_grid_frame_filename_input="");
    ~OPENGL_COMPONENT_SCALAR_FIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Weights_I(std::ostream& stream,const TV_INT& cell) const;
    void Print_Weights_J(std::ostream& stream,const TV_INT& cell) const;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Toggle_Smooth();
    void Toggle_Draw_Mode();
    void Toggle_Color_Map();
    void Toggle_Increase_Color_Map_Range();
    void Toggle_Decrease_Color_Map_Range();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_2D, Toggle_Smooth, "Toggle smooth");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_2D, Toggle_Draw_Mode, "Toggle draw mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_2D, Toggle_Color_Map, "Toggle color map");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_2D, Toggle_Increase_Color_Map_Range, "Toggle increase color map range");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_2D, Toggle_Decrease_Color_Map_Range, "Toggle decrease color map range");
private:
    void Reinitialize();

public:
    OPENGL_SCALAR_FIELD_2D<T,T2>  opengl_scalar_field;

private:
    std::string scalar_field_filename;
    int frame_loaded;
    bool valid;
    bool is_moving_grid;
    std::string rigid_grid_frame_filename;

public:
    std::string weights_i_filename,weights_j_filename;
    bool print_weights_i,print_weights_j;
    ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT> weights_to; //w ij
    ARRAY<ARRAY<PAIR<TV_INT,int > >,TV_INT> weights_from; //w ki
};

}

#endif
