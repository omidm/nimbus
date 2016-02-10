//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D__
#define __OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <string>

namespace PhysBAM
{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    template<class T_SPECIAL> struct WORKAROUND_STRUCT{};
public:
    // Should be able to combine these two constructors into one (with a default arg) but for some reason I can't get it to compile in linux...
    OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D(ARRAY<GRID<TV> > &grid_input,const std::string &scalar_field_filename_set_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool is_moving_grid_input=false,const std::string rigid_grid_frame_filename_set_input="");
    OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D(ARRAY<GRID<TV> > &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,typename OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_MODE draw_mode_input,bool is_moving_grid_input=false,const std::string rigid_grid_frame_filename_set_input="");
    ~OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Display_In_Bool(const int in_color=1) const;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Toggle_Smooth();
    void Toggle_Draw_Mode();
    void Toggle_Color_Map();
    void Next_Grid();
    void Previous_Grid();
    void Toggle_Draw_All_Grids();
    void Toggle_Draw_Ghost_Cells();
    void Set_To_Display_In_Bool_Mode();
    void Turn_On_Dynamically_Set_Scale_Range();
    void Set_Scale_Range(const T2 min,const T2 max);
    void Increase_Scale_Range();
    template<class T3> void Increase_Scale_Range_Helper(WORKAROUND_STRUCT<T3>& ws);
    void Increase_Scale_Range_Helper(WORKAROUND_STRUCT<bool>& ws);
    void Decrease_Scale_Range();
    template<class T3> void Decrease_Scale_Range_Helper(WORKAROUND_STRUCT<T3>& ws);
    void Decrease_Scale_Range_Helper(WORKAROUND_STRUCT<bool>& ws);
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Toggle_Smooth, "Toggle smooth");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Toggle_Draw_Mode, "Toggle draw mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Toggle_Color_Map, "Toggle color map");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Next_Grid, "Next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Previous_Grid, "Previous grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Toggle_Draw_All_Grids, "Toggle draw all grids");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Toggle_Draw_Ghost_Cells,"Toggle to draw ghost cells");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Increase_Scale_Range,"Increase scale range");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_2D, Decrease_Scale_Range,"Decrease scale range");

private:
    void Reinitialize();
    template<class T3> inline void Set_Min_Max_Limits(WORKAROUND_STRUCT<T3>& ws,T2& max_value,T2& min_value);
    inline void Set_Min_Max_Limits(WORKAROUND_STRUCT<bool>& ws,T2& max_value,T2& min_value);
    template<class T3> inline void Set_Min_Max_Values(WORKAROUND_STRUCT<T3>& ws,T2& max_value,T2& min_value,const T2& value_to_check);
    inline void Set_Min_Max_Values(WORKAROUND_STRUCT<bool>& ws,T2& max_value,T2& min_value,const T2& value_to_check);

public:
    OPENGL_SCALAR_FIELD_2D<T,T2>*  opengl_scalar_field;
    OPENGL_SCALAR_FIELD_2D<T,bool>*  opengl_bool_field;
    ARRAY<OPENGL_SCALAR_FIELD_2D<T,T2>* > opengl_scalar_fields;
    ARRAY<OPENGL_SCALAR_FIELD_2D<T,bool>* > opengl_bool_fields;

private:
    std::string scalar_field_filename_set;
    int frame_loaded;
    bool valid;
    int current_grid;
    bool draw_all_grids;
    bool display_in_bool;
    bool dynamically_set_scale_range;
    T scale_maximum;
    bool is_moving_grid;
    std::string rigid_grid_frame_filename_set;
};

}

#endif
