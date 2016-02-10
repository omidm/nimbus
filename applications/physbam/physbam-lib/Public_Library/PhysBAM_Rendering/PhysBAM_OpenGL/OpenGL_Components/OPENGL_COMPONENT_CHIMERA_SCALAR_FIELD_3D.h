//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D__
#define __OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CHIMERA_SLICE.h>
#include <string>

namespace PhysBAM
{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    template<class T_SPECIAL> struct WORKAROUND_STRUCT{};
public:
    // Should be able to combine these two constructors into one (with a default arg) but for some reason I can't get it to compile in linux...
    OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D(ARRAY<GRID<TV> > &grid_input,const std::string &scalar_field_filename_set_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool is_moving_grid_input=false,const std::string rigid_grid_frame_filename_set_input="");
    OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D(ARRAY<GRID<TV> > &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,typename OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_MODE draw_mode_input,bool is_moving_grid_input=false,const std::string rigid_grid_frame_filename_set_input="");
    ~OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Display_In_Bool(const int in_color=1) const;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    
    void Set_Chimera_Slice(OPENGL_CHIMERA_SLICE<T>* slice_input)
    {chimera_slice=slice_input;for(int i=1;i<=opengl_scalar_fields.m;i++){opengl_scalar_fields(i)->Set_Slice(slice_input->opengl_uniform_slices(i));opengl_bool_fields(i)->Set_Slice(slice_input->opengl_uniform_slices(i));}}

    void Slice_Has_Changed() PHYSBAM_OVERRIDE 
    {if(draw) for(int i=1;i<=opengl_scalar_fields.m;i++){opengl_scalar_fields(i)->Slice_Has_Changed();opengl_bool_fields(i)->Slice_Has_Changed();}}

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Toggle_Smooth_Slice();
    void Toggle_Draw_Mode();
    void Toggle_Color_Map();
    void Next_Grid();
    void Previous_Grid();
    void Toggle_Draw_All_Grids();
    void Set_To_Display_In_Bool_Mode();
    void Turn_On_Dynamically_Set_Scale_Range();
    void Set_Scale_Range(const T2 min,const T2 max);
    void Increase_Scale_Range();
    template<class T3> void Increase_Scale_Range_Helper(WORKAROUND_STRUCT<T3>& ws);
    void Increase_Scale_Range_Helper(WORKAROUND_STRUCT<bool>& ws);
    void Decrease_Scale_Range();
    template<class T3> void Decrease_Scale_Range_Helper(WORKAROUND_STRUCT<T3>& ws);
    void Decrease_Scale_Range_Helper(WORKAROUND_STRUCT<bool>& ws);
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Toggle_Smooth_Slice, "Toggle smooth slice");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Toggle_Draw_Mode, "Toggle draw mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Toggle_Color_Map, "Toggle color map");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Next_Grid, "Next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Previous_Grid, "Previous grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Toggle_Draw_All_Grids, "Toggle draw all grids");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Increase_Scale_Range,"Increase scale range");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D, Decrease_Scale_Range,"Decrease scale range");

private:
    void Reinitialize();
    template<class T3> void Set_Min_Max_Limits(WORKAROUND_STRUCT<T3>& ws,T2& max_value,T2& min_value);
    void Set_Min_Max_Limits(WORKAROUND_STRUCT<bool>& ws,T2& max_value,T2& min_value);

public:
    OPENGL_SCALAR_FIELD_3D<T,T2>*  opengl_scalar_field;
    OPENGL_SCALAR_FIELD_3D<T,bool>*  opengl_bool_field;
    ARRAY<OPENGL_SCALAR_FIELD_3D<T,T2>* > opengl_scalar_fields;
    ARRAY<OPENGL_SCALAR_FIELD_3D<T,bool>* > opengl_bool_fields;
    OPENGL_CHIMERA_SLICE<T>* chimera_slice;

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
