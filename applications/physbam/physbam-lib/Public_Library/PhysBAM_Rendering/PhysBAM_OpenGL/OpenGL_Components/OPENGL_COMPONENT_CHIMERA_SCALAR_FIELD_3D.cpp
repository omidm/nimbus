//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
using namespace PhysBAM;

template<class T,class T2,class RW> OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D(ARRAY<GRID<TV> > &grid_input, const std::string &scalar_field_filename_set_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool is_moving_grid_input,const std::string rigid_grid_frame_filename_set_input)
    : OPENGL_COMPONENT("Scalar Field 3D"),scalar_field_filename_set(scalar_field_filename_set_input), frame_loaded(-1), valid(false),current_grid(1), draw_all_grids(true),display_in_bool(false),dynamically_set_scale_range(false),scale_maximum(1),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input)
{
    for(int i=1;i<=grid_input.Size();i++){
        OPENGL_SCALAR_FIELD_3D<T,T2>* tmp_opengl_scalar_field=new OPENGL_SCALAR_FIELD_3D<T,T2>(grid_input(i),*new ARRAY<T2,VECTOR<int,3> >,color_map_input,OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_TEXTURE,is_moving_grid_input);
        OPENGL_SCALAR_FIELD_3D<T,bool>* tmp_opengl_bool_field=new OPENGL_SCALAR_FIELD_3D<T,bool>(grid_input(i),*new ARRAY<bool,VECTOR<int,3> >,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_TEXTURE,is_moving_grid_input);
        opengl_scalar_fields.Append(tmp_opengl_scalar_field);
        opengl_bool_fields.Append(tmp_opengl_bool_field);}
    if(opengl_scalar_fields.Size()>=1) opengl_scalar_field=opengl_scalar_fields(1);
    if(opengl_bool_fields.Size()>=1) opengl_bool_field=opengl_bool_fields(1);
    is_animation = true; 
}

template<class T,class T2,class RW> OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D(ARRAY<GRID<TV> > &grid_input, const std::string &scalar_field_filename_set_input,OPENGL_COLOR_MAP<T2>* color_map_input,typename OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_MODE draw_mode_input,bool is_moving_grid_input,const std::string rigid_grid_frame_filename_set_input)
    : OPENGL_COMPONENT("Scalar Field 3D"), scalar_field_filename_set(scalar_field_filename_set_input), frame_loaded(-1), valid(false), current_grid(1),draw_all_grids(true),display_in_bool(false),dynamically_set_scale_range(false),scale_maximum(1),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input)
{
    for(int i=1;i<=grid_input.Size();i++){
        OPENGL_SCALAR_FIELD_3D<T,T2>* tmp_opengl_scalar_field=new OPENGL_SCALAR_FIELD_3D<T,T2>(grid_input(i),*new ARRAY<T2,VECTOR<int,3> >,color_map_input,draw_mode_input,is_moving_grid_input);
        OPENGL_SCALAR_FIELD_3D<T,bool>* tmp_opengl_bool_field=new OPENGL_SCALAR_FIELD_3D<T,bool>(grid_input(i),*new ARRAY<bool,VECTOR<int,3> >,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS,is_moving_grid_input);
        opengl_scalar_fields.Append(tmp_opengl_scalar_field);
        opengl_bool_fields.Append(tmp_opengl_bool_field);}
    if(opengl_scalar_fields.Size()>=1) opengl_scalar_field=opengl_scalar_fields(1);
    if(opengl_bool_fields.Size()>=1) opengl_bool_field=opengl_bool_fields(1);
    is_animation = true;
}


template<class T,class T2,class RW> OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
~OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D()
{
    for(int i=1;i<=opengl_scalar_fields.Size();i++){
        delete &opengl_scalar_fields(i)->values;
        delete &opengl_bool_fields(i)->values;}
}

template<class T,class T2,class RW> bool OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(scalar_field_filename_set.c_str(),frame_input,1));
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw){
        if(!display_in_bool){
            if(draw_all_grids){
                for(int i=1;i<=opengl_scalar_fields.Size();i++) opengl_scalar_fields(i)->Display(in_color);}
            else opengl_scalar_field->Display(in_color);}
        else Display_In_Bool(in_color);}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Display_In_Bool(const int in_color) const
{
    if(draw_all_grids)
        for(int i=1;i<=opengl_bool_fields.Size();i++) opengl_bool_fields(i)->Display(in_color);
    else opengl_bool_field->Display(in_color);
}

template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_scalar_field->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_scalar_field->Print_Selection_Info(output_stream,current_selection);}
}

//for int, float, double
template<class T,class T2,class RW> template<class T3> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Set_Min_Max_Limits(WORKAROUND_STRUCT<T3>& ws, T2& max_value,T2& min_value)
{
    max_value=-std::numeric_limits<T2>::infinity();
    min_value=std::numeric_limits<T2>::infinity();
}

//for bool only
template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Set_Min_Max_Limits(WORKAROUND_STRUCT<bool>& ws, T2& max_value,T2& min_value)
{
    max_value=0;min_value=1;
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Reinitialize()
{
    if (draw){
        if ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            valid = false;
            if(is_moving_grid){
                for(int i=1;i<=opengl_scalar_fields.Size();i++){
                    std::string grid_filename=STRING_UTILITIES::string_sprintf(rigid_grid_frame_filename_set.c_str(),frame,i);
                    if(FILE_UTILITIES::File_Exists(grid_filename)){
                        FILE_UTILITIES::Read_From_File<T>(grid_filename,opengl_scalar_fields(i)->rigid_grid_frame);
                        FILE_UTILITIES::Read_From_File<T>(grid_filename,opengl_bool_fields(i)->rigid_grid_frame);}}}
            
            std::string filename = STRING_UTILITIES::string_sprintf(scalar_field_filename_set.c_str(),frame,1);
            if (FILE_UTILITIES::File_Exists(filename)){
                for(int i=1;i<=opengl_scalar_fields.Size();i++){
                    std::string filename_i = STRING_UTILITIES::string_sprintf(scalar_field_filename_set.c_str(),frame,i);
                    FILE_UTILITIES::Read_From_File<RW>(filename_i,opengl_scalar_fields(i)->values);}
                    if(dynamically_set_scale_range)
                        for(int i=1;i<=opengl_scalar_fields.Size();i++){
                            WORKAROUND_STRUCT<T2> ws;
                            T2 max_value,min_value;Set_Min_Max_Limits(ws,max_value,min_value);
                            for(CELL_ITERATOR iterator(opengl_scalar_fields(i)->grid);iterator.Valid();iterator.Next()){
                                if(opengl_scalar_fields(i)->values(iterator.Cell_Index())>max_value) max_value=opengl_scalar_fields(i)->values(iterator.Cell_Index());
                                if(opengl_scalar_fields(i)->values(iterator.Cell_Index())<min_value) min_value=opengl_scalar_fields(i)->values(iterator.Cell_Index());}
                            opengl_scalar_fields(i)->Set_Scale_Range(min_value,max_value);}
                    for(int i=1;i<=opengl_scalar_fields.Size();i++) opengl_scalar_fields(i)->Update();}
            else return;

            frame_loaded = frame;
            valid = true;
            
            if(display_in_bool)
                for(int i=1;i<=opengl_scalar_fields.Size();i++){
                    opengl_bool_fields(i)->values.Resize(opengl_scalar_fields(i)->grid.Domain_Indices(3));
                    opengl_scalar_fields(i)->values.Resize(opengl_scalar_fields(i)->grid.Domain_Indices(3));
                    for(CELL_ITERATOR iterator(opengl_scalar_fields(i)->grid,1);iterator.Valid();iterator.Next()){
                        opengl_bool_fields(i)->values(iterator.Cell_Index())=false;
                        opengl_bool_fields(i)->values(iterator.Cell_Index())=(opengl_scalar_fields(i)->values(iterator.Cell_Index())==0?false:true);}
                    opengl_bool_fields(i)->Update();}}}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Smooth_Slice()
{
    opengl_scalar_field->Toggle_Smooth_Slice_Texture();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Draw_Mode()
{
    opengl_scalar_field->Toggle_Draw_Mode();
    opengl_bool_field->Toggle_Draw_Mode();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Color_Map()
{
    opengl_scalar_field->Toggle_Color_Map();
    opengl_bool_field->Toggle_Color_Map();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Next_Grid()
{
    current_grid=min(current_grid+1,opengl_scalar_fields.m);
    std::stringstream ss;
    ss<<"viewing "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
    opengl_scalar_field=opengl_scalar_fields(current_grid);
    opengl_bool_field=opengl_bool_fields(current_grid);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Previous_Grid()
{
    current_grid=max(current_grid-1,1);
    std::stringstream ss;
    ss<<"viewing "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
    opengl_scalar_field=opengl_scalar_fields(current_grid);
    opengl_bool_field=opengl_bool_fields(current_grid);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Draw_All_Grids()
{
    draw_all_grids=!draw_all_grids;
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Set_To_Display_In_Bool_Mode()
{
    display_in_bool=true;
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Turn_On_Dynamically_Set_Scale_Range()
{
    dynamically_set_scale_range=true;
}


template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Increase_Scale_Range()
{
    WORKAROUND_STRUCT<T2> ws;
    Increase_Scale_Range_Helper(ws);
}

//Increase_Scale_Range for T2 = int and float
template<class T,class T2,class RW> template<class T3> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Increase_Scale_Range_Helper(WORKAROUND_STRUCT<T3>& ws)
{
    scale_maximum*=(T)1.1;
    for(int i=1;i<=opengl_scalar_fields.Size();i++){
        opengl_scalar_fields(i)->Set_Scale_Range(0,(T2)scale_maximum);
        opengl_scalar_fields(i)->Update();
    }
}

//Increase_Scale_Range for T2 = bool
template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Increase_Scale_Range_Helper(WORKAROUND_STRUCT<bool>& ws)
{
    scale_maximum*=(T)1.1;
    for(int i=1;i<=opengl_scalar_fields.Size();i++){
        opengl_scalar_fields(i)->Set_Scale_Range(0,scale_maximum!=0);
        opengl_scalar_fields(i)->Update();
    }
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Decrease_Scale_Range()
{
    WORKAROUND_STRUCT<T2> ws;
    Decrease_Scale_Range_Helper(ws);
}

//Decrease_Scale_Range for T2 = int and float
template<class T,class T2,class RW> template<class T3> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Decrease_Scale_Range_Helper(WORKAROUND_STRUCT<T3>& ws)
{
    scale_maximum/=(T)1.1;
    for(int i=1;i<=opengl_scalar_fields.Size();i++){
        opengl_scalar_fields(i)->Set_Scale_Range(0,(T2)scale_maximum);
        opengl_scalar_fields(i)->Update();
    }
}

//Decrease_Scale_Range for T2 = bool
template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Decrease_Scale_Range_Helper(WORKAROUND_STRUCT<bool>& ws)
{
    scale_maximum/=(T)1.1;
    for(int i=1;i<=opengl_scalar_fields.Size();i++){
        opengl_scalar_fields(i)->Set_Scale_Range(0,scale_maximum!=0);
        opengl_scalar_fields(i)->Update();
    }
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<T,T2,RW>::
Set_Scale_Range(const T2 min,const T2 max)
{
    dynamically_set_scale_range=false;
    scale_maximum=(T)max;
    for(int i=1;i<=opengl_scalar_fields.Size();i++){
        opengl_scalar_fields(i)->Set_Scale_Range(min,max);
        opengl_scalar_fields(i)->Update();
    }  
}


template class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<float,int,float>;
template class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<float,bool,float>;
template class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<double,int,double>;
template class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<double,bool,double>;
template class OPENGL_COMPONENT_CHIMERA_SCALAR_FIELD_3D<double,double,double>;
#endif
