//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_2D.h>
using namespace PhysBAM;

template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
OPENGL_COMPONENT_SCALAR_FIELD_2D(GRID<TV> &grid_input, const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool is_moving_grid_input,const std::string rigid_grid_frame_filename_input)
    : OPENGL_COMPONENT("Scalar Field 2D"), opengl_scalar_field(grid_input,*new ARRAY<T2,VECTOR<int,2> >,color_map_input,OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_TEXTURE,is_moving_grid_input),scalar_field_filename(scalar_field_filename_input), frame_loaded(-1), valid(false),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename(rigid_grid_frame_filename_input),print_weights_i(false),print_weights_j(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
}

template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
OPENGL_COMPONENT_SCALAR_FIELD_2D(GRID<TV> &grid_input, const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,typename OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_MODE draw_mode_input,bool is_moving_grid_input,const std::string rigid_grid_frame_filename_input)
    : OPENGL_COMPONENT("Scalar Field 2D"), opengl_scalar_field(grid_input,*new ARRAY<T2,VECTOR<int,2> >,color_map_input,draw_mode_input,is_moving_grid_input),scalar_field_filename(scalar_field_filename_input), frame_loaded(-1), valid(false),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename(rigid_grid_frame_filename_input),print_weights_i(false),print_weights_j(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
}


template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
~OPENGL_COMPONENT_SCALAR_FIELD_2D()
{
    delete &opengl_scalar_field.values;
}

template<class T,class T2,class RW> bool OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename, frame_input);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_scalar_field.Display(in_color);
}

template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Print_Weights_I(std::ostream& stream,const VECTOR<int,2>& cell) const
{
    ARRAY<PAIR<TV_INT,T> >& local_weights=const_cast<ARRAY<PAIR<TV_INT,T> >&>(weights_to(cell));
    for(int i=1;i<=local_weights.m;i++){
        stream<<"    "<<local_weights(i).y<<" to "<<local_weights(i).x<<std::endl;}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Print_Weights_J(std::ostream& stream,const VECTOR<int,2>& cell) const
{
    ARRAY<PAIR<TV_INT,int> >& local_weights=const_cast<ARRAY<PAIR<TV_INT,int> >&>(weights_from(cell));
    for(int i=1;i<=local_weights.m;i++){
        stream<<"    "<<weights_to(local_weights(i).x)(local_weights(i).y).y<<" to "<<local_weights(i).x<<std::endl;}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_scalar_field.Print_Selection_Info(output_stream,current_selection);
        if(print_weights_i && ((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;
            output_stream<<"w_i"<<std::endl;
            Print_Weights_I(output_stream,index);}
        if(print_weights_j && ((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;
            output_stream<<"w_j"<<std::endl;
            Print_Weights_J(output_stream,index);}}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Reinitialize()
{
    if (draw){
        if ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            valid = false;
            if(print_weights_i){
                std::string tmp_i_filename=FILE_UTILITIES::Get_Frame_Filename(weights_i_filename.c_str(), frame);
                if (FILE_UTILITIES::File_Exists(tmp_i_filename)) FILE_UTILITIES::Read_From_File<RW>(tmp_i_filename,weights_to);
                else return;}
            if(print_weights_j){
                std::string tmp_j_filename=FILE_UTILITIES::Get_Frame_Filename(weights_j_filename.c_str(), frame);
                if (FILE_UTILITIES::File_Exists(tmp_j_filename)) FILE_UTILITIES::Read_From_File<RW>(tmp_j_filename,weights_from);
                else return;}
            if(is_moving_grid){
                std::string grid_filename=STRING_UTILITIES::string_sprintf(rigid_grid_frame_filename.c_str(),frame);
                if(FILE_UTILITIES::File_Exists(grid_filename)) FILE_UTILITIES::Read_From_File<T>(grid_filename,opengl_scalar_field.rigid_grid_frame);}
            std::string filename = FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename, frame);
            if (FILE_UTILITIES::File_Exists(filename))
                FILE_UTILITIES::Read_From_File<RW>(filename,opengl_scalar_field.values);
            else return;

            opengl_scalar_field.Update();
            frame_loaded = frame;
            valid = true;}}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Smooth()
{
    opengl_scalar_field.Toggle_Smooth_Texture();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Draw_Mode()
{
    opengl_scalar_field.Toggle_Draw_Mode();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Color_Map()
{
    opengl_scalar_field.Toggle_Color_Map();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Increase_Color_Map_Range()
{
    opengl_scalar_field.Increase_Scale_Range();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Decrease_Color_Map_Range()
{
    opengl_scalar_field.Decrease_Scale_Range();
}

template class OPENGL_COMPONENT_SCALAR_FIELD_2D<float,int,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<float,bool,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<double,int,double>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<double,bool,double>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<double,double,double>;
#endif
