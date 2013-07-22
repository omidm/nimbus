//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
using namespace PhysBAM;
//#####################################################################

template<class T,class T2,class RW> OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D(const ARRAY<GRID<TV> > &grid_input, const std::string &values_filename_set_input, OPENGL_COLOR_MAP<T2>* color_map_input,bool is_moving_grid_input,std::string rigid_grid_frame_filename_set_input)
    :OPENGL_COMPONENT("Face Scalar Field 3D"),
    values_filename_set(values_filename_set_input), frame_loaded(-1), valid(false), current_grid(1),draw_all_grids(true),display_in_bool(false),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input)
{
    for(int i=1;i<=grid_input.Size();i++){
        OPENGL_FACE_SCALAR_FIELD_3D<T,T2>* tmp_opengl_scalar_field=new OPENGL_FACE_SCALAR_FIELD_3D<T,T2>(grid_input(i),*new ARRAY<T2,FACE_INDEX<3> >,color_map_input,is_moving_grid_input);
        OPENGL_FACE_SCALAR_FIELD_3D<T,bool>* tmp_opengl_bool_field=new OPENGL_FACE_SCALAR_FIELD_3D<T,bool>(grid_input(i),*new ARRAY<bool,FACE_INDEX<3> >,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()),is_moving_grid_input);
        opengl_scalar_fields.Append(tmp_opengl_scalar_field);
        opengl_bool_fields.Append(tmp_opengl_bool_field);}
    if(opengl_scalar_fields.Size()>=1) opengl_scalar_field=opengl_scalar_fields(1);
    if(opengl_bool_fields.Size()>=1) opengl_bool_field=opengl_bool_fields(1);
    is_animation = FILE_UTILITIES::Is_Animated(values_filename_set);
    Reinitialize();
}

template<class T,class T2,class RW> OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
~OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D()
{
    for(int i=1;i<=opengl_scalar_fields.Size();i++){
        delete &opengl_scalar_fields(i)->x_face_values;
        delete &opengl_scalar_fields(i)->y_face_values;
        delete &opengl_bool_fields(i)->x_face_values;
        delete &opengl_bool_fields(i)->y_face_values;}
}

template<class T,class T2,class RW> bool OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(values_filename_set.c_str(),frame_input,1));
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Display(const int in_color) const
{
    if(valid && draw){
        if(!display_in_bool){
            if(draw_all_grids) for(int i=1;i<=opengl_scalar_fields.Size();i++) opengl_scalar_fields(i)->Display(in_color);
            else opengl_scalar_field->Display(in_color);}
        else Display_In_Bool(in_color);}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Display_In_Bool(const int in_color) const
{
    if(draw_all_grids)
        for(int i=1;i<=opengl_bool_fields.Size();i++) opengl_bool_fields(i)->Display(in_color);
    else opengl_bool_field->Display(in_color);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_scalar_field->Print_Selection_Info(stream,selection);}
}

template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_scalar_field->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
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

            std::string filename;
            if(!values_filename_set.empty()){
                filename=STRING_UTILITIES::string_sprintf(values_filename_set.c_str(), frame, 1);
                if(FILE_UTILITIES::File_Exists(filename))
                    for(int i=1;i<=opengl_scalar_fields.Size();i++){
                        std::string filename_i=STRING_UTILITIES::string_sprintf(values_filename_set.c_str(), frame, i);
                        FILE_UTILITIES::Read_From_File<RW>(filename_i,opengl_scalar_fields(i)->face_values);
                        opengl_scalar_fields(i)->Update();}
                else return;}
            
            frame_loaded = frame;
            valid = true;

            for(int i=1;i<=opengl_scalar_fields.Size();i++){
                opengl_bool_fields(i)->face_values.Resize(opengl_scalar_fields(i)->grid.Domain_Indices(3));
                for(FACE_ITERATOR iterator(opengl_scalar_fields(i)->grid);iterator.Valid();iterator.Next()){
                    opengl_bool_fields(i)->face_values.Component(iterator.Axis())(iterator.Face_Index())=false;
                    opengl_bool_fields(i)->face_values.Component(iterator.Axis())(iterator.Face_Index())=(opengl_scalar_fields(i)->face_values.Component(iterator.Axis())(iterator.Face_Index())==0?false:true);}
                opengl_bool_fields(i)->Update();}}}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Next_Grid()
{
    current_grid=min(current_grid+1,opengl_scalar_fields.m);
    std::stringstream ss;
    ss<<"viewing "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
    opengl_scalar_field=opengl_scalar_fields(current_grid);
    opengl_bool_field=opengl_bool_fields(current_grid);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Previous_Grid()
{
    current_grid=max(current_grid-1,1);
    std::stringstream ss;
    ss<<"viewing "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
    opengl_scalar_field=opengl_scalar_fields(current_grid);
    opengl_bool_field=opengl_bool_fields(current_grid);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Draw_All_Grids()
{
    draw_all_grids=!draw_all_grids;
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<T,T2,RW>::
Set_To_Display_In_Bool_Mode()
{
    display_in_bool=true;
}

template class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<float,int,float>;
template class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<float,bool,float>;
template class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<double,int,double>;
template class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<double,bool,double>;
template class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D<double,double,double>;
#endif
