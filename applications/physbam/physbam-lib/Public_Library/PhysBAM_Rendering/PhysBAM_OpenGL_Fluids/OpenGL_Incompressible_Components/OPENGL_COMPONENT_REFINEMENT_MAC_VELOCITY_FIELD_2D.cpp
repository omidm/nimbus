//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Fluids/OpenGL_Incompressible_Components/OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D(const GRID<TV> &grid_input,const std::string &velocity_filename_input,const std::string &sub_grids_filename_input)
    : OPENGL_COMPONENT("Refinement MAC Velocity Field 2D"),grid(grid_input),velocity_filename(velocity_filename_input),
     sub_grids_filename(sub_grids_filename_input),level(1),level_loaded(-1),
     valid(false),draw_all_levels(true)
{
    Initialize();
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Initialize()
{
    is_animation=FILE_UTILITIES::Is_Animated(velocity_filename);
    frame_loaded=-1;

    // only 1 level for now
    number_of_levels=1;
    
    opengl_refinement_mac_velocity_fields.Resize(grid.Domain_Indices());
    sub_mac_velocities.Resize(grid.Domain_Indices(),false);
    sub_grids.Resize(grid.Domain_Indices());

    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        sub_grids(iterator.Cell_Index()).MAC_offset=0.5;
        opengl_refinement_mac_velocity_fields(iterator.Cell_Index())=new OPENGL_MAC_VELOCITY_FIELD_2D<T>(sub_grids(iterator.Cell_Index()),sub_mac_velocities(iterator.Cell_Index()));}

    OPENGL_COLOR_RAMP<T>* ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color((T)-1e+2,OPENGL_COLOR::Red());
    ramp->Add_Color(-1,OPENGL_COLOR::Yellow());
    ramp->Add_Color((T)-1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(0,OPENGL_COLOR::Black());
    ramp->Add_Color((T)1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(1,OPENGL_COLOR::Yellow());
    ramp->Add_Color((T)1e+2,OPENGL_COLOR::Red());

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
~OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        delete opengl_refinement_mac_velocity_fields(iterator.Cell_Index());
}

template<class T,class RW> bool OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(velocity_filename, frame_input);
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        //opengl_mac_velocity_field->Print_Selection_Info(stream,selection);
    }
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Display(const int in_color) const
{
    if(valid){
        if(draw){
            //if(draw_all_levels) for(int i=1;i<=opengl_refinement_mac_velocity_fields.m;i++) opengl_refinement_mac_velocity_fields(i)->Display(in_color);
            //else opengl_refinement_mac_velocity_fields(level)->Display(in_color);}
            for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_refinement_mac_velocity_fields(iterator.Cell_Index())->Display(in_color);}}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return RANGE<VECTOR<float,3> >((float)grid.domain.min_corner.x,(float)grid.domain.max_corner.x,(float)grid.domain.min_corner.y,(float)grid.domain.max_corner.y,0,0);
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Reinitialize()
{
    if (draw){
        if ((is_animation && (frame_loaded!=frame || level_loaded!=level)) || (!is_animation && frame_loaded < 0)){
            valid = false;
            std::string tmp_filename=FILE_UTILITIES::Get_Frame_Filename(velocity_filename.c_str(), frame);
            std::string grid_filename = FILE_UTILITIES::Get_Frame_Filename(sub_grids_filename.c_str(), frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename) && FILE_UTILITIES::File_Exists(grid_filename)){
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,sub_mac_velocities);
                FILE_UTILITIES::Read_From_File<RW>(grid_filename,sub_grids);}
            else return;

            frame_loaded=frame;level_loaded=level;valid=true;}}
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Velocity_Mode()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_refinement_mac_velocity_fields(iterator.Cell_Index())->Toggle_Velocity_Mode();
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Velocity_Mode_And_Draw()
{
    if (draw)
    {
        Toggle_Velocity_Mode();
        if ((int)(opengl_refinement_mac_velocity_fields(TV_INT(1,1))->velocity_mode)==0) Toggle_Draw();
    }
    else Toggle_Draw();
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Increase_Vector_Size()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_refinement_mac_velocity_fields(iterator.Cell_Index())->Scale_Vector_Size((T)1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Decrease_Vector_Size()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_refinement_mac_velocity_fields(iterator.Cell_Index())->Scale_Vector_Size(1/(T)1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Arrowhead()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_refinement_mac_velocity_fields(iterator.Cell_Index())->draw_arrowhead = !opengl_refinement_mac_velocity_fields(iterator.Cell_Index())->draw_arrowhead;
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_All_Levels()
{
    draw_all_levels=!draw_all_levels;
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Next_Level(){
    level=min(level+1,number_of_levels);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Previous_Level(){
    level=max(level-1,1);
    Reinitialize();
}

template class OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D<double,double>;
#endif
