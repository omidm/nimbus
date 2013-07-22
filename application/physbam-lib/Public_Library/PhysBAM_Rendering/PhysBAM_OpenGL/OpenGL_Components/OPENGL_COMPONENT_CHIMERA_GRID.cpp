//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CHIMERA_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// OPENGL_COMPONENT_CHIMERA_GRID
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
OPENGL_COMPONENT_CHIMERA_GRID(const std::string filename_set_input,bool is_moving_grid_input,std::string rigid_grid_frame_filename_set_input):OPENGL_COMPONENT("Chimera Area"),filename_set(filename_set_input),frame_loaded(-1),valid(false),current_grid(1),draw_all_grids(true),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input)
{
    int number_of_grids=0;
    while(filename_set!=""){
        std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,number_of_grids+1);
        std::stringstream ss;
        ss<<"Checking "<<filename<<std::endl;
        LOG::filecout(ss.str());
        if(FILE_UTILITIES::File_Exists(filename)) number_of_grids++;else break;}
    std::stringstream ss;
    ss<<"Found "<<number_of_grids<<" levelsets for multiphase"<<std::endl;
    LOG::filecout(ss.str());

    opengl_grids.Resize(number_of_grids);
    for(int j=1;j<=opengl_grids.m;j++)
        opengl_grids(j)=new OPENGL_GRID_2D<T>(*(new GRID<TV>),OPENGL_COLOR::Gray(.5),"",0,is_moving_grid_input,rigid_grid_frame_filename_set_input);

    opengl_grid=opengl_grids(1);

    is_animation=true;
    Reinitialize();
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame_input,1));
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Display(const int in_color) const
{
    if(valid&&draw){
        if(draw_all_grids) for(int i=1;i<=opengl_grids.Size();i++) opengl_grids(i)->Display(in_color);
        else opengl_grid->Display(in_color);}
}

//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Bounding_Box() const
{
    LOG::cout << "bounding box range " << opengl_grid->Bounding_Box() << std::endl;
    if(valid&&draw)return opengl_grid->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    opengl_grid->Print_Selection_Info(output_stream,selection);
}
//#####################################################################
// Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    return opengl_grid->Create_Or_Destroy_Selection_After_Frame_Change(old_selection,delete_selection);
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Reinitialize()
{
    if((is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)){
        valid=false;
        if(is_moving_grid){
            for(int i=1;i<=opengl_grids.m;i++){
                std::string grid_filename=STRING_UTILITIES::string_sprintf(rigid_grid_frame_filename_set.c_str(),frame,i);
                if(FILE_UTILITIES::File_Exists(grid_filename)) FILE_UTILITIES::Read_From_File<T>(grid_filename,opengl_grids(i)->rigid_grid_frame);}}
        std::string filename;
        for(int i=1;i<=opengl_grids.m;i++){
            filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename)) 
                FILE_UTILITIES::Read_From_File<RW>(filename.c_str(),opengl_grids(i)->grid);
            else return;
            //opengl_grids(i)->draw_ghost_values=false;
        }
        frame_loaded=frame;valid=true;}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Next_Grid()
{
    current_grid=min(current_grid+1,opengl_grids.m);
    std::stringstream ss;
    ss<<"viewing grid "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
    opengl_grid=opengl_grids(current_grid);
}
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Previous_Grid()
{
    current_grid=max(current_grid-1,1);
    std::stringstream ss;
    ss<<"viewing grid "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
    opengl_grid=opengl_grids(current_grid);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Toggle_Draw_All_Grids()
{
    draw_all_grids=!draw_all_grids;
}
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_GRID<T,RW>::
Toggle_Draw_Ghost_Cells()
{
    for(int i=1;i<=opengl_grids.Size();i++) opengl_grids(i)->Toggle_Draw_Ghost_Values();
}

//#####################################################################
template class OPENGL_COMPONENT_CHIMERA_GRID<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CHIMERA_GRID<double,double>;
#endif

