//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CHIMERA_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
using namespace PhysBAM;
//#####################################################################

template<class T,class RW> OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
OPENGL_COMPONENT_CHIMERA_LEVELSET_2D(const std::string filename_set_input,bool is_moving_grid_input,const std::string rigid_grid_frame_filename_set_input)
    :OPENGL_COMPONENT("Chimera Levelset 2D"),opengl_levelset(0),filename_set(filename_set_input),frame_loaded(-1),current_grid(1),valid(false),draw_all_grids(true),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input)
{
    int number_of_grids=0;
    while(filename_set!=""){
        std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,number_of_grids+1);
        if(FILE_UTILITIES::File_Exists(filename)) number_of_grids++;else break;}
    LOG::cout<<"Found "<<number_of_grids<<" grids."<<std::endl;

    opengl_levelsets.Resize(number_of_grids);
    for(int j=1;j<=opengl_levelsets.m;j++)
        opengl_levelsets(j)=new OPENGL_LEVELSET_2D<T>(*(new LEVELSET_2D<GRID<TV> >(*(new GRID<TV>),*(new ARRAY<T,VECTOR<int,2> >))),OPENGL_COLOR::Blue(),OPENGL_COLOR::Red((T).5),0,true);
    opengl_levelset=opengl_levelsets(1);

    is_animation=true;
    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
~OPENGL_COMPONENT_CHIMERA_LEVELSET_2D()
{
    for(int j=1;j<=opengl_levelsets.m;j++){
        delete &opengl_levelsets(j)->levelset.grid;
        delete &opengl_levelsets(j)->levelset.phi;
        delete &opengl_levelsets(j)->levelset;}
}

template<class T,class RW> bool OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(filename_set.c_str(),current_grid,frame_input));
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Display(const int in_color) const
{
    if(valid && draw){
        if(draw_all_grids) for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Display(in_color);
        else opengl_levelset->Display(in_color);}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_levelset->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D && opengl_levelset->levelset.grid.Is_MAC_Grid()){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;
            output_stream<<component_name<<": phi="<<opengl_levelset->levelset.phi(index)
                         <<" curvature="<<opengl_levelset->levelset.Compute_Curvature(opengl_levelset->levelset.grid.Center(index))<<std::endl;}
        if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_2D && !opengl_levelset->levelset.grid.Is_MAC_Grid()){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_NODE_2D<T>*)current_selection)->index;
            output_stream<<component_name<<": phi="<<opengl_levelset->levelset.phi(index)<<std::endl;}
        if(current_selection && current_selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_2D){
            OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)current_selection;
            VECTOR<T,2> location=selection->location;
            output_stream<<component_name<<": phi @ particle="<<opengl_levelset->levelset.Phi(location)<<std::endl;}}
}


template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Reinitialize(const bool force_even_if_not_drawn)
{
    if (draw||force_even_if_not_drawn){
        if ((is_animation && (frame_loaded!=frame)) || (!is_animation && frame_loaded<0)){
            if(is_moving_grid){
                for(int i=1;i<=opengl_levelsets.Size();i++){
                    std::string grid_filename=STRING_UTILITIES::string_sprintf(rigid_grid_frame_filename_set.c_str(),frame,i);
                    if(FILE_UTILITIES::File_Exists(grid_filename)){
                        FILE_UTILITIES::Read_From_File<T>(grid_filename,opengl_levelsets(i)->rigid_grid_frame);}}}
            valid=false;std::string filename;
            for(int i=1;i<=opengl_levelsets.m;i++){
                filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,i);
                if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename.c_str(),opengl_levelsets(i)->levelset);
                else return;
                opengl_levelsets(i)->Update();}
            frame_loaded=frame;valid=true;}}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Toggle_Color_Mode()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Color_Map();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Toggle_Smooth()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Smooth_Texture();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Toggle_Normals()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Normals();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Toggle_Draw_Mode()
{
    for(int j=1;j<=opengl_levelsets.m;j++){
        int mask=((int)opengl_levelsets(j)->draw_area<<2) | ((int)opengl_levelsets(j)->draw_curve<<1) | ((int)opengl_levelsets(j)->draw_cells);
        int newmask=(mask%8)+1;
        opengl_levelsets(j)->draw_area=(newmask&4)!=0;
        opengl_levelsets(j)->draw_curve=(newmask&2)!=0;
        opengl_levelsets(j)->draw_cells=(newmask&1)!=0;
        opengl_levelsets(j)->Update();}
}
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Toggle_Draw_Sign()
{
    for(int j=1;j<=opengl_levelsets.m;j++){
        opengl_levelsets(j)->dominant_sign=(opengl_levelsets(j)->dominant_sign==1)?-1:1;
        opengl_levelsets(j)->Update();}
}
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Next_Grid()
{
    current_grid=min(current_grid+1,opengl_levelsets.m);
    LOG::cout<<"viewing levelset "<<current_grid<<std::endl;
    opengl_levelset=opengl_levelsets(current_grid);
}
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Previous_Grid()
{
    current_grid=max(current_grid-1,1);
    LOG::cout<<"viewing levelsetset "<<current_grid<<std::endl;
    opengl_levelset=opengl_levelsets(current_grid);
}
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Toggle_Draw_All_Grids()
{
    draw_all_grids=!draw_all_grids;
}
template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<T,RW>::
Toggle_Draw_Ghost_Values()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Draw_Ghost_Values();
}

template class OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CHIMERA_LEVELSET_2D<double,double>;
#endif
