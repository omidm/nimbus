//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Fluids/OpenGL_Incompressible_Components/OPENGL_COMPONENT_REFINEMENT_GRID_2D.h>
using namespace PhysBAM;
//#####################################################################
// OPENGL_COMPONENT_REFINEMENT_GRID_2D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
OPENGL_COMPONENT_REFINEMENT_GRID_2D(const GRID<TV> &grid_input,const std::string& filename_input):OPENGL_COMPONENT("Sub Grids"),grid(grid_input),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation=true;//FILE_UTILITIES::Is_Animated(filename);
    sub_grids.Resize(grid.Domain_Indices());
    opengl_grids.Resize(grid.Domain_Indices());
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        opengl_grids(iterator.Cell_Index())=new OPENGL_GRID_2D<T>(sub_grids(iterator.Cell_Index()));
        opengl_grids(iterator.Cell_Index())->draw_ghost_values=false;}
    Reinitialize();
}
//#####################################################################
// ~OPENGL_COMPONENT_REFINEMENT_GRID_2D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
~OPENGL_COMPONENT_REFINEMENT_GRID_2D()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) delete opengl_grids(iterator.Cell_Index());
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
Display(const int in_color) const
{
    if(valid&&draw) for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_grids(iterator.Cell_Index())->Display(in_color);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
Bounding_Box() const
{
    if(valid&&draw)return RANGE<VECTOR<float,3> >((float)grid.domain.min_corner.x,(float)grid.domain.max_corner.x,(float)grid.domain.min_corner.y,(float)grid.domain.max_corner.y,0,0);
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_REFINEMENT_GRID_2D<T,RW>::
Reinitialize()
{
    if((is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)){
        valid=false;
        std::string filename_string=FILE_UTILITIES::Get_Frame_Filename(filename, frame);
        if(FILE_UTILITIES::File_Exists(filename_string)){
            FILE_UTILITIES::Read_From_File<RW>(filename_string,sub_grids);}
        else return;
        frame_loaded=frame;valid=true;}
}
//#####################################################################
template class OPENGL_COMPONENT_REFINEMENT_GRID_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_REFINEMENT_GRID_2D<double,double>;
#endif
