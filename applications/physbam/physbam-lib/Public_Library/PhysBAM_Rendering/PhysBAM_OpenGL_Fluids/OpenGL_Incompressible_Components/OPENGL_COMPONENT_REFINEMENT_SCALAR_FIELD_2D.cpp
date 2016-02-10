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
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Fluids/OpenGL_Incompressible_Components/OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D.h>
using namespace PhysBAM;

template<class T,class T2,class RW> OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D(GRID<TV> &grid_input,const std::string &sub_grids_filename_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    : OPENGL_COMPONENT("Sub Scalar Fields 2D"), grid(grid_input),
     scalar_field_filename(scalar_field_filename_input), sub_grids_filename(sub_grids_filename_input), frame_loaded(-1), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
    sub_grids.Resize(grid.Domain_Indices());
    scalar_fields.Resize(grid.Domain_Indices());
    opengl_scalar_field.Resize(grid.Domain_Indices());
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        opengl_scalar_field(iterator.Cell_Index())=new OPENGL_SCALAR_FIELD_2D<T,T2>(sub_grids(iterator.Cell_Index()),scalar_fields(iterator.Cell_Index()),color_map_input);
    Reinitialize();
}

template<class T,class T2,class RW> OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D(GRID<TV> &grid_input,const std::string &sub_grids_filename_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,
                                 typename OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_MODE draw_mode_input)
    : OPENGL_COMPONENT("Sub Scalar Fields 2D"), grid(grid_input),
      scalar_field_filename(scalar_field_filename_input), sub_grids_filename(sub_grids_filename_input), frame_loaded(-1), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
    sub_grids.Resize(grid.Domain_Indices());
    scalar_fields.Resize(grid.Domain_Indices());
    opengl_scalar_field.Resize(grid.Domain_Indices());
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        opengl_scalar_field(iterator.Cell_Index())=new OPENGL_SCALAR_FIELD_2D<T,T2>(sub_grids(iterator.Cell_Index()),scalar_fields(iterator.Cell_Index()),color_map_input);
    Reinitialize();
}

template<class T,class T2,class RW> OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
~OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        delete opengl_scalar_field(iterator.Cell_Index());
}

template<class T,class T2,class RW> bool OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename, frame_input);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw) for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_scalar_field(iterator.Cell_Index())->Display(in_color);
}

template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return RANGE<VECTOR<float,3> >((float)grid.domain.min_corner.x,(float)grid.domain.max_corner.x,(float)grid.domain.min_corner.y,(float)grid.domain.max_corner.y,0,0);
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        //opengl_scalar_field.Print_Selection_Info(output_stream,current_selection);
    }
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Reinitialize()
{
    if (draw)
    {
        if ((is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;
            std::string filename = FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename, frame);
            std::string grid_filename = FILE_UTILITIES::Get_Frame_Filename(sub_grids_filename, frame);
            if (FILE_UTILITIES::File_Exists(filename) && FILE_UTILITIES::File_Exists(grid_filename)){
                FILE_UTILITIES::Read_From_File<RW>(filename,scalar_fields);
                FILE_UTILITIES::Read_From_File<RW>(grid_filename,sub_grids);}
            else
                return;
            Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Set_Uniform_Contour_Values(const T2 min_value,const T2 max_value,const T2 increment)
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_scalar_field(iterator.Cell_Index())->Set_Uniform_Contour_Values(min_value,max_value,increment);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Update()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_scalar_field(iterator.Cell_Index())->Update();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Smooth()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_scalar_field(iterator.Cell_Index())->Toggle_Smooth_Texture();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Draw_Mode()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_scalar_field(iterator.Cell_Index())->Toggle_Draw_Mode();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Color_Map()
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) opengl_scalar_field(iterator.Cell_Index())->Toggle_Color_Map();
}
template class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<float,int,float>;
template class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<float,bool,float>;
template class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<double,int,double>;
template class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<double,bool,double>;
template class OPENGL_COMPONENT_REFINEMENT_SCALAR_FIELD_2D<double,double,double>;
#endif
