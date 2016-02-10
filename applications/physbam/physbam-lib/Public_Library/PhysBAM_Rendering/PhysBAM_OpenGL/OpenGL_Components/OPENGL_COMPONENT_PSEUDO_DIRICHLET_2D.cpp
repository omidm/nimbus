//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D(const GRID<TV> &grid,const std::string &filename_input)
    :OPENGL_COMPONENT("Pseudo Dirichlet"),mac_grid(grid.Get_MAC_Grid()),velocity_scale((T)0.025),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Display(const int in_color) const
{
    OPENGL_COLOR point_color=OPENGL_COLOR::White();
    OPENGL_COLOR velocity_color=OPENGL_COLOR::Yellow();
    if (valid && draw){
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT);
        glDisable(GL_LIGHTING);glPointSize(5.0);
        T x_offset=(T)0.4*mac_grid.dX.x,y_offset=(T)0.4*mac_grid.dX.y;
        point_color.Send_To_GL_Pipeline();
        ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
        for(int i=1;i<=pseudo_dirichlet_cells.m;i++){
            const VECTOR<int,2>& cell_index=pseudo_dirichlet_cells(i).x;char face_mask=pseudo_dirichlet_cells(i).z;
            VECTOR<T,2> pos=mac_grid.X(cell_index);
            if(face_mask&0x01) OpenGL_Vertex(VECTOR<T,2>(pos.x-x_offset,pos.y),vertices);
            if(face_mask&0x02) OpenGL_Vertex(VECTOR<T,2>(pos.x+x_offset,pos.y),vertices);
            if(face_mask&0x04) OpenGL_Vertex(VECTOR<T,2>(pos.x,pos.y-y_offset),vertices);
            if(face_mask&0x08) OpenGL_Vertex(VECTOR<T,2>(pos.x,pos.y+y_offset),vertices);}
        OpenGL_Draw_Arrays(GL_POINTS,2,vertices);
        vertices.Resize(0);
        velocity_color.Send_To_GL_Pipeline();
        for(int i=1;i<=pseudo_dirichlet_cells.m;i++){
            const VECTOR<int,2>& cell_index=pseudo_dirichlet_cells(i).x;const VECTOR<T,2>& velocity=pseudo_dirichlet_cells(i).y;
            VECTOR<T,2> pos=mac_grid.X(cell_index);
            OPENGL_SHAPES::Draw_Arrow(pos,pos+velocity_scale*velocity,vertices);}
        OpenGL_Draw_Arrays(GL_LINES,2,vertices);
        glPopAttrib();
    }
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Reinitialize(bool force)
{
    if (draw||force)
    {
        if (!valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(filename, frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,pseudo_dirichlet_cells);
            else
                return;

            frame_loaded=frame;
            valid=true;
        }
    }
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return RANGE<VECTOR<float,3> >((float)mac_grid.domain.min_corner.x,(float)mac_grid.domain.min_corner.y,0,(float)mac_grid.domain.max_corner.x,(float)mac_grid.domain.max_corner.y,0);
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Set_Vector_Size(const T vector_size)
{
    velocity_scale=vector_size;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Increase_Vector_Size()
{
    velocity_scale*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T,RW>::
Decrease_Vector_Size()
{
    velocity_scale/=(T)1.1;
}
template class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<double,double>;
#endif
