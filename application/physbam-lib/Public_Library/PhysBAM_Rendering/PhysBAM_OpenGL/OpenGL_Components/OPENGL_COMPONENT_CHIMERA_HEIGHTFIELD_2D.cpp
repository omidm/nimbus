//#####################################################################
// Copyright 2004-2011, Eran Guendelman, Linhai Qiu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
using namespace PhysBAM;
//#####################################################################

template<class T,class RW> OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D(ARRAY<GRID<TV> > &chimera_grid_input, 
                                const std::string& height_filename_input,
                                bool height_is_levelset_input)
    :OPENGL_COMPONENT("Heightfield 2D"),chimera_grid(chimera_grid_input),
    triangulated_surface_array(chimera_grid_input.Size()),opengl_triangulated_surface_array(chimera_grid_input.Size()),
    vertical_offset(0), allow_smooth_shading(true), subdivide_surface(false), height_array(chimera_grid_input.Size()), levelset_array(chimera_grid_input.Size()),
    height_filename(height_filename_input),scale(1), displacement_scale(1), valid(false),domain_array(chimera_grid_input.Size()),
    counts_array(chimera_grid_input.Size()),height_is_levelset(height_is_levelset_input),current_grid(1),draw_all_grids(true)
{
    for(int grid_index=1;grid_index<=chimera_grid.Size();grid_index++){
        domain_array(grid_index)=chimera_grid(grid_index).Domain_Indices();
        counts_array(grid_index)=domain_array(grid_index).Edge_Lengths()+1;
        PHYSBAM_ASSERT(counts_array(grid_index).Min()>0);
        
        triangulated_surface_array(grid_index)=TRIANGULATED_SURFACE<T>::Create();
        opengl_triangulated_surface_array(grid_index)=new OPENGL_TRIANGULATED_SURFACE<T>(*triangulated_surface_array(grid_index),false);
        opengl_triangulated_surface_array(grid_index)->Set_Two_Sided();
        opengl_triangulated_surface_array(grid_index)->Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(0,0.6f,0.9f)));
        opengl_triangulated_surface_array(grid_index)->Set_Back_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(0.8f,0,0)));
        
        triangulated_surface_array(grid_index)->mesh.Initialize_Square_Mesh(counts_array(grid_index).x,counts_array(grid_index).y);
        for (int i=domain_array(grid_index).min_corner.x;i<=domain_array(grid_index).max_corner.x;i++) for(int j=domain_array(grid_index).min_corner.y;j<=domain_array(grid_index).max_corner.y;j++){
            triangulated_surface_array(grid_index)->particles.array_collection->Add_Element();}

        height_array(grid_index)=new ARRAY<T,VECTOR<int,2> >(chimera_grid(grid_index).Domain_Indices());
        levelset_array(grid_index)=new LEVELSET_2D<GRID<TV> >(chimera_grid(grid_index),*height_array(grid_index));}

    is_animation = FILE_UTILITIES::Is_Animated(height_filename);
    frame_loaded = -1;

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
~OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D()
{
    for(int grid_index=1;grid_index<=chimera_grid.Size();grid_index++){
        delete &(triangulated_surface_array(grid_index)->mesh);
        delete &(triangulated_surface_array(grid_index)->particles);
        delete triangulated_surface_array(grid_index);}
}

template<class T,class RW> bool OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(height_filename.c_str(),frame_input,current_grid));
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Display(const int in_color) const
{
    if(valid && draw){
        if(draw_all_grids)
            for(int grid_index=1;grid_index<=chimera_grid.Size();grid_index++) opengl_triangulated_surface_array(grid_index)->Display(grid_index);
        else opengl_triangulated_surface_array(current_grid)->Display(in_color);}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_triangulated_surface_array(current_grid)->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Turn_Smooth_Shading_On()
{
    if(allow_smooth_shading) for(int grid_index=1;grid_index<=chimera_grid.Size();grid_index++) opengl_triangulated_surface_array(grid_index)->Turn_Smooth_Shading_On();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Turn_Smooth_Shading_Off()
{
    for(int grid_index=1;grid_index<=chimera_grid.Size();grid_index++) opengl_triangulated_surface_array(grid_index)->Turn_Smooth_Shading_Off();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Reinitialize(bool force)
{
    if(draw){
        if(force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            bool success = true;
            valid = false;

            if(success){
                for(int grid_index=1;grid_index<=chimera_grid.Size();grid_index++){
                    std::string filename_tmp=STRING_UTILITIES::string_sprintf(height_filename.c_str(),frame,grid_index);
                    if(FILE_UTILITIES::File_Exists(filename_tmp)){
                        if(!height_is_levelset) FILE_UTILITIES::Read_From_File<RW>(filename_tmp,*height_array(grid_index));
                        else FILE_UTILITIES::Read_From_File<RW>(filename_tmp,*levelset_array(grid_index));
                        if(((height_array(grid_index)->counts).x < counts_array(grid_index).x) || ((height_array(grid_index)->counts).y < counts_array(grid_index).y)) success = false;}
                    else success=false;}}

            if(success){
                Update_Surface();
                frame_loaded = frame;
                valid = true;}}}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Update_Surface()
{
    for(int grid_index=1;grid_index<=chimera_grid.Size();grid_index++){
        GRID<TV>& grid=chimera_grid(grid_index);
        if(subdivide_surface)
        {
            triangulated_surface_array(grid_index)->mesh.Initialize_Square_Mesh(grid.counts.x,grid.counts.y);
            triangulated_surface_array(grid_index)->particles.array_collection->Clean_Memory();
            triangulated_surface_array(grid_index)->particles.array_collection->Preallocate(grid.counts.x*grid.counts.y);
            for (int i = domain_array(grid_index).min_corner.x; i <= domain_array(grid_index).max_corner.x; i++) for (int j = domain_array(grid_index).min_corner.y; j <= domain_array(grid_index).max_corner.y; j++)
                triangulated_surface_array(grid_index)->particles.array_collection->Add_Element();
        }
        
        for (int i = domain_array(grid_index).min_corner.x; i <= domain_array(grid_index).max_corner.x; i++) for (int j = domain_array(grid_index).min_corner.y; j <= domain_array(grid_index).max_corner.y; j++)
        {
            triangulated_surface_array(grid_index)->particles.X(To_Linear_Index(grid_index,i,j)) = 
                VECTOR<T,3>(grid.Axis_X(i,1), scale*((*height_array(grid_index))(i,j)+vertical_offset), grid.Axis_X(j,2));
        }
        
        if(subdivide_surface) triangulated_surface_array(grid_index)->Loop_Subdivide();
        
        if(opengl_triangulated_surface_array(grid_index)->Is_Smooth_Normals())
            opengl_triangulated_surface_array(grid_index)->Initialize_Vertex_Normals();
        else
            opengl_triangulated_surface_array(grid_index)->Delete_Vertex_Normals();}
}

template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T> *new_selection = 0;
    OPENGL_SELECTION *selection = opengl_triangulated_surface_array(current_grid)->Get_Selection(buffer,buffer_size);
    if(selection)
    {
        if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX)
        {
            int index = ((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T> *)selection)->index;
            new_selection = new OPENGL_SELECTION_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T>(this);
            new_selection->index = From_Linear_Index(current_grid,index);
        }
        delete selection;
    }
    return new_selection;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type != OPENGL_SELECTION::COMPONENT_HEIGHTFIELD_2D) return;
    OPENGL_SELECTION_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T>*)selection;
    int vertex_index = To_Linear_Index(current_grid,real_selection->index.x, real_selection->index.y);
    OPENGL_SELECTION *surface_selection = opengl_triangulated_surface_array(current_grid)->Get_Vertex_Selection(vertex_index);
    opengl_triangulated_surface_array(current_grid)->Highlight_Selection(surface_selection);
    delete surface_selection; // Highlight_Selection made a copy of it
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Clear_Highlight()
{
    opengl_triangulated_surface_array(current_grid)->Clear_Highlight();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Set_Scale(T scale_input)
{
    scale = scale_input;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Increase_Scale()
{
    scale *= (T)1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Decrease_Scale()
{
    scale *= 1/(T)1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Increase_Displacement_Scale()
{
    displacement_scale *= (T)1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Decrease_Displacement_Scale()
{
    displacement_scale *= 1/(T)1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Toggle_Subdivision()
{
    subdivide_surface = !subdivide_surface;
    Reinitialize(true);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Toggle_Draw_All_Grids()
{
    draw_all_grids=!draw_all_grids;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Next_Grid()
{
    current_grid=min(current_grid+1,chimera_grid.m);
    std::stringstream ss;
    ss<<"viewing "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T,RW>::
Previous_Grid()
{
    current_grid=max(current_grid-1,1);
    std::stringstream ss;
    ss<<"viewing "<<current_grid<<std::endl;
    LOG::filecout(ss.str());
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_CHIMERA_HEIGHTFIELD_2D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Centered_Box();
}

template class OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CHIMERA_HEIGHTFIELD_2D<double,double>;
#endif
