//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D.h>
using namespace PhysBAM;
    
template<class T,class RW> OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D(const ARRAY<GRID<TV> > &grid_input,const std::string &velocity_filename_set_input,bool is_moving_grid_input,std::string rigid_grid_frame_filename_set_input)
    : OPENGL_COMPONENT("MAC Velocity Field 3D"),draw_vorticity(false),
     velocity_filename_set(velocity_filename_set_input),current_grid(1),valid(false),draw_all_grids(true),min_vorticity(FLT_MAX),max_vorticity(FLT_MIN),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input)
{
    Initialize(grid_input);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Initialize(const ARRAY<GRID<TV> > &grid_input)
{
    is_animation=true;
    frame_loaded=-1;

    int number_of_grids=grid_input.Size();

    opengl_mac_velocity_fields.Resize(number_of_grids);
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) 
        opengl_mac_velocity_fields(i)=new OPENGL_MAC_VELOCITY_FIELD_3D<T>(*(new GRID<TV>(grid_input(i))),is_moving_grid);
    opengl_mac_velocity_field=opengl_mac_velocity_fields(1);
    
    opengl_vorticity_magnitude=new OPENGL_SCALAR_FIELD_3D<T>(opengl_mac_velocity_field->grid,*(new ARRAY<T,VECTOR<int,3> >),OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1),OPENGL_SCALAR_FIELD_3D<T>::DRAW_TEXTURE,is_moving_grid);
    opengl_vorticity_magnitude->Set_Scale_Range(0,100);

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
~OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++)
        delete &opengl_mac_velocity_fields(i)->grid;
    delete opengl_vorticity_magnitude;
}

template<class T,class RW> bool OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(velocity_filename_set.c_str(),frame_input,1));
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": ";
        opengl_mac_velocity_field->Print_Selection_Info(stream,selection);
        if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_3D && opengl_mac_velocity_field->grid.Is_MAC_Grid()){
            VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)selection)->index;
            if(draw_vorticity) stream<<"vorticity magnitude = "<<opengl_vorticity_magnitude->values(index)<<std::endl;}}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Display(const int in_color) const
{
    if(valid){
        if(draw){
            if(draw_all_grids) 
                for(int i=1;i<=opengl_mac_velocity_fields.m;i++) 
                    opengl_mac_velocity_fields(i)->Display(in_color);
            else 
                opengl_mac_velocity_fields(current_grid)->Display(in_color);}
        if(draw_vorticity) opengl_vorticity_magnitude->Display(in_color);}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_mac_velocity_field->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Reinitialize()
{
    if (draw){
        if ((is_animation && (frame_loaded!=frame)) || (!is_animation && frame_loaded < 0)){
            valid = false;
            if(is_moving_grid){
                for(int i=1;i<=opengl_mac_velocity_fields.m;i++){
                    std::string grid_filename=STRING_UTILITIES::string_sprintf(rigid_grid_frame_filename_set.c_str(),frame,i);
                    if(FILE_UTILITIES::File_Exists(grid_filename)) FILE_UTILITIES::Read_From_File<T>(grid_filename,opengl_mac_velocity_fields(i)->rigid_grid_frame);}}
            for(int i=1;i<=opengl_mac_velocity_fields.m;i++){
                std::string tmp_filename=STRING_UTILITIES::string_sprintf(velocity_filename_set.c_str(),frame,i);
                if(FILE_UTILITIES::File_Exists(tmp_filename)){
                    FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_mac_velocity_fields(i)->face_velocities);}
                else return;
                opengl_mac_velocity_fields(i)->Update();}
            frame_loaded=frame;
            valid=true;
            if(draw_vorticity){Update_Vorticity();opengl_vorticity_magnitude->Update();}}}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Velocity_Mode()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Toggle_Velocity_Mode();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Velocity_Mode_And_Draw()
{
    if (draw)
    {
        Toggle_Velocity_Mode();
        if ((int)opengl_mac_velocity_field->velocity_mode==0) Toggle_Draw();
    }
    else Toggle_Draw();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Increase_Vector_Size()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Scale_Vector_Size((T)1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Decrease_Vector_Size()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Scale_Vector_Size(1/(T)1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Arrowhead()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->draw_arrowhead = !opengl_mac_velocity_fields(i)->draw_arrowhead;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Draw_All_Grids()
{
    draw_all_grids=!draw_all_grids;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Next_Grid(){
    current_grid=min(current_grid+1,opengl_mac_velocity_fields.m);
    opengl_mac_velocity_field=opengl_mac_velocity_fields(current_grid);
    opengl_vorticity_magnitude->Set_Slice(chimera_slice->opengl_uniform_slices(current_grid));
    if(draw_vorticity){Update_Vorticity();opengl_vorticity_magnitude->Update();}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Previous_Grid(){
    current_grid=max(current_grid-1,1);
    opengl_mac_velocity_field=opengl_mac_velocity_fields(current_grid);
    opengl_vorticity_magnitude->Set_Slice(chimera_slice->opengl_uniform_slices(current_grid));
    if(draw_vorticity){Update_Vorticity();opengl_vorticity_magnitude->Update();}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Draw_Vorticity()
{
    draw_vorticity=!draw_vorticity;
    if(draw_vorticity) valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Normalize_Vorticity_Color_Map()
{
    if(!draw_vorticity) return;
    opengl_vorticity_magnitude->Set_Scale_Range(min_vorticity,max_vorticity);
    valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Update_Vorticity()
{
    GRID<TV>& grid=opengl_mac_velocity_field->grid;
    RANGE<VECTOR<int,3> > domain_indices(grid.Domain_Indices());domain_indices.Change_Size(-VECTOR<int,3>::All_Ones_Vector());
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup(opengl_mac_velocity_field->face_velocities);
    opengl_vorticity_magnitude->values.Resize(grid.Domain_Indices());
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,domain_indices);iterator.Valid();iterator.Next()){VECTOR<int,3> index=iterator.Cell_Index();
        T vorticity_magnitude=VORTICITY_UNIFORM<TV>::Vorticity(grid,lookup,index).Magnitude();
        opengl_vorticity_magnitude->values(index)=vorticity_magnitude;
        min_vorticity=min(min_vorticity,vorticity_magnitude);
        max_vorticity=max(max_vorticity,vorticity_magnitude);}
}

template<class T, class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Draw_Divergence()
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}

template<class T, class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Update_Divergence()
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}

template<class T, class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<T,RW>::
Update_Streamlines()
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}

template class OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D<double,double>;
#endif
