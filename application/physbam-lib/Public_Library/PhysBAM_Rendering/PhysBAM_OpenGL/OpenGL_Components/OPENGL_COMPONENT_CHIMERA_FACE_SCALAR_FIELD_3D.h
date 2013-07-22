//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D__
#define __OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D__
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CHIMERA_SLICE.h>
#include <string>

namespace PhysBAM{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
public:
    OPENGL_FACE_SCALAR_FIELD_3D<T,T2>*  opengl_scalar_field;
    ARRAY<OPENGL_FACE_SCALAR_FIELD_3D<T,T2>* > opengl_scalar_fields;
    OPENGL_FACE_SCALAR_FIELD_3D<T,bool>*  opengl_bool_field;
    ARRAY<OPENGL_FACE_SCALAR_FIELD_3D<T,bool>* > opengl_bool_fields;
    OPENGL_CHIMERA_SLICE<T>* chimera_slice;
private:
    std::string values_filename_set;
    int frame_loaded;
    bool valid;
    int current_grid;
    bool draw_all_grids;
    bool display_in_bool;
    bool is_moving_grid;
    std::string rigid_grid_frame_filename_set;

public:
    OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D(const ARRAY<GRID<TV> > &grid_input,const std::string &values_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool is_moving_grid_input=false,std::string rigid_grid_frame_filename_set_input="");
    ~OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D();

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE
    {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw && valid;}

    void Set_Chimera_Slice(OPENGL_CHIMERA_SLICE<T>* slice_input)
    {chimera_slice=slice_input;for(int i=1;i<=opengl_scalar_fields.m;i++){opengl_scalar_fields(i)->Set_Slice(slice_input->opengl_uniform_slices(i));opengl_bool_fields(i)->Set_Slice(slice_input->opengl_uniform_slices(i));}}
    
    void Slice_Has_Changed() PHYSBAM_OVERRIDE 
    {if(draw) for(int i=1;i<=opengl_scalar_fields.m;i++){opengl_scalar_fields(i)->Slice_Has_Changed();opengl_bool_fields(i)->Slice_Has_Changed();}}

//#####################################################################
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Display_In_Bool(const int in_color=1) const;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Toggle_Draw_All_Grids();
    void Set_To_Display_In_Bool_Mode();
    
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D,Next_Grid,"Switch to next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D,Previous_Grid,"Switch to previous grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_FACE_SCALAR_FIELD_3D,Toggle_Draw_All_Grids,"Toggle draw all grids");
private:
    void Reinitialize();
    void Next_Grid();
    void Previous_Grid();
//#####################################################################
};
}
#endif
