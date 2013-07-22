//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D__
#define __OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CHIMERA_SLICE.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T,class RW=T>
class OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef ARRAY<T,FACE_INDEX<3> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::LINEAR_INTERPOLATION_SCALAR::template REBIND<TV>::TYPE T_LINEAR_INTERPOLATION_VECTOR;
    OPENGL_MAC_VELOCITY_FIELD_3D<T>* opengl_mac_velocity_field;
    ARRAY<OPENGL_MAC_VELOCITY_FIELD_3D<T>* > opengl_mac_velocity_fields;
    OPENGL_SCALAR_FIELD_3D<T>* opengl_vorticity_magnitude;
    bool draw_vorticity;
    OPENGL_CHIMERA_SLICE<T>* chimera_slice;
private:
    std::string velocity_filename_set;
    int frame_loaded;
    int current_grid;
    bool valid;
    bool draw_all_grids;
    T min_vorticity,max_vorticity;
    bool is_moving_grid;
    std::string rigid_grid_frame_filename_set;

public:
    OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D(const ARRAY<GRID<TV> > &grid_input,const std::string &velocity_filename_set_input,bool is_moving_grid_input=false,std::string rigid_grid_frame_filename_set_input="");
    
    void Initialize(const ARRAY<GRID<TV> > &grid_array_input);
    ~OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D();

    void Set_Chimera_Slice(OPENGL_CHIMERA_SLICE<T>* slice_input)
    {chimera_slice=slice_input;for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Set_Slice(slice_input->opengl_uniform_slices(i));}

    void Slice_Has_Changed() PHYSBAM_OVERRIDE 
    {if(draw) for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Slice_Has_Changed();opengl_vorticity_magnitude->Slice_Has_Changed();}

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Toggle_Velocity_Mode();
    void Toggle_Velocity_Mode_And_Draw();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Toggle_Draw_Divergence();
    void Toggle_Draw_All_Grids();
    void Next_Grid();
    void Previous_Grid();
    void Toggle_Draw_Vorticity();
    void Normalize_Vorticity_Color_Map();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D,Next_Grid,"Switch to next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D,Previous_Grid,"Switch to previous grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D,Toggle_Draw_All_Grids,"Toggle mutliple/single grid draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D, Toggle_Velocity_Mode, "Toggle velocity mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D, Toggle_Velocity_Mode_And_Draw, "Toggle velocity mode and draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D, Toggle_Arrowhead, "Toggle arrowhead");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D, Toggle_Draw_Vorticity, "Toggle draw vorticity");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_3D, Normalize_Vorticity_Color_Map, "Normalize vorticity map based on current frame");

private:
    void Reinitialize();
    void Update_Divergence();
    void Update_Streamlines();
    void Update_Vorticity();
};

}

#endif
