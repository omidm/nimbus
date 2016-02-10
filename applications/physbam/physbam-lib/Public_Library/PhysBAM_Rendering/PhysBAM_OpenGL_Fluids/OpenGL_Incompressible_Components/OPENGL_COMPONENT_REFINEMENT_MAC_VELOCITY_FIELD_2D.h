//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D__
#define __OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T,class RW=T>
class OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
public:
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef ARRAY<T,FACE_INDEX<2> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    GRID<TV> grid;
    ARRAY<GRID<TV>,TV_INT> sub_grids;
    ARRAY<OPENGL_MAC_VELOCITY_FIELD_2D<T>*,TV_INT> opengl_refinement_mac_velocity_fields;
private:
    std::string velocity_filename,sub_grids_filename;
    int frame_loaded;
    int level,number_of_levels;
    int level_loaded;
    bool valid;
    bool draw_all_levels;
    ARRAY<ARRAY<T,FACE_INDEX<2> >,TV_INT> sub_mac_velocities;
    T size;
    OPENGL_COLOR vector_color;
    typename OPENGL_MAC_VELOCITY_FIELD_2D<T>::VELOCITY_MODE velocity_mode;

public:
    OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D(const GRID<TV> &grid,const std::string &velocity_filename_input,const std::string &sub_grids_filename_input);
    void Initialize();
    ~OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Set_Size(T s)
    {size=s;for(int i=1;i<=opengl_refinement_mac_velocity_fields.array.m;i++) if(opengl_refinement_mac_velocity_fields.array(i)) opengl_refinement_mac_velocity_fields.array(i)->size=s;}
    void Set_Vector_Color(OPENGL_COLOR c)
    {vector_color=c;for(int i=1;i<=opengl_refinement_mac_velocity_fields.array.m;i++) if(opengl_refinement_mac_velocity_fields.array(i)) opengl_refinement_mac_velocity_fields.array(i)->vector_color=c;}
    void Set_Velocity_Mode(typename OPENGL_MAC_VELOCITY_FIELD_2D<T>::VELOCITY_MODE m)
    {velocity_mode=m;for(int i=1;i<=opengl_refinement_mac_velocity_fields.array.m;i++) if(opengl_refinement_mac_velocity_fields.array(i)) opengl_refinement_mac_velocity_fields.array(i)->Set_Velocity_Mode(m);}

    void Toggle_Velocity_Mode();
    void Toggle_Velocity_Mode_And_Draw();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Toggle_Draw_All_Levels();
    void Next_Level();
    void Previous_Level();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D, Toggle_Velocity_Mode, "Toggle velocity mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D, Toggle_Velocity_Mode_And_Draw, "Toggle velocity mode and draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_REFINEMENT_MAC_VELOCITY_FIELD_2D, Toggle_Arrowhead, "Toggle arrowhead");

private:
    void Reinitialize();
};

}

#endif
