//#####################################################################
// Copyright 2004-2009, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_CHIMERA_GRID_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_CHIMERA_GRID_3D__
#define __OPENGL_COMPONENT_CHIMERA_GRID_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CHIMERA_SLICE.h>

namespace PhysBAM
{
template<class T,class RW=T>
class OPENGL_COMPONENT_CHIMERA_GRID_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
public:
    ARRAY<OPENGL_GRID_3D<T>* > opengl_grids;
    OPENGL_GRID_3D<T>* opengl_grid;
    OPENGL_CHIMERA_SLICE<T>* chimera_slice;
private:
    std::string filename_set;
    int frame_loaded;
    bool valid;
    int current_grid;
    bool draw_all_grids;
    bool is_moving_grid;
    std::string rigid_grid_frame_filename_set;
public:

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size)
    {return opengl_grid->Get_Selection(buffer,buffer_size);}

    virtual void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE
    {opengl_grid->Highlight_Selection(selection);}

    virtual void Clear_Highlight() PHYSBAM_OVERRIDE
    {opengl_grid->Clear_Highlight();}

    void Set_Chimera_Slice(OPENGL_CHIMERA_SLICE<T>* slice_input)
    {chimera_slice=slice_input;for(int i=1;i<=opengl_grids.m;i++) opengl_grids(i)->Set_Slice(slice_input->opengl_uniform_slices(i));}

    void Slice_Has_Changed() PHYSBAM_OVERRIDE 
    {for(int i=1;i<=opengl_grids.m;i++) opengl_grids(i)->Slice_Has_Changed();}

//#####################################################################
    OPENGL_COMPONENT_CHIMERA_GRID_3D(const std::string filename_set_input,bool is_moving_grid_input=false,std::string rigid_grid_frame_filename_set_input="");
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw && valid;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    void Toggle_Draw_All_Grids();
    void Toggle_Draw_Ghost_Values();
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_GRID_3D,Next_Grid,"Switch to next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_GRID_3D,Previous_Grid,"Switch to previous grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_GRID_3D,Toggle_Draw_All_Grids,"Toggle to draw all the grids or not");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CHIMERA_GRID_3D,Toggle_Draw_Ghost_Values,"Toggle to draw ghost cell values");

private:
    void Reinitialize();
    void Next_Grid();
    void Previous_Grid();
//#####################################################################
};
}

#endif
