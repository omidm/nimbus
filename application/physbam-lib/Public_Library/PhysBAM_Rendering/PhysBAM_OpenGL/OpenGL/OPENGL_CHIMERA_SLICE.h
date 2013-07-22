//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_CHIMERA_SLICE
//##################################################################### 
#ifndef __OPENGL_CHIMERA_SLICE__
#define __OPENGL_CHIMERA_SLICE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>

namespace PhysBAM{
class OPENGL_WORLD;
template<class T>
class OPENGL_CHIMERA_SLICE:public OPENGL_SLICE
{
    typedef GRID<VECTOR<T,3> > T_GRID;
public:
    OPENGL_WORLD &world;
    ARRAY<OPENGL_UNIFORM_SLICE*> opengl_uniform_slices;
    OPENGL_SLICE::SLICE_MODE mode;
    int current_grid;
    ARRAY<PLANE<T> > plane1_array;
    ARRAY<PLANE<T> > plane2_array;

    OPENGL_CHIMERA_SLICE(OPENGL_WORLD &world_input)
        : world(world_input), current_grid(1)
    {
    }

    bool Is_Slice_Mode() PHYSBAM_OVERRIDE
    {
        return mode==NODE_SLICE || mode==CELL_SLICE;
    }

    void Initialize(ARRAY<GRID<VECTOR<T,3> > > grids_input)
    {
        mode=NO_SLICE;
        for(int i=1;i<=grids_input.Size();i++){
            opengl_uniform_slices.Append(new OPENGL_UNIFORM_SLICE(world));
            opengl_uniform_slices(opengl_uniform_slices.m)->Initialize(grids_input(i));}
        plane1_array.Resize(grids_input.Size());
        plane2_array.Resize(grids_input.Size());
        Update_Chimera_Clip_Planes();
    }

    void Set_Slice_Mode(OPENGL_SLICE::SLICE_MODE mode_input) PHYSBAM_OVERRIDE
    {
        mode=mode_input;
        for(int i=1;i<=opengl_uniform_slices.Size();i++){
            opengl_uniform_slices(i)->Set_Slice_Mode(mode);}
        Update_Chimera_Clip_Planes();
    }

    void Toggle_Slice_Mode() PHYSBAM_OVERRIDE
    {
        Set_Slice_Mode((OPENGL_SLICE::SLICE_MODE)(((int)mode+1)%3));
        Update_Chimera_Clip_Planes();
    }
    
    void Toggle_Slice_Axis() PHYSBAM_OVERRIDE
    {
        for(int i=1;i<=opengl_uniform_slices.Size();i++){
            opengl_uniform_slices(i)->Toggle_Slice_Axis();}
        Update_Chimera_Clip_Planes();
    }

    void Increment_Slice() PHYSBAM_OVERRIDE
    {
        opengl_uniform_slices(current_grid)->Increment_Slice();
        Update_Chimera_Clip_Planes();
    }

    void Decrement_Slice() PHYSBAM_OVERRIDE
    {
        opengl_uniform_slices(current_grid)->Decrement_Slice();
        Update_Chimera_Clip_Planes();
    }
    
    void Enable_Clip_Planes() PHYSBAM_OVERRIDE
    {
        for(int i=1;i<=opengl_uniform_slices.Size();i++){
            opengl_uniform_slices(i)->Enable_Clip_Planes();}
    }

    void Update_Chimera_Clip_Planes()
    {
        for(int grid_index=1;grid_index<=opengl_uniform_slices.m;grid_index++){
            T_GRID grid=opengl_uniform_slices(grid_index)->grid;
            int index=opengl_uniform_slices(grid_index)->index;
            int axis=opengl_uniform_slices(grid_index)->axis;
            T pos=(mode==NODE_SLICE)?grid.Node(index,index,index)[axis]:grid.Center(index,index,index)[axis];
            PLANE<T> plane1(VECTOR<T,3>(0,0,0),VECTOR<T,3>(0,0,0)),plane2(VECTOR<T,3>(0,0,0),VECTOR<T,3>(0,0,0));
            plane1.normal[axis]=-1;plane1.x1[axis]=pos-grid.dX[axis]/(T)1.9;
            plane2.normal[axis]=1;plane2.x1[axis]=pos+grid.dX[axis]/(T)1.9;
            plane1_array(grid_index)=plane1;plane2_array(grid_index)=plane2;}
    }

    void Print_Slice_Info(std::ostream& output_stream) PHYSBAM_OVERRIDE;
private:
    void Update_Clip_Planes() PHYSBAM_OVERRIDE 
    {}
};
}
#endif
