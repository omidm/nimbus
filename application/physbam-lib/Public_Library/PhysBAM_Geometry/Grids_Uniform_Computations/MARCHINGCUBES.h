//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHINGCUBES
//#####################################################################
#ifndef __MARCHINGCUBES__
#define __MARCHINGCUBES__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
namespace PhysBAM
{
template<class T> class TRIANGULATED_SURFACE;

template<class T>
class MARCHINGCUBES:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV; typedef VECTOR<int,3> TV_INT;
    typedef GRID<TV> T_GRID;
private:
    LEVELSET_3D<T_GRID>& levelset;
    T contour_value;
    T_GRID grid;
    bool is_distance_field;

    T_GRID edge_grid[3];
    ARRAY<int,TV_INT> edge_grid_array[3];

    ARRAY<TV> vertices;
    ARRAY<TV_INT> triangles;
    ARRAY<TV> normals;  ////vertex normals

public:
    MARCHINGCUBES(LEVELSET_3D<T_GRID>& levelset_input,const T contour_value_input=0,const bool is_distance_field_input=false)
        :levelset(levelset_input),contour_value(contour_value_input),is_distance_field(is_distance_field_input)
    {
        if(levelset_input.grid.MAC_offset==(T).5)
            grid=levelset_input.grid.Get_Regular_Grid_At_MAC_Positions();
        else grid=levelset_input.grid;
        Initialize_Marching_Cubes();
    }

    void Initialize_Marching_Cubes();
    void Marching_Cubes();
    static TRIANGULATED_SURFACE<T>* Create_Triangulated_Surface_From_Levelset(LEVELSET_3D<GRID<TV> >& levelset,const T contour_value_input=0,const bool is_distance_field_input=false)
    {MARCHINGCUBES<T> mc(levelset,contour_value_input,is_distance_field_input);mc.Marching_Cubes();return mc.Get_Triangulated_Surface();}
    TRIANGULATED_SURFACE<T>* Get_Triangulated_Surface();
    void Get_Triangulated_Surface(TRIANGULATED_SURFACE<T>& surface);

    ////common interfaces of MARCHINGCUBES,MARCHINGSQUARES,MARCHINGTRIANGLES
    static TRIANGULATED_SURFACE<T>* Create_Surface_From_Levelset(LEVELSET_3D<T_GRID>& levelset,const T contour_value_input=0)
    {return Create_Triangulated_Surface_From_Levelset(levelset,contour_value_input);}
    TRIANGULATED_SURFACE<T>* Get_Surface(){return Get_Triangulated_Surface();}
    void Get_Surface(TRIANGULATED_SURFACE<T>& surface){Get_Triangulated_Surface(surface);}

private:
    int Get_Particle_Index_On_Edge(const TV_INT index,const int edge_index);
    void Set_Particle_Index_On_Edge(TV_INT cell_index,int edge_index,int value_input);
    unsigned int Get_Cell_Type(T v0,T v1,T v2,T v3,T v4,T v5,T v6,T v7,T iso_value);
    bool Has_Particle_On_Edge(TV_INT cell_index,int edge_index/*0-11*/);
    void Get_Edge_Vertex_Index(TV_INT cell_index,int edge_index,/*result*/TV_INT& v1,/*result*/TV_INT& v2);
};
}
#endif
