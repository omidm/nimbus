//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHINGTRIANGLES
//#####################################################################
#ifndef __MARCHINGTRIANGLES__
#define __MARCHINGTRIANGLES__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>

namespace PhysBAM
{
template<class T> class SEGMENTED_CURVE_2D;

template<class T>
class MARCHINGTRIANGLES:public NONCOPYABLE
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef GRID<TV> T_GRID;
private:
    LEVELSET_2D<T_GRID>& levelset;
    T contour_value;
    GRID<TV> grid;
    bool is_distance_field;

    ARRAY<TV> vertices;
    ARRAY<TV_INT> edges;

    T_GRID edge_grid_x,edge_grid_y,edge_grid_c;
    ARRAY<int,TV_INT> edge_grid_array_x,edge_grid_array_y,edge_grid_array_c;
public:
    MARCHINGTRIANGLES(LEVELSET_2D<T_GRID>& levelset_input,const T contour_value_input=0,const bool is_distance_field_input=false)
        :levelset(levelset_input),contour_value(contour_value_input),is_distance_field(is_distance_field_input)
    {
        if(levelset_input.grid.MAC_offset==(T).5)
            grid=levelset_input.grid.Get_Regular_Grid_At_MAC_Positions();
        else grid=levelset_input.grid;
        Initialize_Marching_Triangles();
    }
    void Initialize_Marching_Triangles();
    void Marching_Triangles();
    static SEGMENTED_CURVE_2D<T>* Create_Segmented_Curve_From_Levelset(LEVELSET_2D<T_GRID>& levelset,const T contour_value_input=0,const bool is_distance_field_input=true)
    {
        MARCHINGTRIANGLES<T> ms(levelset,contour_value_input,is_distance_field_input);
        ms.Marching_Triangles();
        return ms.Get_Segmented_Curve();
    }
    SEGMENTED_CURVE_2D<T>* Get_Segmented_Curve();
    void Get_Segmented_Curve(SEGMENTED_CURVE_2D<T>& curve);

    ////common interfaces of MARCHINGCUBES,MARCHINGSQUARES,MARCHINGTRIANGLES
    static SEGMENTED_CURVE_2D<T>* Create_Surface_From_Levelset(LEVELSET_2D<T_GRID>& levelset,const T contour_value_input=0)
    {return Create_Segmented_Curve_From_Levelset(levelset,contour_value_input);}
    SEGMENTED_CURVE_2D<T>* Get_Surface(){return Get_Segmented_Curve();}
    void Get_Surface(SEGMENTED_CURVE_2D<T>& curve){Get_Segmented_Curve(curve);}
private:
    unsigned int Get_Triangle_Type(T t0,T t1,T t2,T iso_value)
    {unsigned int type=0;if(t0<iso_value)type|=1;if(t1<iso_value)type|=2;if(t2<iso_value)type|=4;return type;}
    void Process_Triangle(TV_INT v0,TV_INT v1,TV_INT v2);
    int Get_Particle_Index_On_Edge(TV_INT v0,TV_INT v1);
    void Set_Particle_Index_On_Edge(TV_INT v0,TV_INT v1,int input_value);
    int Add_Particle_To_Edge(TV_INT v0,TV_INT v1);
    bool Add_Edge(TV_INT v0,TV_INT v1){PHYSBAM_NOT_IMPLEMENTED();};
};
}
#endif
