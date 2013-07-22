//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHINGSQUARES
//#####################################################################
#ifndef __MARCHINGSQUARES__
#define __MARCHINGSQUARES__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
namespace PhysBAM{

template<class T> class SEGMENTED_CURVE_2D;

template<class T>
class MARCHINGSQUARES:public NONCOPYABLE
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef GRID<TV> T_GRID;
private:
    LEVELSET_2D<T_GRID>& levelset;
    T contour_value;
    GRID<TV> grid;
    bool is_distance_field;

    ARRAY<TV> vertices;
    ARRAY<TV_INT> edges;

public:
    MARCHINGSQUARES(LEVELSET_2D<T_GRID>& levelset_input,const T contour_value_input=0,const bool is_distance_field_input=false)
        :levelset(levelset_input),contour_value(contour_value_input),is_distance_field(is_distance_field_input)
    {
        if(levelset_input.grid.MAC_offset==(T).5)
            grid=levelset_input.grid.Get_Regular_Grid_At_MAC_Positions();
        else grid=levelset_input.grid;
    }

    void Marching_Squares();
    static SEGMENTED_CURVE_2D<T>* Create_Segmented_Curve_From_Levelset(LEVELSET_2D<T_GRID>& levelset,const T contour_value_input=0,const bool is_distance_field_input=true)
    {
        MARCHINGSQUARES<T> ms(levelset,contour_value_input,is_distance_field_input);
        ms.Marching_Squares();
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
    unsigned int Square_Type(T v0,T v1,T v2,T v3,T iso_value=0);
};
}

#endif
