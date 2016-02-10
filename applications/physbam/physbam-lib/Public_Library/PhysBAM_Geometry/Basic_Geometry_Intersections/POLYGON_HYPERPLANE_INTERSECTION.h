//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __POLYGON_HYPERPLANE_INTERSECTION__
#define __POLYGON_HYPERPLANE_INTERSECTION__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Cut_Convex_Polygon_With_Hyperplane_And_Discard_Outside_Polygon(POINT_SIMPLEX_1D<T>& polygon,const POINT_SIMPLEX_1D<T>& hyperplane)
{
    if(hyperplane.Signed_Distance(polygon.x1)>0)
        return false;
    return true;
}
template<class T> bool Cut_Convex_Polygon_With_Hyperplane_And_Discard_Outside_Polygon(SEGMENT_2D<T>& polygon,const LINE_2D<T>& hyperplane)
{
    ARRAY<SEGMENT_2D<T> > inside_polygons;
    SEGMENT_2D<T>::Cut_With_Hyperplane_And_Discard_Outside_Simplices(polygon,hyperplane,inside_polygons);
    if(!inside_polygons.Size())
        return false;
    polygon=inside_polygons(1);
    return true;
}
template<class T> bool Cut_Convex_Polygon_With_Hyperplane_And_Discard_Outside_Polygon(POLYGON<VECTOR<T,3> >& polygon,const PLANE<T>& hyperplane)
{
    typedef VECTOR<T,3> TV;
    //const T interpolation_tolerance=1e-2;
    //LOG::cout << "clipping polygon " << polygon.X << " " << hyperplane.x1 << " " << hyperplane.normal << std::endl;

    int n=polygon.X.Size();

    //T maximum_length=0;
    int minimum_index=0;
    T minimum_distance=0;
    int maximum_index=0;
    T maximum_distance=0;
    for(int index=1;index<=n;index++){
        T distance=hyperplane.Signed_Distance(polygon.X(index));
        if(distance<minimum_distance){
            minimum_distance=distance;
            minimum_index=index;}
        if(distance>maximum_distance){
            maximum_distance=distance;
            maximum_index=index;}}
    
    if(!minimum_index) //ALL OUTSIDE
        return false;
    if(!maximum_index) //ALL INSIDE
        return true;
    
    int lower_index_1=0,upper_index_1=0;
    for(int index=maximum_index;index!=minimum_index;index=(index%n)+1)
        if(hyperplane.Signed_Distance(polygon.X(index))>=0)
            lower_index_1=index;
    for(int index=minimum_index;index!=maximum_index;index=(index%n)+1)
        if(hyperplane.Signed_Distance(polygon.X(index))<0)
            upper_index_1=index;

    int lower_index_2=lower_index_1%n+1;
    T lower_distance_1=hyperplane.Signed_Distance(polygon.X(lower_index_1));
    T lower_distance_2=hyperplane.Signed_Distance(polygon.X(lower_index_2)); //THIS IS NEGATIVE
    T lower_interpolation_fraction=clamp(lower_distance_1/(lower_distance_1-lower_distance_2),(T)0,(T)1);
    TV lower_location=lower_interpolation_fraction*polygon.X(lower_index_2)+(1-lower_interpolation_fraction)*polygon.X(lower_index_1);

    int upper_index_2=upper_index_1%n+1;
    T upper_distance_1=hyperplane.Signed_Distance(polygon.X(upper_index_1)); //THIS IS NEGATIVE
    T upper_distance_2=hyperplane.Signed_Distance(polygon.X(upper_index_2));
    T upper_interpolation_fraction=clamp(upper_distance_2/(upper_distance_2-upper_distance_1),(T)0,(T)1);
    TV upper_location=upper_interpolation_fraction*polygon.X(upper_index_1)+(1-upper_interpolation_fraction)*polygon.X(upper_index_2);

    //LOG::cout << "lower " << lower_index_1 << " " << lower_index_2 << " " << lower_interpolation_fraction << " - " << lower_distance_1 << " " << lower_distance_2 << std::endl;
    //LOG::cout << "upper " << upper_index_1 << " " << upper_index_2 << " " << upper_interpolation_fraction << " - " << upper_distance_1 << " " << upper_distance_2 << std::endl;

    //EETODO: REMOVE REDUNDANT VERTICES

    if(upper_index_2<=lower_index_2){
        for(int index=lower_index_1;index>=upper_index_2;index--)
            polygon.X.Remove_Index(index);
        polygon.X.Insert(lower_location,upper_index_2);
        polygon.X.Insert(upper_location,upper_index_2);
    }else{
        for(int index=n;index>=upper_index_2;index--)
            polygon.X.Remove_Index(index);
        for(int index=lower_index_2-1;index>=1;index--)
            polygon.X.Remove_Index(index);
        polygon.X.Append(upper_location);
        polygon.X.Append(lower_location);}

    //LOG::cout << "clipped polygon " << inside_polygon.X << std::endl;
    
    return true;
}
//#####################################################################
};
};
#endif
