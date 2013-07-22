//#####################################################################
// Copyright 2005-2007, Kevin Der, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLEX_INTERACTIONS
//##################################################################### 
#ifndef __SIMPLEX_INTERACTIONS__
#define __SIMPLEX_INTERACTIONS__    

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T>
class SIMPLEX_INTERACTIONS
{
public:
//#####################################################################
    static bool Intersection(const VECTOR<VECTOR<T,2>,2>& segment_1,const VECTOR<VECTOR<T,2>,2>& segment_2);
    static bool Two_Segment_Intersection_Barycentric_Coordinates(const VECTOR<VECTOR<T,2>,2>& segment1,const VECTOR<VECTOR<T,2>,2>& segment2,VECTOR<T,1>& weights1,VECTOR<T,1>& weights2);
    static bool Intersection(const VECTOR<VECTOR<T,2>,3>& triangle,const VECTOR<T,2>& point);
    static bool Intersection(const VECTOR<VECTOR<T,2>,3>& triangle,const VECTOR<VECTOR<T,2>,2>& segment);
    static bool Intersection(const VECTOR<VECTOR<T,2>,3>& triangle_1,const VECTOR<VECTOR<T,2>,3>& triangle_2);
    static bool Intersection(const VECTOR<VECTOR<T,3>,3>& triangle_1,const VECTOR<VECTOR<T,3>,3>& triangle_2);
    static bool Intersection(const VECTOR<VECTOR<T,3>,3>& triangle,const VECTOR<VECTOR<T,3>,2>& segment);
    static T Intersection_Fraction(const VECTOR<VECTOR<T,3>,3>& triangle,const VECTOR<VECTOR<T,3>,2>& segment);
    static bool Intersection(const VECTOR<VECTOR<T,3>,4>& tetrahedron,const VECTOR<VECTOR<T,3>,2>& segment);
    static bool Intersection(const VECTOR<VECTOR<T,3>,4>& tetrahedron,const VECTOR<VECTOR<T,3>,3>& triangle);
    static bool Intersection(const VECTOR<VECTOR<T,3>,4>& tetrahedron_1,const VECTOR<VECTOR<T,3>,4>& tetrahedron_2);
    static bool Three_Triangle_Intersection_Barycentric_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
        VECTOR<T,2>& weights1,VECTOR<T,2>& weights2,VECTOR<T,2>& weights3);
    static void Triangle_Segment_Intersection_Barycentric_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle,const VECTOR<VECTOR<T,3>,2>& segment,VECTOR<T,2>& triangle_weights,T& segment_weight);
//#####################################################################
};   
}
#endif

