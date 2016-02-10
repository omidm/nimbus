//#####################################################################
// Copyright 2006, Kevin Der.
// Copyright 2008, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGONAL_TRIANGULATION
//##################################################################### 
#ifndef __POLYGONAL_TRIANGULATION__
#define __POLYGONAL_TRIANGULATION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Polygon_Orientation.h>
namespace PhysBAM{

template<class T>
class POLYGONAL_TRIANGULATION
{
private:
    POLYGONAL_TRIANGULATION()
    {}
    
    ~POLYGONAL_TRIANGULATION()
    {}

//#####################################################################
public:
    template<class VECTORT2_ARRAY,class INT_ARRAY>
    static int Triangulate_Nonconvex_Simple_Polygon(const VECTORT2_ARRAY& coordinates,const INT_ARRAY& polygon_input,ARRAY< VECTOR<int,3> >& triangles,bool keep_degenerate_triangles=false);
    template<class VECTORT2_ARRAY,class INT_ARRAY_ARRAY>
    static int Triangulate_Nonconvex_Nonsimple_Polygon(const VECTORT2_ARRAY& coordinates,const INT_ARRAY_ARRAY& polygon_input,ARRAY< VECTOR<int,3> >& triangles,bool keep_degenerate_triangles=false);
    static void Triangulate_Convex_Planar_Polygon(const ARRAY<VECTOR<T,2> >& positions,ARRAY<VECTOR<int,2> > segments,ARRAY<VECTOR<int,3> >& triangles);
    static void Triangulate_Nonconvex_Planar_Connected_Polygon(const ARRAY<VECTOR<T,2> >& positions,const ARRAY<ARRAY<VECTOR<int,2> > >& segments_input,ARRAY<VECTOR<int,3> >& triangles);
private:
    static bool Segment_Intersects(const VECTOR<int,2>& candidate_segment,const ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<VECTOR<T,2> >& positions);
    static T Find_Sharpest_Angle(const ARRAY<VECTOR<T,2> >& positions,const ARRAY<VECTOR<int,2> >& all_segments,const VECTOR<int,2>& given_segment,const ARRAY<int>& test_segments);
    static int Orientation(const VECTOR<T,2>& x1,const VECTOR<T,2>& x2,const VECTOR<T,2>& x3);
    static bool Point_Triangle_Intersects(const VECTOR<T,2>& y,const VECTOR<T,2>& x1,const VECTOR<T,2>& x2,const VECTOR<T,2>& x3);
    static bool Segment_Segment_Intersects(const VECTOR<T,2>& x1,const VECTOR<T,2>& x2,const VECTOR<T,2>& y1,const VECTOR<T,2>& y2);
//#####################################################################
};
}
#endif
