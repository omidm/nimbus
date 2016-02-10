//#####################################################################
// Copyright 2007, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXACT_SIMPLEX_INTERACTIONS
//##################################################################### 
#ifndef __EXACT_SIMPLEX_INTERACTIONS__
#define __EXACT_SIMPLEX_INTERACTIONS__    

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class EXACT_SIMPLEX_INTERACTIONS
{
public:
//#####################################################################
    static bool Exact_Segment_Intersection(const VECTOR<float,2>& p1,const VECTOR<float,2>& p2,const VECTOR<float,2>& q3,const VECTOR<float,2>& q4,
        const int rank1,const int rank2,const int rank3,const int rank4);
    static bool Exact_Segment_Intersection(const VECTOR<float,2>& p1,const VECTOR<float,2>& p2,const VECTOR<float,2>& q3,const VECTOR<float,2>& q4);
    static bool Positive_Signed_Area(VECTOR<float,2> x1,VECTOR<float,2> x2,VECTOR<float,2> x3,int rank1,int rank2,int rank3);
    static bool Positive_Signed_Area(VECTOR<float,2> x1,VECTOR<float,2> x2,VECTOR<float,2> x3);
    static double Exact_Signed_Area(const VECTOR<float,2>& x1,const VECTOR<float,2>& x2,const VECTOR<float,2>& x3);
//#####################################################################
};   
}
#endif

