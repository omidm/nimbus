//#####################################################################
// Copyright 2006-2007, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROBUST_SIMPLEX_INTERACTIONS
//##################################################################### 
#ifndef __ROBUST_SIMPLEX_INTERACTIONS__
#define __ROBUST_SIMPLEX_INTERACTIONS__    

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <limits>
namespace PhysBAM{

template<class TV>
class ROBUST_SIMPLEX_INTERACTIONS
{
    // DUMMY CLASS Everything Real is specialized
};

template<class T>
class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    HASHTABLE<VECTOR<TV,3>,VECTOR<T,2> > triple_products; // sum of positive and negative terms for each triple of points
    HASHTABLE<VECTOR<TV,4>,PAIR<T,bool> > signed_volumes_times_six; // volume and indication whether the result is robust or not
    T tolerance;
public:

    ROBUST_SIMPLEX_INTERACTIONS()
    {Set_Tolerance();}

    void Flush_Intersection_Cache()
    {triple_products.Remove_All();signed_volumes_times_six.Remove_All();}

    void Set_Tolerance(const T tolerance_input=(T)1e2*std::numeric_limits<T>::epsilon())
    {tolerance=tolerance_input;}

//#####################################################################
private:
    VECTOR<T,2> Triple_Product(const VECTOR<TV,3>& locations);
    PAIR<T,bool> Signed_Volume_Times_Six(const VECTOR<TV,4>& locations);
public:
    void Triangle_Segment_Intersection_Weights(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment,VECTOR<T,2>& triangle_weights,T& segment_weight,bool *is_robust_input=0);
    bool Intersection(const VECTOR<TV,4>& tetrahedron,const VECTOR<TV,2>& segment,bool *is_robust_input=0);
    bool Intersection(const VECTOR<TV,4>& tetrahedron,const VECTOR<TV,3>& triangle,bool *is_robust_input=0);
//#####################################################################
};   

template<class T>
class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    HASHTABLE<VECTOR<TV,2>,VECTOR<T,2> > cross_products; // sum of positive and negative terms for each triple of points
    HASHTABLE<VECTOR<TV,3>,PAIR<T,bool> > signed_areas_times_two; // area and indication whether the result is robust or not
    T tolerance;
public:

    ROBUST_SIMPLEX_INTERACTIONS()
    {Set_Tolerance();}

    void Flush_Intersection_Cache()
    {signed_areas_times_two.Remove_All();}

    void Set_Tolerance(const T tolerance_input=(T)1e2*std::numeric_limits<T>::epsilon()) // TODO: change this a more aggressive bound as we have less terms in 2D
    {tolerance=tolerance_input;}

//#####################################################################
private:
    PAIR<T,bool> Signed_Area_Times_Two(const VECTOR<TV,3>& locations);
    VECTOR<T,2> Cross_Product(const VECTOR<TV,2>& locations);
public:
    bool Intersection(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment,bool *is_robust_input=0);
    bool Intersection(const VECTOR<TV,3>& triangle,const VECTOR<TV,1>& point,bool *is_robust_input=0);
    VECTOR<bool,2> Intersection_Test(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment);
    VECTOR<bool,2> Intersection_Test(const VECTOR<TV,3>& triangle,const VECTOR<TV,1>& point);
//#####################################################################
};   

}
#endif

