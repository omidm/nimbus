//#####################################################################
// Copyright 2006-2007, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROBUST_SIMPLEX_INTERACTIONS
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/permutation.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/ROBUST_SIMPLEX_INTERACTIONS.h>
using namespace PhysBAM;
//#####################################################################
// Function Triple_Product
//#####################################################################
template<class T> VECTOR<T,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Triple_Product(const VECTOR<TV,3>& locations)
{
    VECTOR<T,2> result; // (non_negative_sum,negative_sum)
    for(int i=1;i<=6;i++) if(triple_products.Get(permute_three(locations,i),result)){if(!permutation_of_three_is_even(i)) exchange(result.x,result.y);return result;}
    T product;
    product=locations.x.x*locations.y.y*locations.z.z;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.y.x*locations.z.y*locations.x.z;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.z.x*locations.x.y*locations.y.z;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.x.x*locations.z.y*locations.y.z;if(product<0) result.x-=product;else result.y+=product;
    product=locations.y.x*locations.x.y*locations.z.z;if(product<0) result.x-=product;else result.y+=product;
    product=locations.z.x*locations.y.y*locations.x.z;if(product<0) result.x-=product;else result.y+=product;
    triple_products.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Signed_Volume_Times_Six
//#####################################################################
template<class T> PAIR<T,bool> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Signed_Volume_Times_Six(const VECTOR<TV,4>& locations)
{
    PAIR<T,bool> result;
    for(int i=1;i<=24;i++) if(signed_volumes_times_six.Get(permute_four(locations,i),result)){if(!permutation_of_four_is_even(i)) result.x=-result.x;return result;}
    VECTOR<T,2> determinant=Triple_Product(VECTOR<TV,3>(locations[2],locations[3],locations[4]))+Triple_Product(VECTOR<TV,3>(locations[1],locations[4],locations[3]))+
        Triple_Product(VECTOR<TV,3>(locations[1],locations[2],locations[4]))+Triple_Product(VECTOR<TV,3>(locations[1],locations[3],locations[2]));
    result.x=TV::Triple_Product(locations(2)-locations(1),locations(3)-locations(1),locations(4)-locations(1));
    if(abs(determinant.x-determinant.y)<tolerance*(determinant.x+determinant.y)) result.y=false;
    else{result.y=true;
        bool positive1=(determinant.x-determinant.y>0),positive2=(result.x>0);
        if(positive1^positive2) PHYSBAM_FATAL_ERROR();}
    signed_volumes_times_six.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Triangle_Segment_Intersection_Weights
//#####################################################################
template<class T> void ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::    
Triangle_Segment_Intersection_Weights(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment,VECTOR<T,2>& triangle_weights,T& segment_weight,bool *is_robust_input)
{
    MATRIX<T,3> matrix(triangle(1)-triangle(3),triangle(2)-triangle(3),segment(2)-segment(1));
    VECTOR<T,3> weights=matrix.Robust_Solve_Linear_System(segment(2)-triangle(3));
    for(int i=1;i<=2;i++) triangle_weights(i)=weights(i);segment_weight=weights(3);
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Intersection(const VECTOR<TV,4>& tetrahedron,const VECTOR<TV,2>& segment,bool *is_robust_input)
{
    bool is_robust=true;
    // Separator plane must contain two tetrahedron vertices and one segment vertex
    for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++) for(int k=1;k<=2;k++){
        PAIR<T,bool> volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),segment(k),segment(3-k)));
        if(!volume1.y){is_robust=false;continue;}
        for(int l=1;l<=4;l++) if(l!=i && l!=j){
            PAIR<T,bool> volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),segment(k),tetrahedron(l)));
            if(!volume2.y){is_robust=false;goto Next_Plane;}
            if((volume1.x>0)^(volume2.x<0)) goto Next_Plane;}
        if(is_robust_input) *is_robust_input=true;return false;
        Next_Plane:;}
    // No separator plane found, objects are intersecting
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Intersection(const VECTOR<TV,4>& tetrahedron,const VECTOR<TV,3>& triangle,bool *is_robust_input)
{
    bool is_robust=true;
    PAIR<T,bool> volume1,volume2;
    // Case 1 : Entire tetrahedron on same half-space of triangle plane
    volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(1),triangle(2),triangle(3),tetrahedron(1)));
    if(!volume1.y) is_robust=false;
    else for(int i=2;i<=4;i++){
        volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(1),triangle(2),triangle(3),tetrahedron(i)));if(!volume2.y){is_robust=false;break;}
        if((volume1.x>0)^(volume2.x>0)) break;
        if(i==4){if(is_robust_input) *is_robust_input=true;return false;}}
    // Case 2 : Separator plane contains two triangle vertices
    for(int i=1;i<=2;i++) for(int j=i+1;j<=3;j++) for(int k=1;k<=4;k++){
        volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(i),triangle(j),tetrahedron(k),triangle(6-i-j)));if(!volume1.y){is_robust=false;continue;}
        for(int l=1;l<=4;l++) if(l!=k){
            volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(i),triangle(j),tetrahedron(k),tetrahedron(l)));
            if(!volume2.y){is_robust=false;goto Next_Plane1;}
            if((volume1.x>0)^(volume2.x<0)) goto Next_Plane1;}
        if(is_robust_input) *is_robust_input=true;return false;
        Next_Plane1:;}
    // Case 3 : Separator plane contains two tetraherdon vertices
    for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++) for(int k=1;k<=3;k++){
        VECTOR<int,2> indices=VECTOR<int,4>(1,2,3,4).Remove_Index(j).Remove_Index(i);
        volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),triangle(k),tetrahedron(indices(1))));if(!volume1.y){is_robust=false;continue;}
        volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),triangle(k),tetrahedron(indices(2))));if(!volume2.y){is_robust=false;continue;}
        if((volume1.x>0)^(volume2.x>0)) continue;
        for(int l=1;l<=3;l++) if(l!=k){
            volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),triangle(k),triangle(l)));if(!volume2.y){is_robust=false;goto Next_Plane2;}
            if((volume1.x>0)^(volume2.x<0)) goto Next_Plane2;}
        if(is_robust_input) *is_robust_input=true;return false;
        Next_Plane2:;}
    // No separator plane found, objects are intersecting
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Cross_Product
//#####################################################################
template<class T> VECTOR<T,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Cross_Product(const VECTOR<TV,2>& locations)
{
    VECTOR<T,2> result; // (non_negative_sum,negative_sum)
    for(int i=1;i<=2;i++) if(cross_products.Get(permute_two(locations,i),result)){if(!permutation_of_two_is_even(i)) exchange(result.x,result.y);return result;}
    T product;
    product=locations.x.x*locations.y.y;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.x.y*locations.y.x;if(product<0) result.x-=product;else result.y+=product;
    cross_products.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Signed_Area_Times_Two
//#####################################################################
template<class T> PAIR<T,bool> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Signed_Area_Times_Two(const VECTOR<TV,3>& locations)
{
    PAIR<T,bool> result;
    for(int i=1;i<=6;i++) if(signed_areas_times_two.Get(permute_three(locations,i),result)){
        if(!permutation_of_three_is_even(i)) result.x=-result.x;return result;}
    VECTOR<T,2> determinant=Cross_Product(VECTOR<TV,2>(locations[2],locations[3]))+Cross_Product(VECTOR<TV,2>(locations[3],locations[1]))+Cross_Product(VECTOR<TV,2>(locations[1],locations[2]));
    result.x=TV::Cross_Product(locations[2]-locations[1],locations[3]-locations[1]).x;
    if(abs(determinant.x-determinant.y)<tolerance*(determinant.x+determinant.y)) result.y=false;
    else{result.y=true;
        bool positive1=(determinant.x-determinant.y>0),positive2=(result.x>0);
        if(positive1^positive2) PHYSBAM_FATAL_ERROR();}
    signed_areas_times_two.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment,bool *is_robust_input)
{
    bool is_robust=true;
    PAIR<T,bool> area1,area2;
    // Case 1 : Entire triangle on same half-space of segment line
    area1=Signed_Area_Times_Two(VECTOR<TV,3>(segment(1),segment(2),triangle(1)));
    if(!area1.y) is_robust=false;
    else for(int i=2;i<=3;i++){
        area2=Signed_Area_Times_Two(VECTOR<TV,3>(segment(1),segment(2),triangle(i)));if(!area2.y){is_robust=false;break;}
        if((area1.x>0)^(area2.x>0)) break;
        if(i==4){if(is_robust_input) *is_robust_input=true;return false;}}
    // Case 2 : try combinations of one point on triangle and one point on segment
    for(int segment_i=1;segment_i<=2;segment_i++) for(int triangle_i=1;triangle_i<=3;triangle_i++){
        VECTOR<int,2> other_triangle_indices(triangle_i%3+1,(triangle_i%3+1)%3+1);
        int other_segment_i=3-segment_i;
        area1=Signed_Area_Times_Two(VECTOR<TV,3>(segment(segment_i),triangle(triangle_i),segment(other_segment_i)));if(!area1.y){is_robust=false;continue;}
        for(int dummy=1;dummy<=2;dummy++){int other_triangle_i=other_triangle_indices(dummy);
            area2=Signed_Area_Times_Two(VECTOR<TV,3>(segment(segment_i),triangle(triangle_i),triangle(other_triangle_i)));if(!area2.y){is_robust=false;goto Next_Line;}
            if((area1.x>0)^(area2.x<0)) goto Next_Line;}
        if(is_robust_input) *is_robust_input=true;return false;
      Next_Line:;}
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection(const VECTOR<TV,3>& triangle,const VECTOR<TV,1>& point,bool *is_robust_input)
{
    bool is_robust=true;
    PAIR<T,bool> area1,area2;
    // Iterate over all edges
    for(int triangle_1=1;triangle_1<=2;triangle_1++) for(int triangle_2=triangle_1+1;triangle_2<=3;triangle_2++){
        area1=Signed_Area_Times_Two(VECTOR<TV,3>(triangle(triangle_1),triangle(triangle_2),point(1)));if(!area1.y){is_robust=false;continue;}
        int other_triangle_index=VECTOR<int,3>(1,2,3).Remove_Index(triangle_2).Remove_Index(triangle_1)[1];
        area2=Signed_Area_Times_Two(VECTOR<TV,3>(triangle(triangle_1),triangle(triangle_2),triangle(other_triangle_index)));if(!area2.y){is_robust=false;continue;}
        if((area1.x>0)^(area2.x<0)) continue;
        if(is_robust_input) *is_robust_input=true;return false;}
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Intersection_Test
//#####################################################################
template<class T> VECTOR<bool,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection_Test(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment)
{
    bool robust;
    bool value=Intersection(triangle,segment,&robust);
    return VECTOR<bool,2>(value,robust);
}
//#####################################################################
// Function Intersection_Test
//#####################################################################
template<class T> VECTOR<bool,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection_Test(const VECTOR<TV,3>& triangle,const VECTOR<TV,1>& point)
{
    bool robust;
    bool value=Intersection(triangle,point,&robust);
    return VECTOR<bool,2>(value,robust);
}
//####################################################################
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<float,2> >;
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,2> >;
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,3> >;
#endif
