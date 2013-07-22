//#####################################################################
// Copyright 2006, Kevin Der, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERSECTION_REGISTRY
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Dynamics/Fracture/INTERSECTION_REGISTRY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> INTERSECTION_REGISTRY<T,d>::
INTERSECTION_REGISTRY(CUTTING_SIMPLICES<T,d>& cutting_simplices_input)
    :cutting_simplices(cutting_simplices_input),index_for_last_old_intersection(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> int INTERSECTION_REGISTRY<T,d>::
Intersection(const VECTOR<int,d>& simplices)
{
    ARRAY<int> intersection_list;
    if(!Intersection_List(simplices,intersection_list)) return 0;
    if(intersection_list.m>1) PHYSBAM_FATAL_ERROR();
    return intersection_list(1);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class T,int d> template<int d2> bool INTERSECTION_REGISTRY<T,d>:: // must have d2 <= d
Intersection_List(const VECTOR<int,d2>& simplices,ARRAY<int>& intersection_list,typename ENABLE_IF<(d2<=d),UNUSABLE>::TYPE unusable)
{
    VECTOR<int,d2> converted_simplices=cutting_simplices.Convert_Simplex_Indices_To_Global_Indices(simplices); // TODO: make a version that takes already converted simplex indices
    // use the intersection list that has the fewest entries to compare the the others agains
    for(int i=2;i<=d2;i++)
        if(intersections_on_simplex(converted_simplices[i]).m<intersections_on_simplex(converted_simplices[1]).m)
            exchange(converted_simplices[1],converted_simplices[i]);
    intersection_list=intersections_on_simplex(converted_simplices[1]);
    for(int i=2;i<=d2;i++)
        for(int j=intersection_list.m;j>=1;j--)
            if(!simplices_on_intersection(intersection_list(j)).Contains(converted_simplices[i]))
                intersection_list.Remove_Index_Lazy(j);
    return intersection_list.m!=0;
}
//#####################################################################
// Function Intersection_List_For_Cuts
//#####################################################################
template<class T,int d> template<int d2> bool INTERSECTION_REGISTRY<T,d>:: // must have d2 <= d
Intersection_List_For_Cuts(const VECTOR<int,d2>& simplices,const VECTOR<int,d+1>& element_nodes,ARRAY<int>& intersection_list,typename ENABLE_IF<(d2<=d),UNUSABLE>::TYPE unusable)
{
    VECTOR<int,d2> simplices_copy=simplices;
    for(int i=2;i<=d2;i++) if(intersections_on_simplex(simplices_copy[i]).m<intersections_on_simplex(simplices_copy[1]).m) exchange(simplices_copy[1],simplices_copy[i]);
    intersection_list=intersections_on_simplex(simplices_copy[1]);
    for(int i=2;i<=d2;i++) for(int j=intersection_list.m;j>=1;j--) if(!simplices_on_intersection(intersection_list(j)).Contains(simplices_copy[i])) intersection_list.Remove_Index_Lazy(j);
    for(int i=intersection_list.m;i>=1;i--){
        const ARRAY<int>& all_simplices=simplices_on_intersection(intersection_list(i));
        bool ok=false;
        for(int j=1;j<=all_simplices.m;j++){ // TODO: why is this necessary?
            VECTOR<int,d> simplex_nodes=cutting_simplices.simplices(all_simplices(j)).nodes;
            if(element_nodes.Contains_All(simplex_nodes)){ok=true;break;}}
        if(!ok) intersection_list.Remove_Index(i);}
    return intersection_list.m!=0;
}
//#####################################################################
// Function Register_Intersection
//#####################################################################
template<class T,int d> template<class T_ARRAY> void INTERSECTION_REGISTRY<T,d>::
Register_Intersection(const T_ARRAY& simplices,const typename REBIND<T_ARRAY,VECTOR<T,d-1> >::TYPE& weights,const int particle)
{
    if(particle>simplices_on_intersection.m) Resize_Intersections(particle);
    int max_simplex=0;for(int i=1;i<=simplices.m;i++) max_simplex=max(max_simplex,simplices(i));
    if(max_simplex>intersections_on_simplex.m) Resize_Simplices(max_simplex);
    if(!simplices_on_intersection(particle).m){
        simplices_on_intersection(particle)=simplices;
        simplex_weights_on_intersection(particle)=weights;
        for(int s=1;s<=simplices.m;s++) intersections_on_simplex(simplices(s)).Append(particle);
        return;}
    for(int i=1;i<=simplices.m;i++){
        if(!simplices_on_intersection(particle).Contains(simplices(i))){
            simplices_on_intersection(particle).Append(simplices(i));
            simplex_weights_on_intersection(particle).Append(weights(i));
            intersections_on_simplex(simplices(i)).Append(particle);}
        else assert(intersections_on_simplex(simplices(i)).Contains(particle));}
}
//#####################################################################
// Function Get_Simplex_Weights_Of_Intersection
//#####################################################################
template<class T,int d> VECTOR<T,d-1> INTERSECTION_REGISTRY<T,d>::
Get_Simplex_Weights_Of_Intersection(const int intersection,const int simplex)
{
    int index=simplices_on_intersection(intersection).Find(simplex);
    if(!index) index=simplices_on_intersection(intersection).Find(cutting_simplices.simplices(simplex).parent);assert(index>0);
    return simplex_weights_on_intersection(intersection)(index);
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
template class INTERSECTION_REGISTRY<T,2>; \
template class INTERSECTION_REGISTRY<T,3>;
INSTANTIATION_HELPER(float)
template bool INTERSECTION_REGISTRY<float,3>::Intersection_List<2>(VECTOR<int,2> const&,ARRAY<int,int>&,ENABLE_IF<(2) <= (3),INTERSECTION_REGISTRY<float,3>::UNUSABLE>::TYPE);
template bool INTERSECTION_REGISTRY<float,2>::Intersection_List_For_Cuts<2>(VECTOR<int,2> const&,VECTOR<int,3> const&,ARRAY<int,int>&,INTERSECTION_REGISTRY<float,2>::UNUSABLE);
template bool INTERSECTION_REGISTRY<float,3>::Intersection_List_For_Cuts<2>(VECTOR<int,2> const&,VECTOR<int,4> const&,ARRAY<int,int>&,INTERSECTION_REGISTRY<float,3>::UNUSABLE);
template void INTERSECTION_REGISTRY<float,2>::Register_Intersection<ARRAY<int,int> >(ARRAY<int,int> const&,REBIND<ARRAY<int,int>,VECTOR<float,1> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<float,2>::Register_Intersection<VECTOR<int,2> >(VECTOR<int,2> const&,REBIND<VECTOR<int,2>,VECTOR<float,1> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<float,3>::Register_Intersection<ARRAY<int,int> >(ARRAY<int,int> const&,REBIND<ARRAY<int,int>,VECTOR<float,2> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<float,3>::Register_Intersection<VECTOR<int,1> >(VECTOR<int,1> const&,REBIND<VECTOR<int,1>,VECTOR<float,2> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<float,3>::Register_Intersection<VECTOR<int,3> >(VECTOR<int,3> const&,REBIND<VECTOR<int,3>,VECTOR<float,2> >::TYPE const&,int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
template bool INTERSECTION_REGISTRY<double,3>::Intersection_List<2>(VECTOR<int,2> const&,ARRAY<int,int>&,ENABLE_IF<(2) <= (3),INTERSECTION_REGISTRY<double,3>::UNUSABLE>::TYPE);
template bool INTERSECTION_REGISTRY<double,2>::Intersection_List_For_Cuts<2>(VECTOR<int,2> const&,VECTOR<int,3> const&,ARRAY<int,int>&,INTERSECTION_REGISTRY<double,2>::UNUSABLE);
template bool INTERSECTION_REGISTRY<double,3>::Intersection_List_For_Cuts<2>(VECTOR<int,2> const&,VECTOR<int,4> const&,ARRAY<int,int>&,INTERSECTION_REGISTRY<double,3>::UNUSABLE);
template void INTERSECTION_REGISTRY<double,2>::Register_Intersection<ARRAY<int,int> >(ARRAY<int,int> const&,REBIND<ARRAY<int,int>,VECTOR<double,1> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<double,2>::Register_Intersection<VECTOR<int,2> >(VECTOR<int,2> const&,REBIND<VECTOR<int,2>,VECTOR<double,1> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<double,3>::Register_Intersection<ARRAY<int,int> >(ARRAY<int,int> const&,REBIND<ARRAY<int,int>,VECTOR<double,2> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<double,3>::Register_Intersection<VECTOR<int,1> >(VECTOR<int,1> const&,REBIND<VECTOR<int,1>,VECTOR<double,2> >::TYPE const&,int);
template void INTERSECTION_REGISTRY<double,3>::Register_Intersection<VECTOR<int,3> >(VECTOR<int,3> const&,REBIND<VECTOR<int,3>,VECTOR<double,2> >::TYPE const&,int);
#endif
