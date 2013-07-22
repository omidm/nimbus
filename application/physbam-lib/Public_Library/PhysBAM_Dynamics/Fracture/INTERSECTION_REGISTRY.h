//#####################################################################
// Copyright 2006-2007, Kevin Der, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERSECTION_REGISTRY
//##################################################################### 
#ifndef __INTERSECTION_REGISTRY__
#define __INTERSECTION_REGISTRY__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_SIMPLICES.h>
namespace PhysBAM{

template<class T,int d>
class INTERSECTION_REGISTRY
{
    struct UNUSABLE{};
public:
    CUTTING_SIMPLICES<T,d>& cutting_simplices;
    ARRAY<ARRAY<int> > intersections_on_simplex;
    ARRAY<ARRAY<int> > simplices_on_intersection;
    ARRAY<ARRAY<VECTOR<T,d-1> > > simplex_weights_on_intersection;
    int index_for_last_old_intersection;

    INTERSECTION_REGISTRY(CUTTING_SIMPLICES<T,d>& cutting_simplices_input);

    template<int d2> // must have d2 <= d
    bool Intersection_List(const VECTOR<int,d2>& simplices,ARRAY<int>& intersection_list,typename ENABLE_IF<(d2<=d),UNUSABLE>::TYPE unusable=UNUSABLE());
    template<int d2> // must have d2 <= d
    bool Intersection_List_For_Cuts(const VECTOR<int,d2>& simplices,const VECTOR<int,d+1>& element_nodes,
        ARRAY<int>& intersection_list,typename ENABLE_IF<(d2<=d),UNUSABLE>::TYPE unusable=UNUSABLE());
    template<class T_ARRAY>
    void Register_Intersection(const T_ARRAY& simplices,const typename REBIND<T_ARRAY,VECTOR<T,d-1> >::TYPE& weights,const int particle);
    VECTOR<T,d-1> Get_Simplex_Weights_Of_Intersection(const int intersection,const int simplex);
    bool Intersection_Is_On_Boundary(const int intersection)
    {return cutting_simplices.simplices(simplices_on_intersection(intersection)(1)).parent==0;} // boundary intersection has all children or all parents

    void Resize_Intersections(const int size)
    {simplices_on_intersection.Resize(size);simplex_weights_on_intersection.Resize(size);}

    void Resize_Simplices(const int size)
    {intersections_on_simplex.Resize(size);}

    void Preallocate_Intersections(const int size)
    {simplices_on_intersection.Preallocate(size);simplex_weights_on_intersection.Preallocate(size);}

    void Preallocate_Simplices(const int size)
    {intersections_on_simplex.Preallocate(size);}

    int Number_Of_Intersections(){return simplices_on_intersection.m;}

    void Set_Index_For_Last_Old_Intersection()
    {index_for_last_old_intersection=Number_Of_Intersections();}

    void Print() const
    {std::stringstream ss;
    for(int i=1;i<=simplices_on_intersection.m;i++) ss << "Intersection " << i << " has simplices " << simplices_on_intersection(i) << std::endl;
    LOG::filecout(ss.str());}

    template<class RW>
    void Read(std::istream& input_stream)
    {Read_Binary<RW>(input_stream,cutting_simplices,intersections_on_simplex,simplices_on_intersection,simplex_weights_on_intersection,index_for_last_old_intersection);}

    template<class RW>
    void Write(std::ostream& output_stream) const
    {Write_Binary<RW>(output_stream,cutting_simplices,intersections_on_simplex,simplices_on_intersection,simplex_weights_on_intersection,index_for_last_old_intersection);}

//#####################################################################
    int Intersection(const VECTOR<int,d>& simplices);
//#####################################################################    
};
}
#include <PhysBAM_Dynamics/Read_Write/Fracture/READ_WRITE_INTERSECTION_REGISTRY.h>
#endif
