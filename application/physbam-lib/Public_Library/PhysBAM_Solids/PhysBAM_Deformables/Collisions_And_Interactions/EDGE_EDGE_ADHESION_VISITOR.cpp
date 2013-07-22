//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EDGE_EDGE_ADHESION_VISITOR
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/EDGE_EDGE_ADHESION_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_ADHESION.h>
using namespace PhysBAM;
//#####################################################################
// Function Find_Place_For_Spring
//#####################################################################
template<class TV> bool EDGE_EDGE_ADHESION_VISITOR<TV>::
Find_Place_For_Spring(const int reference_segment_index,const int other_segment_index,const T distance)
{
    if(adhesion.segments_with_springs(reference_segment_index).x.m>adhesion.max_connections){ // out of room or just ran out
        ARRAY<PAIR<T,int> >& heap=adhesion.segments_with_springs(reference_segment_index).x;
        if(adhesion.segments_with_springs(reference_segment_index).y==false){ARRAYS_COMPUTATIONS::Heapify(heap);adhesion.segments_with_springs(reference_segment_index).y=true;} // make into heap
        if(heap(1).x<=distance) return false; // not allowed, biggest is smaller than our candidate
        else{ // allowed, so kick out some existing pair
            adhesion.springs->Delete(VECTOR<int,2>(reference_segment_index,heap(1).y).Sorted());

            const VECTOR<int,2>& segment1_nodes=adhesion.curve.mesh.elements(reference_segment_index),&segment2_nodes=adhesion.curve.mesh.elements(heap(1).y);
            adhesion.intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]));
            adhesion.intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment2_nodes[1],segment2_nodes[2],segment1_nodes[1],segment1_nodes[2]));

            heap(1)=PAIR<T,int>(distance,other_segment_index);ARRAYS_COMPUTATIONS::Heapify(heap,1,heap.m);
            return true;}}
    else{ // plenty of room
        adhesion.segments_with_springs(reference_segment_index).x.Append(PAIR<T,int>(distance,other_segment_index));
        return true;}
}
//#####################################################################
// Function Store
//#####################################################################
template<class TV> void EDGE_EDGE_ADHESION_VISITOR<TV>::
Store(const int segment1_local_index,const int segment2_local_index)
{
    typedef typename SEGMENT_ADHESION<TV>::SPRING_STATE T_SPRING_STATE;
    int segment1_index=mesh1_indices(segment1_local_index),segment2_index=mesh2_indices(segment2_local_index);
    VECTOR<int,2> segment_indices=VECTOR<int,2>(segment1_index,segment2_index).Sorted();
    const VECTOR<int,2> &segment1_nodes=adhesion.curve.mesh.elements(segment_indices[1]),&segment2_nodes=adhesion.curve.mesh.elements(segment_indices[2]);

    if(adhesion.existing_pairs.Contains(segment_indices)) return; // cull initially close pairs
    if(!local && segment1_index>segment2_index) return; // processor that has the lowest segment # owns it
    if(segment1_nodes.Contains(segment2_nodes[1]) || segment1_nodes.Contains(segment2_nodes[2])) return; // neighbor segments
    if(adhesion.springs->Contains(segment_indices)) return; // already have adhesion pair

    // check for interaction
    T_SPRING_STATE state;
    SEGMENT_3D<T> segment1(adhesion.curve.particles.X.Subset(segment1_nodes)),segment2(adhesion.curve.particles.X.Subset(segment2_nodes));
    if(!segment1.Edge_Edge_Interaction(segment2,adhesion.on_distance,state.distance,state.normal,state.weights,true)) return;

    // TODO: this branch should be unnecessary, we should always be able to use segment_indices[1] as the reference.
    PHYSBAM_ASSERT(local || segment_indices[1]==segment1_index); // TODO: make this a assert()
    //int reference_segment_index,other_segment_index;
    //if(!local){reference_segment_index=segment1_index;other_segment_index=segment2_index;} // if this pair strides processors, use local segment to decide room
    //else{reference_segment_index=segment_indices[1];other_segment_index=segment_indices[2];}
        
    //T distance=state.distance;
    T distance=state.distance*(2-abs(TV::Dot_Product((adhesion.curve.particles.X(segment1_nodes[2])-adhesion.curve.particles.X(segment1_nodes[1])).Normalized(),
                                                         (adhesion.curve.particles.X(segment2_nodes[2])-adhesion.curve.particles.X(segment2_nodes[1])).Normalized())));
    {std::stringstream ss;ss<<"FOUND SPRING "<<segment1_index<<" segment2 "<<segment2_index<<std::endl;LOG::filecout(ss.str());}
    //if(Find_Place_For_Spring(reference_segment_index,other_segment_index,state.distance)){
    if(Find_Place_For_Spring(segment_indices[1],segment_indices[2],distance)){
        // compute rest of spring parameters
        state.external=!local;
        state.nodes=segment1_nodes.Append_Elements(segment2_nodes);
        T harmonic_mass=Pseudo_Inverse((T).25*(one_over_effective_mass(state.nodes[1])+one_over_effective_mass(state.nodes[2])
                +one_over_effective_mass(state.nodes[3])+one_over_effective_mass(state.nodes[4])));
        state.damping=adhesion.overdamping_fraction*2*sqrt(adhesion.youngs_modulus*adhesion.restlength*harmonic_mass);
        adhesion.springs->Insert(segment_indices,state);
        // add to collision ignore
        adhesion.intersecting_edge_edge_pairs.Set(VECTOR<int,4>(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]));
        adhesion.intersecting_edge_edge_pairs.Set(VECTOR<int,4>(segment2_nodes[1],segment2_nodes[2],segment1_nodes[1],segment1_nodes[2]));
    }
}
//####################################################################
#define INSTANTIATION_HELPER(T) \
    template void EDGE_EDGE_ADHESION_VISITOR<VECTOR<T,3> >::Store(int,int); \
    template void BOX_HIERARCHY<VECTOR<T,3> >::Intersection_List<EDGE_EDGE_ADHESION_VISITOR<VECTOR<T,3> >,ZERO>(BOX_HIERARCHY<VECTOR<T,3> > const&, \
        EDGE_EDGE_ADHESION_VISITOR<VECTOR<T,3> >&,ZERO) const;

INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif
