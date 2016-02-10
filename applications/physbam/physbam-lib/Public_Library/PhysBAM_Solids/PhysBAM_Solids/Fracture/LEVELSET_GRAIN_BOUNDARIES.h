//#####################################################################
// Copyright 2007-2008, Frank Losasso, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_GRAIN_BOUNDARIES
//##################################################################### 
#ifndef __LEVELSET_GRAIN_BOUNDARIES__
#define __LEVELSET_GRAIN_BOUNDARIES__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_GRAIN_BOUNDARIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
namespace PhysBAM{

template<class TV,int d>
class LEVELSET_GRAIN_BOUNDARIES:public FRACTURE_GRAIN_BOUNDARIES<TV,d>
{
    typedef typename TV::SCALAR T;

public:
    using FRACTURE_GRAIN_BOUNDARIES<TV,d>::mesh;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::seed_positions;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::seed_weakness_multipliers;
    using FRACTURE_GRAIN_BOUNDARIES<TV,d>::seed_weakness_multipliers_callback;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::node_smallest_distance;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::is_breakable;
    using FRACTURE_GRAIN_BOUNDARIES<TV,d>::node_region;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::fracture_callbacks;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::is_levelset_grain_boundary;

    LEVELSET_GRAIN_BOUNDARIES(const PARTICLES<TV>& particles,const SIMPLEX_MESH<d>& mesh_input,const ARRAY<TV>& seed_positions_input,const ARRAY<T>& seed_weakness_multipliers_input,const FRAME<TV> frame,
                              const FRACTURE_CALLBACKS<TV>* fracture_callbacks_input=0)
        :FRACTURE_GRAIN_BOUNDARIES<TV,d>(particles,mesh_input,seed_positions_input,seed_weakness_multipliers_input,frame,fracture_callbacks_input)
    {
        is_levelset_grain_boundary=true;
        Calculate_Grain_Boundaries(particles.X,seed_positions,frame,node_smallest_distance,node_region,fracture_callbacks);
    }

    static void Calculate_Grain_Boundaries(ARRAY_VIEW<const TV> input_positions,ARRAY_VIEW<const TV> seed_positions,const FRAME<TV> frame,ARRAY<T>& node_smallest_distance,ARRAY<int>& node_region,
                                           const FRACTURE_CALLBACKS<TV>* fracture_callbacks=0)
    {if(seed_positions.Size()<2) return;

    ARRAY<TV> positions(input_positions);
    int number_of_positions=positions.Size();
    ARRAY<ARRAY<T> > distances(seed_positions.Size());
    ARRAY<TV> perturbation(number_of_positions);if(fracture_callbacks) fracture_callbacks->Perturb(positions,perturbation,frame);
    for(int i=1;i<=seed_positions.Size();i++){
        distances(i).Resize(number_of_positions);
        for(int node=1;node<=number_of_positions;node++){
            distances(i)(node)=(perturbation(node)+seed_positions(i)-positions(node)).Magnitude();}}
    node_smallest_distance.Resize(number_of_positions);node_region.Resize(number_of_positions);
    for(int node=1;node<=number_of_positions;node++){
        int smallest=1;int second_smallest=2;
        for(int i=2;i<=seed_positions.Size();i++){
            if(distances(i)(node)<distances(smallest)(node)){second_smallest=smallest;smallest=i;}
            else if(distances(i)(node)<distances(second_smallest)(node)){second_smallest=i;}}
        node_smallest_distance(node)=max((T)1e-6,distances(second_smallest)(node)-(distances(second_smallest)(node)+distances(smallest)(node))/2);
        node_region(node)=smallest;}}

    void Phi_For_Region_In_Element(const int element,const int region,VECTOR<T,d+1>& phi)
    {for(int i=1;i<=d+1;i++){int node=mesh.elements(element)[i];if(node_region(node)==region)phi[i]=-node_smallest_distance(node);else phi[i]=node_smallest_distance(node);}}

    int Number_Of_Regions()
    {return seed_positions.m;}
//#####################################################################
};  
}
#endif
