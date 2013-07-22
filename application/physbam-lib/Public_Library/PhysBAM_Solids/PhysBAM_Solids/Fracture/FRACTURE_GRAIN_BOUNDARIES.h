//#####################################################################
// Copyright 2007, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_GRAIN_BOUNDARIES
//##################################################################### 
#ifndef __FRACTURE_GRAIN_BOUNDARIES__
#define __FRACTURE_GRAIN_BOUNDARIES__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
namespace PhysBAM{

template<class TV,int d>
class FRACTURE_GRAIN_BOUNDARIES:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;

public:
    const SIMPLEX_MESH<d>& mesh;
    ARRAY<TV> seed_positions;
    ARRAY<T> seed_weakness_multipliers;
    ARRAY<T> seed_weakness_multipliers_callback;
    ARRAY<T> node_smallest_distance;
    ARRAY<bool> is_breakable;
    ARRAY<int> node_region;
    const FRACTURE_CALLBACKS<TV>* fracture_callbacks;
    bool is_levelset_grain_boundary;

    FRACTURE_GRAIN_BOUNDARIES(const PARTICLES<TV>& particles,const SIMPLEX_MESH<d>& mesh_input,const ARRAY<TV>& seed_positions_input,const ARRAY<T>& seed_weakness_multipliers_input,const FRAME<TV> frame,
                              const FRACTURE_CALLBACKS<TV>* fracture_callbacks_input=0)
        :mesh(mesh_input),seed_positions(seed_positions_input),seed_weakness_multipliers(seed_weakness_multipliers_input),fracture_callbacks(fracture_callbacks_input),is_levelset_grain_boundary(false)
    {}

    virtual ~FRACTURE_GRAIN_BOUNDARIES()
    {}

    void Set_Fracture_Callbacks(const FRACTURE_CALLBACKS<TV>* fracture_callbacks_input)
    {
        fracture_callbacks=fracture_callbacks_input;
    }

    bool Element_Spanning_Regions(const int element)
    {for(int i=2;i<=d+1;i++)if(node_region(mesh.elements(element)[1])!=node_region(mesh.elements(element)[i]))return true;return false;}

    int Regions_Intersecting_Element(const int element,VECTOR<int,d+1>& regions)
    {int number_of_regions=0;
    for(int i=1;i<=d+1;i++){
        bool add_region=true;int node=mesh.elements(element)[i];
        for(int j=1;j<=number_of_regions;j++)if(node_region(node)==regions(j)){add_region=false;break;}
        if(add_region){number_of_regions++;regions(number_of_regions)=node_region(node);}}
    return number_of_regions;}

    void Initialize_Element_Weakness_Multiplier(const ARRAY<TV>& tet_centers)
    {if(!fracture_callbacks) return;
    seed_weakness_multipliers_callback.Resize(tet_centers.m);
    fracture_callbacks->Grain_Boundary_Weakness_Multiplier(tet_centers,seed_weakness_multipliers_callback);}

    T Element_Weakness_Multiplier(const int tet_index, const int number_of_regions,const VECTOR<int,d+1>& regions)
    {T weakness=-1;
    if(fracture_callbacks) weakness=seed_weakness_multipliers_callback(tet_index);
    if(weakness != -1) return weakness;
    weakness=0;
    for(int i=1;i<=number_of_regions;i++)weakness+=seed_weakness_multipliers(regions[i]);return weakness/number_of_regions;}

    void Initialize_Breakability(const ARRAY<TV>& tet_centers) // might be more efficient to combine with Initialize function above
    {if(!fracture_callbacks) return;
    is_breakable.Resize(tet_centers.m);
    fracture_callbacks->Get_Breakability(tet_centers,is_breakable);}

    virtual void Phi_For_Region_In_Element(const int element,const int region,VECTOR<T,d+1>& phi){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    int Number_Of_Nodes_In_Region(const int region)
    {int count=0;for(int i=1;i<=mesh.number_nodes;i++)if(node_region(i)==region)count++;return count;}

    virtual int Number_Of_Regions(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}

//#####################################################################
};  
}
#endif
