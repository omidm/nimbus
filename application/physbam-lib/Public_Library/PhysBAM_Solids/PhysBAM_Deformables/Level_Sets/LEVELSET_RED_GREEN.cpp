#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Level_Sets/LEVELSET_RED_GREEN.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_RED_GREEN<TV>::
LEVELSET_RED_GREEN(T_RED_GREEN_GRID& grid_input,ARRAY<T>& phi_input)
    :grid(grid_input),phi(phi_input),overlay_levelset(0),overlay_grid(0),overlay_phi(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_RED_GREEN<TV>::
~LEVELSET_RED_GREEN()
{}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV> void LEVELSET_RED_GREEN<TV>::
Euler_Step(const ARRAY<TV>& velocity,const T dt,const T time)
{
    // use the overlay grid to do the euler step
    LOG::Time("overlay levelset creation");
    Lazy_Update_Overlay_Levelset(&velocity);
    LOG::Time("overlay levelset euler step");
    //overlay_levelset->Euler_Step(overlay_velocity,dt,time); // TODO: update this to work with new octree velocity
    // copy the phi values back to this phi array
    LOG::Time("overlay levelset copy to red green");
    ARRAY<TV>& node_locations=grid.Node_Locations();
    for(int i=1;i<=grid.number_of_nodes;i++)phi(i)=overlay_levelset->Phi(node_locations(i));
    LOG::Stop_Time();
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class TV> void LEVELSET_RED_GREEN<TV>::
Fast_Marching_Method()
{
    LOG::Time("overlay levelset fast march");
    Lazy_Update_Overlay_Levelset();
    overlay_levelset->Fast_Marching_Method();
    
    // copy the phi values back to this phi array
    LOG::Time("overlay levelset copy to red green");
    ARRAY<TV>& node_locations=grid.Node_Locations();
    for(int i=1;i<=grid.number_of_nodes;i++)phi(i)=overlay_levelset->Phi(node_locations(i));
    LOG::Stop_Time();
}
//#####################################################################
// Function Lazy_Update_Overlay_Levelset
//#####################################################################
template<class TV> void LEVELSET_RED_GREEN<TV>::
Lazy_Update_Overlay_Levelset(const ARRAY<TV>* velocity)
{
    if(!overlay_levelset){
        overlay_grid=new T_DYADIC_GRID();
        grid.Create_Overlay_Grid(*overlay_grid,3,true,true);
        overlay_phi=new ARRAY<T>(overlay_grid->number_of_cells);
        overlay_levelset=new T_LEVELSET_DYADIC(*overlay_grid,*overlay_phi);}
    for(int i=1;i<=overlay_grid->number_of_nodes;i++)(*overlay_phi)(i)=Phi(overlay_grid->Node_Location(i));
    if(velocity){
        overlay_velocity.Resize(overlay_grid->number_of_nodes);
        for(int i=1;i<=overlay_grid->number_of_nodes;i++)overlay_velocity(i)=grid.Interpolate_Nodes(*velocity,overlay_grid->Node_Location(i));}
}
//#####################################################################
template class LEVELSET_RED_GREEN<VECTOR<float,2> >;
template class LEVELSET_RED_GREEN<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_RED_GREEN<VECTOR<double,2> >;
template class LEVELSET_RED_GREEN<VECTOR<double,3> >;
#endif
#endif
