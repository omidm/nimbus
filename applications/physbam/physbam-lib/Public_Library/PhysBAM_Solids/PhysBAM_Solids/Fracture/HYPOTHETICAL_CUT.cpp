//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYPOTHETICAL_CUT
//##################################################################### 
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_CUT.h>
using namespace PhysBAM;
//#####################################################################
// Function Contains_Embedded_Particle
//##################################################################### 
template<class TV,int d> bool HYPOTHETICAL_CUT<TV,d>::
Contains_Embedded_Particle(const EMBEDDED_OBJECT<TV,d>& embedded_object,const int emb_particle) const
{
    for(int hp=1;hp<=hypothetical_nodes.m;hp++)
        if(embedded_object.Are_Parents(hypothetical_nodes(hp).parents,emb_particle)) return true;
    return false;
}
//#####################################################################
// Function Position
//##################################################################### 
template<class TV,int d> TV HYPOTHETICAL_CUT<TV,d>::
Position(const int hypothetical_node_index) const
{
    return LINEAR_INTERPOLATION<T,TV>::Linear(embedded_object.particles.X(hypothetical_nodes(hypothetical_node_index).parents[1]),
        embedded_object.particles.X(hypothetical_nodes(hypothetical_node_index).parents[2]),hypothetical_nodes(hypothetical_node_index).interpolation_fraction);
}
//#####################################################################
// Function Add_Hypothetical_Nodes_To_Embedded_Object
//##################################################################### 
template<class TV,int d> void HYPOTHETICAL_CUT<TV,d>::
Add_Hypothetical_Nodes_To_Embedded_Object(EMBEDDED_OBJECT<TV,d>& embedded_object)
{
    for(int n=1;n<=hypothetical_nodes.m;n++)
        hypothetical_nodes(n).index_in_embedded_particles=embedded_object.Add_Embedded_Particle_If_Not_Already_There(hypothetical_nodes(n).parents,hypothetical_nodes(n).interpolation_fraction);
}
//#####################################################################
template class HYPOTHETICAL_CUT<VECTOR<float,2>,2>;
template class HYPOTHETICAL_CUT<VECTOR<float,3>,2>;
template class HYPOTHETICAL_CUT<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HYPOTHETICAL_CUT<VECTOR<double,2>,2>;
template class HYPOTHETICAL_CUT<VECTOR<double,3>,2>;
template class HYPOTHETICAL_CUT<VECTOR<double,3>,3>;
#endif
