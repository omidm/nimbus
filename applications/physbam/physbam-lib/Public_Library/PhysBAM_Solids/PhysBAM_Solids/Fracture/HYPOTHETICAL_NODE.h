//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYPOTHETICAL_NODE
//##################################################################### 
#ifndef __HYPOTHETICAL_NODE__
#define __HYPOTHETICAL_NODE__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_OBJECT.h>
namespace PhysBAM{

template<class TV,int d>
class HYPOTHETICAL_NODE
{
    typedef typename TV::SCALAR T;
public:
    VECTOR<int,2> parents;
    T interpolation_fraction;
    int index_in_embedded_particles;

    HYPOTHETICAL_NODE()
        :interpolation_fraction(0),index_in_embedded_particles(0)
    {}

    HYPOTHETICAL_NODE(const EMBEDDED_OBJECT<TV,d>& embedded_object,const int parent1_input,const int parent2_input,const T interpolation_fraction_input)
        :parents(parent1_input,parent2_input)
    {
        if((index_in_embedded_particles=embedded_object.Embedded_Particle_On_Segment(parents))){
            interpolation_fraction=embedded_object.interpolation_fraction(index_in_embedded_particles);
            if(embedded_object.parent_particles(index_in_embedded_particles)!=parents) exchange(parents[1],parents[2]);}
        else interpolation_fraction=embedded_object.Clamp_Interpolation_Fraction(interpolation_fraction_input);
    }

    bool Is_Parent(const int parent_node) const
    {return parents.Contains(parent_node);}

//#####################################################################    
};
}
#endif
