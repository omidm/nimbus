//#####################################################################
// Copyright 2006-2012, Saket Patkar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BUBBLE_PARTICLES
//#####################################################################
#ifndef __BUBBLE_PARTICLES__
#define __BUBBLE_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>


namespace PhysBAM{

template<class TV>
class BUBBLE_PARTICLES: public CLONEABLE<BUBBLE_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef CLONEABLE<BUBBLE_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> > BASE;
public:
    using BASE::array_collection;using BASE::X;using BASE::V;
    ARRAY_VIEW<T> radius;
    ARRAY_VIEW<T> mass;
    ARRAY_VIEW<T> pressure;
    ARRAY_VIEW<T> radial_V;
    ARRAY_VIEW<ARRAY<PAIR<TV_INT,T> > > weights;
    BUBBLE_PARTICLES():radius(0,0),mass(0,0),pressure(0,0),radial_V(0,0),weights(0,0){
        array_collection->Add_Array(ATTRIBUTE_ID_BRADIUS,&radius);
        array_collection->Add_Array(ATTRIBUTE_ID_BMASS,&mass);
        array_collection->Add_Array(ATTRIBUTE_ID_PRESSURE,&pressure);
        array_collection->Add_Array(ATTRIBUTE_ID_RADIAL_VELOCITY,&radial_V);
        array_collection->Add_Array(ATTRIBUTE_ID_WEIGHTS,&weights);
    }
    void Resize(const int new_size){
        array_collection->Resize(new_size);
    }
    
};
template class BUBBLE_PARTICLES<VECTOR<float,1> >;
template class BUBBLE_PARTICLES<VECTOR<float,2> >;
template class BUBBLE_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BUBBLE_PARTICLES<VECTOR<double,1> >;
template class BUBBLE_PARTICLES<VECTOR<double,2> >;
template class BUBBLE_PARTICLES<VECTOR<double,3> >;
#endif
}
#endif
