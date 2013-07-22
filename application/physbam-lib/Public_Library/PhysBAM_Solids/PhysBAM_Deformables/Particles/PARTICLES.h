//#####################################################################
// Copyright 2004-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLES
//#####################################################################
#ifndef __PARTICLES__
#define __PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/CENTER.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/EULER_STEP.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>

namespace PhysBAM{

template<class TV> class SOFT_BINDINGS; 

template<class TV>
class PARTICLES:public CLONEABLE<PARTICLES<TV>,GEOMETRY_PARTICLES<TV> > // X, V
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<PARTICLES<TV>,GEOMETRY_PARTICLES<TV> > BASE;
public:
    using BASE::array_collection;using BASE::X;using BASE::V;

    ARRAY_VIEW<T> mass;
    ARRAY_VIEW<T> one_over_mass;
    ARRAY_VIEW<T> effective_mass;
    ARRAY_VIEW<T> one_over_effective_mass;
    bool store_mass;

    PARTICLES(ARRAY_COLLECTION* array_collection_input);
    PARTICLES();
    virtual ~PARTICLES();

    void Store_Mass(bool store=true)
    {store_mass=store;if(store) array_collection->Add_Array(ATTRIBUTE_ID_MASS,&mass);else array_collection->Remove_Array(ATTRIBUTE_ID_MASS);}

    T Min_Mass() const 
    {return mass.Size()?ARRAYS_COMPUTATIONS::Min(mass):FLT_MAX;}

    T Max_Mass() const 
    {return mass.Size()?ARRAYS_COMPUTATIONS::Max(mass):0;}

    void Euler_Step_Position(const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(X,V,dt);}

    template<class T_INDICES>
    void Euler_Step_Position(const T_INDICES& indices,const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(indices,X,V,dt);}

    TV Center_Of_Mass()
    {return POINT_CLOUDS_COMPUTATIONS::Weighted_Center(X,mass);}

//#####################################################################
    TV Center_Of_Mass() const;
    void Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings);
    template<class T_INDICES> void Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings,const T_INDICES& indices,const bool copy_existing_elements=true);
    void Initialize_Array_Collection();
//#####################################################################
};
}
#endif
