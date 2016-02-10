//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLES<TV>::
PARTICLES(ARRAY_COLLECTION* array_collection_input)
    :mass(0,0),one_over_mass(0,0),effective_mass(0,0),one_over_effective_mass(0,0),store_mass(false)
{
    delete array_collection;array_collection=array_collection_input;
    Initialize_Array_Collection();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLES<TV>::
PARTICLES()
    :mass(0,0),one_over_mass(0,0),effective_mass(0,0),one_over_effective_mass(0,0),store_mass(false)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLES<TV>::
~PARTICLES()
{}
//#####################################################################
// Function Center_Of_Mass
//#####################################################################
template<class TV> TV PARTICLES<TV>::
Center_Of_Mass() const
{
    if(this->store_mass) return POINT_CLOUDS_COMPUTATIONS::Weighted_Center(this->X,this->mass);
    return this->array_collection->Size()?ARRAYS_COMPUTATIONS::Sum(this->X)/(T)this->array_collection->Size():TV(); // default to treating mass as one
}
//#####################################################################
// Function Compute_Auxiliary_Attributes
//#####################################################################
template<class TV> void PARTICLES<TV>::
Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings)
{
    Compute_Auxiliary_Attributes(soft_bindings,IDENTITY_ARRAY<>(soft_bindings.particles.array_collection->Size()),false);
}
//#####################################################################
// Function Compute_Auxiliary_Attributes
//#####################################################################
template<class TV> template<class T_INDICES> void PARTICLES<TV>::
Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings,const T_INDICES& indices,const bool copy_existing_elements)
{
    if(array_collection->template Get_Array<T>(ATTRIBUTE_ID_ONE_OVER_MASS))
        array_collection->Remove_Array_Using_Index(array_collection->Get_Attribute_Index(ATTRIBUTE_ID_ONE_OVER_MASS));
    array_collection->template Add_Array<T>(ATTRIBUTE_ID_ONE_OVER_MASS,&one_over_mass);
    if(array_collection->template Get_Array<T>(ATTRIBUTE_ID_EFFECTIVE_MASS))
        array_collection->Remove_Array_Using_Index(array_collection->Get_Attribute_Index(ATTRIBUTE_ID_EFFECTIVE_MASS));
    array_collection->template Add_Array<T>(ATTRIBUTE_ID_EFFECTIVE_MASS,&effective_mass);
    if(array_collection->template Get_Array<T>(ATTRIBUTE_ID_ONE_OVER_EFFECTIVE_MASS))
        array_collection->Remove_Array_Using_Index(array_collection->Get_Attribute_Index(ATTRIBUTE_ID_ONE_OVER_EFFECTIVE_MASS));
    array_collection->template Add_Array<T>(ATTRIBUTE_ID_ONE_OVER_EFFECTIVE_MASS,&one_over_effective_mass);
    for(int i=1;i<=indices.Size();i++){int p=indices(i);
        one_over_mass(p)=Robust_Inverse(mass(p));}
    for(int i=1;i<=indices.Size();i++){int p=indices(i);
        one_over_effective_mass(p)=soft_bindings.One_Over_Effective_Mass(p);
        effective_mass(p)=Robust_Inverse(one_over_effective_mass(p));}
}
//#####################################################################
// Function Initialize_Array_Collection
//#####################################################################
template<class TV> void PARTICLES<TV>::
Initialize_Array_Collection()
{
    GEOMETRY_PARTICLES<TV>::Initialize_Array_Collection();
}
//#####################################################################
template class PARTICLES<VECTOR<float,1> >;
template class PARTICLES<VECTOR<float,2> >;
template class PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLES<VECTOR<double,1> >;
template class PARTICLES<VECTOR<double,2> >;
template class PARTICLES<VECTOR<double,3> >;
#endif
}
