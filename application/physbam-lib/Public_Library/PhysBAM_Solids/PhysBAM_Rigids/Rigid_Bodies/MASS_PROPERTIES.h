//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Eran Guendelman, Don Hatch, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MASS_PROPERTIES
//#####################################################################
#ifndef __MASS_PROPERTIES__
#define __MASS_PROPERTIES__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
namespace PhysBAM{

template <class TV> class FRAME;

template<class TV,int d=TV::m-1>
class MASS_PROPERTIES:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
    struct DUMMY_IMPLEMENTATION{};struct NORMAL_IMPLEMENTATION{};
    typedef typename IF<INTS_EQUAL<d,0>::value,DUMMY_IMPLEMENTATION,NORMAL_IMPLEMENTATION>::TYPE ACCESS_IMPLEMENTATION;
private:
    const T_SIMPLICIAL_OBJECT& object;
    T mass,density;
    bool use_mass;
    T volume;
    TV center;
    T_WORLD_SPACE_INERTIA_TENSOR inertia_tensor_over_density;
public:

    MASS_PROPERTIES(const T_SIMPLICIAL_OBJECT& object,const bool thin_shell);

    void Set_Mass(const T mass_input)
    {mass=mass_input;use_mass=true;}

    void Set_Density(const T density_input)
    {density=density_input;use_mass=false;}

    T Volume() const
    {return volume;}

    const TV& Center() const
    {return center;}

    void Set_Center(const TV& center_input)
    {center = center_input;}
    
    T Mass() const
    {return use_mass?mass:density*volume;}

    T Density() const
    {return use_mass?mass/volume:density;}

    T_WORLD_SPACE_INERTIA_TENSOR Inertia_Tensor() const
    {return Density()*inertia_tensor_over_density;}

//#####################################################################
    static T Thin_Shell_Volume(const T_SIMPLICIAL_OBJECT& object);
    void Transform_To_Object_Frame(FRAME<TV>& frame,T_INERTIA_TENSOR& object_space_inertia_tensor) const;
    void Transform_To_Object_Frame(FRAME<TV>& frame,T_INERTIA_TENSOR& object_space_inertia_tensor,POINT_CLOUD<TV>& point_cloud) const;
private:
    template<bool thin_shell> void Compute_Properties(MASS_PROPERTIES<TV,d>&,NORMAL_IMPLEMENTATION);
    template<bool thin_shell> void Compute_Properties(MASS_PROPERTIES<TV,d>&,DUMMY_IMPLEMENTATION);
//#####################################################################
};
}
#endif
