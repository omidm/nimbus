//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_PARTICLES
//#####################################################################
#ifndef __RIGID_BODY_PARTICLES__
#define __RIGID_BODY_PARTICLES__

#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_PARTICLES:public CLONEABLE<RIGID_BODY_PARTICLES<TV>,RIGID_GEOMETRY_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<RIGID_BODY_PARTICLES<TV>,RIGID_GEOMETRY_PARTICLES<TV> > BASE;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
public:
    using BASE::rigid_geometry;using BASE::array_collection;using BASE::Delete_All_Particles;

    ARRAY_VIEW<typename TV::SPIN> angular_momentum;
    ARRAY_VIEW<T> mass;
    ARRAY_VIEW<T_INERTIA_TENSOR> inertia_tensor;
    ARRAY_VIEW<bool> kinematic;

    RIGID_BODY_PARTICLES(ARRAY_COLLECTION* array_collection);
    RIGID_BODY_PARTICLES();
    virtual ~RIGID_BODY_PARTICLES();

private:
    void Delete_Particle(const int index)
    {PHYSBAM_FATAL_ERROR();}

    void Delete_Particles_On_Deletion_List(const bool preserve_order=false)
    {PHYSBAM_FATAL_ERROR();}

    void Add_To_Deletion_List(const int index)
    {PHYSBAM_FATAL_ERROR();}
public:
    void Remove_Body(const int p)
    {BASE::Remove_Geometry(p);}

//#####################################################################
    void Initialize_Array_Collection();
//#####################################################################
};
}
#endif
