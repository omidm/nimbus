#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_SIMPLICIAL_OBJECT
//#####################################################################
#ifndef __LEVELSET_SIMPLICIAL_OBJECT__
#define __LEVELSET_SIMPLICIAL_OBJECT__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_SUBSET.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Level_Sets/LEVELSET_RED_GREEN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV,int d>
class LEVELSET_SIMPLICIAL_OBJECT:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename RED_GREEN_POLICY<VECTOR<T,d> >::GRID_T T_RED_GREEN_GRID;
    typedef typename MESH_POLICY<d>::MESH T_MESH;
public:
    typedef int HAS_TYPED_READ_WRITE;

    T_RED_GREEN_GRID grid;
    EMBEDDED_MATERIAL_SURFACE<TV,d>& embedding;
    POINT_CLOUD_SUBSET<TV,GEOMETRY_PARTICLES<TV> > particles;
    ARRAY<int> cell_to_simplex_mapping;
    ARRAY<int> node_to_particle_mapping;
    ARRAY<T> phi;
    LEVELSET_RED_GREEN<VECTOR<T,d> > levelset;
    ARRAY<VECTOR<T,d> > V;
private:
    bool need_destroy_embedding;
public:

    LEVELSET_SIMPLICIAL_OBJECT(EMBEDDED_MATERIAL_SURFACE<TV,d>& embedding_input);
    ~LEVELSET_SIMPLICIAL_OBJECT();

//#####################################################################
    static LEVELSET_SIMPLICIAL_OBJECT* Create();
    void Initialize();
    void Build_Embedded_Object(const bool verbose);
    ARRAY<TV> Get_Material_Coordinates();
    void Euler_Step(const T& dt,const T& time);
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
//#####################################################################
};
}
#endif
#endif
