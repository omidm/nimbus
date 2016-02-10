//#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDING
//#####################################################################
#ifndef __EMBEDDING__
#define __EMBEDDING__

#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV> class BINDING_LIST;

template<class TV>
class EMBEDDING:public STRUCTURE<TV>
{
    STATIC_ASSERT((TV::m>1));

    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
public:
    GEOMETRY_PARTICLES<TV>& particles;
    TRIANGLE_MESH material_surface_mesh;
    T_TRIANGULATED_OBJECT material_surface;
private:
    bool need_destroy_particles;
public:

    EMBEDDING(GEOMETRY_PARTICLES<TV>& particles_input);

    virtual ~EMBEDDING()
    {if(need_destroy_particles) delete &particles;}

    static EMBEDDING<TV>* Create()
    {EMBEDDING<TV>* embedding=new EMBEDDING<TV>(*new GEOMETRY_PARTICLES<TV>);embedding->need_destroy_particles=true;return embedding;}

    static EMBEDDING<TV>* Create(GEOMETRY_PARTICLES<TV>& new_particles)
    {EMBEDDING<TV>* embedding=new EMBEDDING<TV>(new_particles);return embedding;}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("EMBEDDING<VECTOR<T,%d> >",TV::dimension);}

    void Update_Number_Nodes() PHYSBAM_OVERRIDE
    {material_surface.Update_Number_Nodes();}

    void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const PHYSBAM_OVERRIDE
    {material_surface.Mark_Nodes_Referenced(marks,mark);}

public:

//#####################################################################
};

template<class T_INPUT>
class EMBEDDING<VECTOR<T_INPUT,1> >:public STRUCTURE<VECTOR<T_INPUT,1> >
{
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T_INPUT,1> >::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
public:
    TRIANGLE_MESH material_surface_mesh;
    T_TRIANGULATED_OBJECT material_surface;
//#####################################################################
};
}
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Fracture/READ_WRITE_EMBEDDING.h>
#endif
