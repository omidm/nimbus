//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_MATERIAL_SURFACE
//#####################################################################
#ifndef __EMBEDDED_MATERIAL_SURFACE__
#define __EMBEDDED_MATERIAL_SURFACE__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
namespace PhysBAM{

template<class TV> class DEFORMABLE_BODY_COLLECTION;

template<class TV,int d>
class EMBEDDED_MATERIAL_SURFACE:public EMBEDDING<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
    typedef EMBEDDING<TV> BASE;
    using BASE::particles;using BASE::material_surface_mesh;using BASE::material_surface;

    T_EMBEDDED_OBJECT& embedded_object;
    POINT_CLOUD_SUBSET<TV,GEOMETRY_PARTICLES<TV> >& embedded_particles;
    ARRAY<VECTOR<int,2> >& parent_particles;
    ARRAY<T>& interpolation_fraction;
    ARRAY<bool> previously_perturbed;
    int number_of_soft_bound_particles;
    ARRAY<int> soft_particles;
private:
    bool need_destroy_embedded_object;
public:

    EMBEDDED_MATERIAL_SURFACE(T_EMBEDDED_OBJECT& embedded_object);
    virtual ~EMBEDDED_MATERIAL_SURFACE();

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("EMBEDDED_MATERIAL_SURFACE<VECTOR<T,%d>,%d>",TV::dimension,d);}

    void Update_Number_Nodes() PHYSBAM_OVERRIDE
    {embedded_object.Update_Number_Nodes();material_surface.Update_Number_Nodes();}

public:

//#####################################################################
    static typename EMBEDDING_POLICY<TV,d>::EMBEDDING* Create();
    static typename EMBEDDING_POLICY<TV,d>::EMBEDDING* Create(GEOMETRY_PARTICLES<TV>& new_particles);
    STRUCTURE<TV>* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& particles,ARRAY<int>* particle_indices=0) const PHYSBAM_OVERRIDE;
    virtual void Create_Material_Surface(const bool verbose=true);
    void Update_Binding_List_From_Embedding(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection);
    virtual void Perturb_Nodes_For_Collision_Freeness(const T perturb_amount)=0;
private:
    virtual void Construct_Material_Surface_Mesh()=0;
protected:
//#####################################################################
};
}
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Fracture/READ_WRITE_EMBEDDED_MATERIAL_SURFACE.h>
#endif
