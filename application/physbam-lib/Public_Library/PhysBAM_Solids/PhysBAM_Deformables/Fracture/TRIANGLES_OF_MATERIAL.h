//#####################################################################
// Copyright 2004-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLES_OF_MATERIAL
//##################################################################### 
#ifndef __TRIANGLES_OF_MATERIAL__
#define __TRIANGLES_OF_MATERIAL__

#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
namespace PhysBAM{

template<class TV>
class TRIANGLES_OF_MATERIAL:public EMBEDDED_MATERIAL_SURFACE<TV,2>
{
    typedef typename TV::SCALAR T;

    typedef EMBEDDED_MATERIAL_SURFACE<TV,2> BASE;
    typedef typename BASE::T_EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
public:
    using BASE::material_surface_mesh;using BASE::material_surface;using BASE::embedded_object;using BASE::previously_perturbed;using BASE::particles;
    using BASE::embedded_particles;

    TRIANGLES_OF_MATERIAL(T_EMBEDDED_OBJECT& embedded_object_input);

//#####################################################################
    void Perturb_Nodes_For_Collision_Freeness(const T perturb_amount=1e-6) PHYSBAM_OVERRIDE;
private:
    void Construct_Material_Surface_Mesh() PHYSBAM_OVERRIDE;
    void Add_To_Material_Surface_Mesh_Triangle(const int material_node1,const int material_node2,const int material_node3);
    void Add_To_Material_Surface_Mesh_Quad(const int x1,const int x2,const int x3,const int x4);
    void Add_To_Material_Surface_Mesh_Face_Triangle(const int triangle);
    void Add_To_Material_Surface_Mesh_Subquadrilateral_Containg_Diamond_Node(const int triangle);
    void Add_To_Material_Surface_Mesh_Diamond_Quad(const int tri_node,const int emb_node1,const int emb_node2,const int emb_node3);
    void Add_To_Material_Surface_Mesh_Corner_Triangles(const int triangle);
    void Add_To_Material_Surface_Mesh_Corner_Triangle(const int emb_node1,const int tri_node,const int emb_node2);
    void Add_To_Material_Surface_Mesh_Isolated_Node_Subtriangle(const int triangle);
    void Add_To_Material_Surface_Mesh_Subtriangle(const int curve_particle1,const int curve_particle2,const int triangle_particle);
    void Add_To_Material_Surface_Mesh_Subquadrilateral_Opposite_Isolated_Node(const int triangle);
    void Add_To_Material_Surface_Mesh_Subquadrilateral(const int curve_particle1,const int curve_particle2,const int triangle_particle1,const int triangle_particle2);
//#####################################################################
};
}
#endif
