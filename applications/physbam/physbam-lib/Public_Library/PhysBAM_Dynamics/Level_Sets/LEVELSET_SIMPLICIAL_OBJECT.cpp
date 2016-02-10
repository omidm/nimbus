#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_SIMPLICIAL_OBJECT
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Level_Sets/READ_WRITE_LEVELSET_RED_GREEN.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_SIMPLICIAL_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
LEVELSET_SIMPLICIAL_OBJECT(EMBEDDED_MATERIAL_SURFACE<TV,d>& embedding_input)
    :embedding(embedding_input),particles(embedding.particles),levelset(grid,phi),need_destroy_embedding(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
~LEVELSET_SIMPLICIAL_OBJECT()
{
    if(need_destroy_embedding) delete &embedding;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> LEVELSET_SIMPLICIAL_OBJECT<TV,d>* LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
Create()
{
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    LEVELSET_SIMPLICIAL_OBJECT* levelset=new LEVELSET_SIMPLICIAL_OBJECT(*EMBEDDED_MATERIAL_SURFACE<TV,d>::Create());
    levelset->need_destroy_embedding=true;return levelset;
#endif
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int d> void LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
Initialize()
{
    phi.Resize(grid.number_of_nodes);V.Resize(grid.number_of_nodes);
}
//#####################################################################
// Function Build_Embedded_Object
//#####################################################################
template<class TV,int d> void LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
Build_Embedded_Object(const bool verbose)
{
    if(verbose) LOG::cout<<"Building embedded object"<<std::endl;
    embedding.embedded_object.Clean_Memory();
    T_MESH& mesh=embedding.embedded_object.simplicial_object.mesh;
    grid.Build_Mesh(mesh,&phi,cell_to_simplex_mapping,&node_to_particle_mapping);
    // extend particle subset if necessary
    if(particles.array_collection->Size()<mesh.number_nodes) particles.Add_Elements(mesh.number_nodes-particles.array_collection->Size());
    // fix mesh to point at global particles
    for(int t=1;t<=mesh.elements.m;t++) mesh.elements(t)=VECTOR<int,d+1>::Map(particles.active_indices,mesh.elements(t));
    mesh.number_nodes=particles.point_cloud.array_collection->Size();
    // remap phi
    Initialize();
    ARRAY<T> particle_based_phi(particles.array_collection->Size());
    for(int i=1;i<=grid.number_of_nodes;i++)if(node_to_particle_mapping(i))particle_based_phi(node_to_particle_mapping(i))=phi(i);
    if(verbose){
        LOG::cout<<"  particles: "<<particles.array_collection->Size()<<std::endl;
        LOG::cout<<"  elements: "<<mesh.elements.m<<std::endl;}
    embedding.embedded_object.Calculate_Boundary_From_Levelset_On_Nodes(particle_based_phi,false,true);
}
//#####################################################################
// Function Get_Material_Coordinates
//#####################################################################
template<class TV,int d> ARRAY<TV> LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
Get_Material_Coordinates()
{
    ARRAY<TV> Xm(particles.array_collection->Size());
    ARRAY<VECTOR<T,d> >& cell_based_Xm=grid.Node_Locations();
    for(int i=1;i<=grid.number_of_nodes;i++)if(node_to_particle_mapping(i)) Xm(node_to_particle_mapping(i))=TV(cell_based_Xm(i));
    return Xm;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV,int d> void LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
Euler_Step(const T& dt,const T& time)
{
    V.Resize(grid.number_of_nodes);levelset.Euler_Step(V,dt,time);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,int d> void LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
Read(TYPED_ISTREAM& input)
{
    Read_Binary(input,grid,cell_to_simplex_mapping,node_to_particle_mapping,phi);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,int d> void LEVELSET_SIMPLICIAL_OBJECT<TV,d>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,grid,cell_to_simplex_mapping,node_to_particle_mapping,phi);
}
//#####################################################################
template class LEVELSET_SIMPLICIAL_OBJECT<VECTOR<float,2>,2>;
template class LEVELSET_SIMPLICIAL_OBJECT<VECTOR<float,3>,2>;
template class LEVELSET_SIMPLICIAL_OBJECT<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_SIMPLICIAL_OBJECT<VECTOR<double,2>,2>;
template class LEVELSET_SIMPLICIAL_OBJECT<VECTOR<double,3>,2>;
template class LEVELSET_SIMPLICIAL_OBJECT<VECTOR<double,3>,3>;
#endif
}
#endif
