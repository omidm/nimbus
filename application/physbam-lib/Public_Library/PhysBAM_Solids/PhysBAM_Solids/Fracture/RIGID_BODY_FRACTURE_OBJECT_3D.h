//#####################################################################
// Copyright 2006, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_FRACTURE_OBJECT_3D
//#####################################################################
#ifndef __RIGID_BODY_FRACTURE_OBJECT_3D__
#define __RIGID_BODY_FRACTURE_OBJECT_3D__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_GRAIN_BOUNDARIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/LEVELSET_GRAIN_BOUNDARIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/PAINTED_GRAIN_BOUNDARIES.h>
namespace PhysBAM{

template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class SOLID_BODY_COLLECTION;

template<class T_input>
class RIGID_BODY_FRACTURE_OBJECT_3D:public RIGID_BODY<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef RIGID_BODY<TV> BASE;
    using BASE::Angular_Momentum;using BASE::simplicial_object;using BASE::implicit_object;using BASE::thin_shell;using BASE::particle_index;using BASE::axis_aligned_bounding_box;
    using BASE::Frame;using BASE::Twist;using BASE::Mass;using BASE::Inertia_Tensor;using BASE::Update_Angular_Velocity;using BASE::CFL_initialized;using BASE::Initialize_CFL;using BASE::bounding_box_radius;

    PARTICLES<TV> particles; // object space particles
    SOLID_BODY_COLLECTION<TV>& solid_body_collection; // this is the global deformable object
    ARRAY<int>& particle_to_rigid_body_id;
    ARRAY<int>& deformable_to_rigid_particles; // which rigid_body_particles each deformable particle matches to (internal fractures may cause many to 1 mapping)
    ARRAY<int> rigid_to_deformable_particles; // these number should be the id numbers in the deformable particles
    ARRAY<int> rigid_to_deformable_tets; // these number should be the id numbers in the deformable particles
    int parent_rigid_body_id;
    ARRAY<TV> average_dX;
    T fracture_threshold;
    ARRAY<FRACTURE_GRAIN_BOUNDARIES<TV,3>*> grain_boundaries;
    ARRAY<int> new_vertices_to_parent_vertices;
    ARRAY<int> new_triangles;
    ARRAY<int> old_triangles;
    ARRAY<int> partial_triangles;

    RIGID_BODY_FRACTURE_OBJECT_3D(PARTICLES<TV>& deformable_body_particles,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,
        ARRAY<int>& particle_to_rigid_body_id_input,ARRAY<int>& deformable_to_rigid_particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
        :RIGID_BODY<TV>(rigid_body_collection_input,true),solid_body_collection(solid_body_collection_input),
        particle_to_rigid_body_id(particle_to_rigid_body_id_input),deformable_to_rigid_particles(deformable_to_rigid_particles_input),parent_rigid_body_id(0),fracture_threshold(0)
    {}
    ~RIGID_BODY_FRACTURE_OBJECT_3D(){}

    void Update_Particle_To_Rigid_Body_Id_Mapping()
    {for(int p=1;p<=rigid_to_deformable_particles.m;p++) particle_to_rigid_body_id(rigid_to_deformable_particles(p))=particle_index;}

    void Initialize_Grain_Boundaries(const ARRAY<TV>& seed_positions,const ARRAY<T>& seed_weakness_multipliers,const FRACTURE_CALLBACKS<TV>* fracture_callbacks,const bool levelset_grain_boundaries=true) // seed positions are specified in world space
    {
        ARRAY<TV> object_space_seed_positions(seed_positions.m);
        SIMPLEX_MESH<3>& mesh=this->template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh;
        FRAME<TV> inverse_frame=Frame().Inverse();
        for(int i=1;i<=seed_positions.m;i++) object_space_seed_positions(i)=inverse_frame*seed_positions(i);
        if(levelset_grain_boundaries){
            grain_boundaries.Append(new LEVELSET_GRAIN_BOUNDARIES<TV,3>(particles,mesh,object_space_seed_positions,seed_weakness_multipliers,Frame(),fracture_callbacks));
        }
        else{
            //ARRAY<TV>* seed_positions=new ARRAY<TV>();ARRAY<T>* seed_weaknesses=new ARRAY<T>();
            {std::stringstream ss;ss<<"creating painted boundaries"<<std::endl;LOG::filecout(ss.str());}
            grain_boundaries.Append(new PAINTED_GRAIN_BOUNDARIES<TV,3>(particles,mesh,object_space_seed_positions,seed_weakness_multipliers,Frame(),fracture_callbacks));}
        {std::stringstream ss;ss<<"done with init grain boundaries"<<std::endl;LOG::filecout(ss.str());}
        // put bool for type of grain boundaries?
    }

    void Initialize_Rigid_Body_From_Fragment(const T density,const T cell_size,const int subdivision_loops,const bool use_implicit_surface_maker,const int levels_of_octree)
    {
        EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& embedding=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>&>();
        TETRAHEDRALIZED_VOLUME<T>& deformable_tetrahedralized_volume=embedding.embedded_object.simplicial_object;
        if(!embedding.material_surface_mesh.elements.m) embedding.Create_Material_Surface();
        assert(particles.array_collection->Size()==0&&rigid_to_deformable_particles.m==0);

        particles.array_collection->Add_Elements(solid_body_collection.deformable_body_collection.dynamic_particles.m);
        rigid_to_deformable_particles.Resize(solid_body_collection.deformable_body_collection.dynamic_particles.m);

        particle_to_rigid_body_id.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        deformable_to_rigid_particles.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size());

        // get all the particles of the material surface
        int index=1;
        for(int particle_id=1;particle_id<=solid_body_collection.deformable_body_collection.dynamic_particles.m;particle_id++){
            int i=solid_body_collection.deformable_body_collection.dynamic_particles(particle_id);
            rigid_to_deformable_particles(index)=i;
            deformable_to_rigid_particles(i)=index;
            particles.X(index)=solid_body_collection.deformable_body_collection.particles.X(i);
            index++;}

        ARRAY<VECTOR<int,3> > tri_elements;tri_elements.Preallocate(embedding.material_surface.mesh.elements.m);
        for(int t=1;t<=embedding.material_surface_mesh.elements.m;t++){
            int index1=embedding.material_surface_mesh.elements(t)[1];
            int index2=embedding.material_surface_mesh.elements(t)[2],index3=embedding.material_surface_mesh.elements(t)[3];
            tri_elements.Append(VECTOR<int,3>(deformable_to_rigid_particles(index1),deformable_to_rigid_particles(index2),deformable_to_rigid_particles(index3)));}

        ARRAY<VECTOR<int,4> > tet_elements;tet_elements.Preallocate(deformable_tetrahedralized_volume.mesh.elements.m);rigid_to_deformable_tets.Preallocate(deformable_tetrahedralized_volume.mesh.elements.m);
        for(int t=1;t<=deformable_tetrahedralized_volume.mesh.elements.m;t++){
            int index1=deformable_tetrahedralized_volume.mesh.elements(t)[1];
            int index2=deformable_tetrahedralized_volume.mesh.elements(t)[2],index3=deformable_tetrahedralized_volume.mesh.elements(t)[3],index4=deformable_tetrahedralized_volume.mesh.elements(t)[4];
            tet_elements.Append(VECTOR<int,4>(deformable_to_rigid_particles(index1),deformable_to_rigid_particles(index2),deformable_to_rigid_particles(index3),deformable_to_rigid_particles(index4)));
            rigid_to_deformable_tets.Append(t);}

        TRIANGULATED_SURFACE<T>* material_surface=new TRIANGULATED_SURFACE<T>(*new TRIANGLE_MESH(particles.array_collection->Size(),tri_elements),particles);
        material_surface->Update_Number_Nodes();
        TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=new TETRAHEDRALIZED_VOLUME<T>(*new TETRAHEDRON_MESH(particles.array_collection->Size(),tet_elements),particles);
        tetrahedralized_volume->Update_Number_Nodes();
        tetrahedralized_volume->Initialize_Hierarchy();

        if(thin_shell) PHYSBAM_FATAL_ERROR("thin shell case currently not handled");
        MASS_PROPERTIES<TV> mass_properties(*material_surface,false); // TODO: should this handle the thin shell case as well?
        mass_properties.Set_Density(density);Mass()=mass_properties.Mass();
        FRAME<TV> frame_local;mass_properties.Transform_To_Object_Frame(frame_local,Inertia_Tensor());Set_Frame(frame_local);
        FRAME<TV> inverse_frame=Frame().Inverse();
        for(int i=1;i<=particles.array_collection->Size();i++) particles.X(i)=inverse_frame*particles.X(i);

        RANGE<TV> box=RANGE<TV>::Empty_Box();for(int i=1; i <= particles.array_collection->Size(); i++) box.Enlarge_To_Include_Point(particles.X(i));
        T uniform_levelset_cell_size=box.Edge_Lengths().Max_Abs()/30; // 30 chosen based on average, don't think this actually used so wanted something reasonable        {std::stringstream ss;ss<<"Bounding box: "<<box<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"cell_size: "<<uniform_levelset_cell_size<<" and edge length: "<<box.Edge_Lengths().Max_Abs()<<std::endl;LOG::filecout(ss.str());}

        // pass false for create_levelset_test since we always want to use the material surface which is more detailed
        Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface(*tetrahedralized_volume,*material_surface,uniform_levelset_cell_size,subdivision_loops,0,use_implicit_surface_maker,levels_of_octree);
    }

    void Align_Deformable_Object_With_Rigid_Body()
    {
        // currently contains all particles
        PARTICLES<TV>& deformable_particles=solid_body_collection.deformable_body_collection.particles;
        for(int p=1;p<=rigid_to_deformable_particles.m;p++) {
            deformable_particles.array_collection->Copy_Element(*particles.array_collection,deformable_to_rigid_particles(rigid_to_deformable_particles(p)),rigid_to_deformable_particles(p));
            deformable_particles.X(rigid_to_deformable_particles(p))=Frame()*particles.X(deformable_to_rigid_particles(rigid_to_deformable_particles(p)));
            deformable_particles.V(rigid_to_deformable_particles(p))=Pointwise_Object_Velocity(deformable_particles.X(deformable_to_rigid_particles(rigid_to_deformable_particles(p))));}
        solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    }

    void Populate_Triangle_Lists(RIGID_BODY_FRACTURE_OBJECT_3D<T>* parent_rigid_body_fracture_object)
    {}

//#####################################################################
};
}
#endif
