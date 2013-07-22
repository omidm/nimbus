//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_ANIMATED_VISUALIZATION
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_ANIMATED_VISUALIZATION__
#define __OPTIX_ANIMATED_VISUALIZATION__

#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_ANIMATED_VISUALIZATION_INTERFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_SMOKE.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL_PHONG.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW_GLUT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/SELECTION_RIGID_BODY.h>

#include <climits>

namespace PhysBAM
{

// for interactive visualization, assume we have smoke object and just ONE rigid body in the form of triangulated surface
template<class T>
class OPTIX_ANIMATED_VISUALIZATION : public OPTIX_ANIMATED_VISUALIZATION_INTERFACE<VECTOR<T,3> > {
private:
    typedef VECTOR<T,3> TV;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    ARRAY<T,TV_INT>* my_densities_simulated;
    RIGID_GEOMETRY_COLLECTION<TV>* my_rigid_geometry_collection_simulated;
    ARRAY<T,FACE_INDEX<TV::dimension> >* my_face_velocities_simulated;
    GRID<TV> my_grid;

    int current_frame;
    OPTIX_RENDERING_SMOKE<T> *rendering_smoke;
    OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION<T> *rigid_geometry;

public:
    GRID<TV>& Get_Grid() { return my_grid;}

    // void Set_Scalar_Values_Simulated(ARRAY<T,TV_INT>* densities_simulated = NULL) {return rendering_smoke->Set_Scalar_Values_Simulated(densities_simulated);}
    // ARRAY<T,TV_INT>* Get_Scalar_Values_Simulated() {return my_densities_simulated; }

    void Set_Scalar_Values_Upsample_Scale(int upsample_scale) {}
    int Get_Upsample_Scale() {return 1;}

    // void Set_Rigid_Bodies_Simulated(RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection_simulated) {return rendering_smoke->Set_Rigid_Bodies_Simulated(rigid_geometry_collection_simulated); }
    RIGID_GEOMETRY_COLLECTION<TV>* Get_Rigid_Bodies_Simulated() {return my_rigid_geometry_collection_simulated;}

    OPTIX_ANIMATED_VISUALIZATION() : current_frame(0) {
        OPTIX_RENDER_WORLD<float>::Instance()->window = new OPENGL_WINDOW_GLUT(*OPTIX_RENDER_WORLD<float>::Instance(), "PhysBAM OptiX",
                           OPTIX_RENDER_WORLD<float>::Instance()->camera->screen_width,
                           OPTIX_RENDER_WORLD<float>::Instance()->camera->screen_height);

        rendering_smoke = new OPTIX_RENDERING_SMOKE<T>();
        OPTIX_RENDER_WORLD<T>::Instance()->Add_Object(rendering_smoke);
        my_densities_simulated = new ARRAY<T,TV_INT>();
    }
    ~OPTIX_ANIMATED_VISUALIZATION() {
        delete rendering_smoke;
    }

    void Initialize() {
        OPTIX_RENDER_WORLD<T>::Instance()->Initialize();
    }

    void Set_Frame(int frame_input) {
        current_frame = frame_input;
        // std::cout << "Frame: " << current_frame << "\n";

        /*
        for(typename GRID<TV>::CELL_ITERATOR iterator(my_grid);iterator.Valid();iterator.Next())
            std::cout << (*my_densities_simulated)(iterator.Cell_Index()) << " ";

            std::cout << "\n\n***********************\n\n";
            */
        rendering_smoke->ReinitializeFromSimulation();

         /*       rigid_geometry_collection->particles.Resize(rigid_geometry_collection_simulation->particles.array_collection->Size());
        for(int i=1;i<=rigid_geometry_collection->particles.array_collection->Size();i++) {
            rigid_geometry_collection->particles.X(i)=rigid_geometry_collection_simulation->particles.X(i);
            rigid_geometry_collection->particles.rigid_geometry(i)=rigid_geometry_collection_simulation->particles.rigid_geometry(i);
            rigid_geometry_collection->particles.rotation(i)=rigid_geometry_collection_simulation->particles.rotation(i);
            rigid_geometry_collection->particles.V(i)=rigid_geometry_collection_simulation->particles.V(i);
            rigid_geometry_collection->particles.V(i)=rigid_geometry_collection_simulation->particles.V(i);
            rigid_geometry_collection->particles.angular_velocity(i)=rigid_geometry_collection_simulation->particles.angular_velocity(i);
            rigid_geometry_collection->particles.structure_ids(i)=rigid_geometry_collection_simulation->particles.structure_ids(i);}

        int max_number_of_bodies=max(opengl_triangulated_surface.Size(),rigid_geometry_collection->particles.array_collection->Size());
        opengl_colors.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(opengl_colors,OPENGL_COLOR::Cyan());
        Resize_Structures(max_number_of_bodies);

        for(int i=1;i<=max_number_of_bodies;i++) if(rigid_geometry_collection->Is_Active(i)) Create_Geometry(i);

        // Update active bodies / remove inactive bodies
        for(int id(1);id<=rigid_geometry_collection->particles.array_collection->Size();id++){
            if(rigid_geometry_collection->Is_Active(id)){
                Update_Geometry(id);
                RIGID_GEOMETRY<TV>& body=rigid_geometry_collection->Rigid_Geometry(id);
                IMPLICIT_OBJECT<TV>* object_space_implicit_object=body.implicit_object?body.implicit_object->object_space_implicit_object:0;
                if(body.name=="ground" || (object_space_implicit_object && typeid(*object_space_implicit_object)==typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >)))
                    Set_Object_Material(id,OPENGL_COLOR::Ground_Tan((T)1));
                else{
                    if(one_sided) Set_Object_Material(id,front_color_map->Lookup(Value(id)-1));
                    else Set_Object_Material(id,front_color_map->Lookup(Value(id)-1),back_color_map->Lookup(Value(id)-1));}}
            else Destroy_Geometry(id);}
        for(int id=rigid_geometry_collection->particles.array_collection->Size()+1;id<=opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);
        */
    }

    void Set_MAC_Velocities_Simulated(ARRAY<T,FACE_INDEX<TV::dimension> >* face_velocities_simulated) {
        my_face_velocities_simulated = face_velocities_simulated;
        // std::cout << "Velocities copied\n";
    }

    void Set_Grid(GRID<TV> grid) {
        my_grid = grid;
    }

    void Set_Scalar_Values_Simulated(ARRAY<T,TV_INT>* densities_simulated = NULL) {
        if (densities_simulated != NULL) {
            if (my_densities_simulated != NULL) {
                delete my_densities_simulated;
            }
            my_densities_simulated = densities_simulated;
        } else {
            my_densities_simulated->Resize(my_grid.Domain_Indices());
        }
        rendering_smoke->Set_Densities(my_densities_simulated, &my_grid);
        // std::cout << "Scalar values copied\n";
        // OOPS! my_grid might not be initialized
    }

    ARRAY<T,TV_INT>* Get_Scalar_Values_Simulated() {
        return my_densities_simulated;
    }

    void Set_Rigid_Bodies_Simulated(RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection_simulated) {
        my_rigid_geometry_collection_simulated = rigid_geometry_collection_simulated;
        ARRAY<OPTIX_MATERIAL<T> *, int> materials;
        materials.Resize(1);
        materials(1) = new OPTIX_MATERIAL_PHONG<T>(VECTOR<float,3>(0.1f, 0.1f, 0.1f),
                             VECTOR<float,3>(0.7f, 0.4f, 0.f),
                             VECTOR<float,3>(0.1f, 0.1f, 0.1f),
                             VECTOR<float,3>(0.f, 0.f, 0.f),
                             VECTOR<float,3>(0.5f, 0.5f, 0.5f),
                             32.f);

        rigid_geometry = new OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION<T>(rigid_geometry_collection_simulated, materials);
        OPTIX_RENDER_WORLD<T>::Instance()->Add_Object(rigid_geometry);
        // std::cout << "Rigid bodies copied\n";
    }

    void Run() {
        OPTIX_RENDER_WORLD<float>::Instance()->Render_Frame();
    }

    bool Has_Current_Selection() { return my_rigid_geometry_collection_simulated->particles.array_collection->Size() != 0; }
    SELECTION_RIGID_BODY<TV>* Get_Current_Selection() {return rigid_geometry->getCurrentSelection();}

    void Rigid_Geometry_Updated() {}
};

}

#endif
#endif
