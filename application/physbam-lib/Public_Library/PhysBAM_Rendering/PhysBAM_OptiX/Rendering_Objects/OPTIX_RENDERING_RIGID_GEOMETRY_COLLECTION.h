//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION__
#define __OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD_KEY_LISTENER.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD_MOUSE_LISTENER.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>

#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>

#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/SELECTION_RIGID_BODY.h>

#include <iostream>
#include <string>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Utilities/OPTIX_UTILITIES.h>

namespace PhysBAM{

using namespace optix;

template<class TV>
class OPTIX_SELECTION_RIGID_BODY : public SELECTION_RIGID_BODY<TV> {
    TV my_object_location;
    TV my_old_object_location;
    bool first;
public:
    int body_id;
    OPTIX_SELECTION_RIGID_BODY() : first(true) {
    }

    void Set_Object_Location(TV object_location) {
        my_object_location = object_location;
    }

    TV Get_Object_Location() {
        return my_object_location;
    }

    bool Was_Dragged() {
        bool result = first || (my_object_location != my_old_object_location);
        my_old_object_location = my_object_location;
        first = false;
        return result;
    }

    int Get_Body_ID() {
        return body_id;
    }
};

template<class T> class OPTIX_RENDERING_OBJECT;
template<class T> class OPTIX_RENDER_WORLD;

template<class T> class OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION : public OPTIX_RENDERING_GEOMETRY<T>, OPTIX_RENDER_WORLD_KEY_LISTENER, OPTIX_RENDER_WORLD_MOUSE_LISTENER {
    typedef VECTOR<T,3> TV;
    GENERIC_PARSER<T> *parser;
    Geometry rt_object;
    Program  rt_object_intersection_program;
    Program  rt_object_bounding_box_program;
    Acceleration rt_acceleration;
    OPTIX_RENDER_WORLD<T>* world_input;
    bool enable_fps_display,write_images;
    std::string image_path;

    // OPTIX_RENDERING_TRIANGULATED_SURFACE<T> **rendering_surfaces;
    ARRAY<OPTIX_RENDERING_TRIANGULATED_SURFACE<float> *, int> triangulated_surfaces;

    unsigned current_frame;
    unsigned start_frame, last_frame;
    Group rt_group;
    ARRAY<Transform> rt_transform_array;
    ARRAY<ARRAY<MATRIX<T, 4>, int>, int> animation_matrices;
    bool isPause;
    bool interactive;
    RIGID_GEOMETRY_COLLECTION<TV>* my_rigid_geometry_collection;
    ARRAY<OPTIX_SELECTION_RIGID_BODY<TV>, int> my_selections;

    void FetchOutTriangulatedSurfaces(RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection,
                                      ARRAY<OPTIX_MATERIAL<T>*, int> &materials) {
        //TV v = rigid_geometry_collection->particles.X(1);
        // std::cout << "rigid_geometry_collection->particles.X(1) = " << v.x << " " << v.y << " " << v.z << "\n";
        int number_of_geometry_particles = rigid_geometry_collection->particles.array_collection->Size();
        triangulated_surfaces.Resize(number_of_geometry_particles);
        if (interactive) {
            my_selections.Resize(number_of_geometry_particles);
        }

        /*OPTIX_MATERIAL<float> **materials (VECTOR<float,3>(0.2f, 0.2f, 0.2f),
                         VECTOR<float,3>(0.0f, 0.0f, 0.4f),
                         VECTOR<float,3>(0.2f, 0.2f, 0.2f),
                         VECTOR<float,3>(0.f, 0.f, 0.f),
                         VECTOR<float,3>(0.5f, 0.5f, 0.5f),
                         32.f);*/

        // std::cout << "Loaded Objects:\n";
        for(int id=1; id<=number_of_geometry_particles; id++) {
            RIGID_GEOMETRY<TV>& rigid_geometry = rigid_geometry_collection->Rigid_Geometry(id);
            // std::cout << rigid_geometry.Name() << "\n";
            // IMPLICIT_OBJECT<TV>* object_space_implicit_object=rigid_geometry.implicit_object->object_space_implicit_object;

            TRIANGULATED_SURFACE<T>* triangulated_surface=rigid_geometry.simplicial_object;
            /*
            int i,j,k;
            triangulated_surface->mesh.elements(1).Get(i,j,k);
            */
            // TV v1 = triangulated_surface->particles.X(i);
            // std::cout << "v = " << v1.x << ", " << v1.y << ", " << v1.z << "\n";


            triangulated_surfaces(id) = new OPTIX_RENDERING_TRIANGULATED_SURFACE<float>(triangulated_surface, materials(id % materials.m + 1));

            if (!interactive) {
                float angle;
                VECTOR<float,3> axis;
                rigid_geometry.Frame().r.Get_Angle_Axis(angle,axis);

                triangulated_surfaces(id)->Set_Transform(MATRIX<float,4>::Translation_Matrix(rigid_geometry.Frame().t)*
                                                         MATRIX<float,4>::Rotation_Matrix(axis, angle));
            } else {
                // my_selections(id) = new SELECTION_RIGID_BODY();
                my_selections(id).body_id = id;
                my_selections(id).Set_Object_Location(rigid_geometry_collection->particles.X(id));
                triangulated_surfaces(id)->Set_Transform(MATRIX<float,4>::Translation_Matrix(my_selections(id).Get_Object_Location()));
            }
        }
    }
public:
    // static rigid geometry
    OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION(RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection, ARRAY<OPTIX_MATERIAL<T>*, int> &materials):OPTIX_RENDERING_GEOMETRY<T>("Rigid Geometry Collection") {
        interactive = true;enable_fps_display=true;
        FetchOutTriangulatedSurfaces(rigid_geometry_collection, materials);
        my_rigid_geometry_collection = rigid_geometry_collection;
        current_frame = start_frame = last_frame = 0;
    }

    // rigid geometry animated from file
    OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION(std::string _basedir, ARRAY<OPTIX_MATERIAL<T>*, int> &materials, unsigned int frame_number=0, std::string _image_path="", bool isPaused=true):OPTIX_RENDERING_GEOMETRY<T>("Rigid Geometry Collection") {
        interactive = false;enable_fps_display=false;image_path=_image_path;
        if(_image_path.length()==0){write_images=false;}
        else{write_images=true;FILE_UTILITIES::Create_Directory(image_path);}
        std::string basedir = _basedir;
        isPause=isPaused;
        Initialize_Geometry_Particle();
        Initialize_Read_Write_General_Structures();

        FILE_UTILITIES::Read_From_Text_File(basedir+"/common/first_frame",start_frame);
        FILE_UTILITIES::Read_From_Text_File(basedir+"/common/last_frame",last_frame);
        // std::cout << "start_frame="<<start_frame<<"\n";
        // std::cout << "last_frame="<<last_frame<<"\n";
        current_frame = (frame_number==0 || !(start_frame<=frame_number && frame_number<=last_frame)) ? start_frame : frame_number;

        RIGID_BODY_COLLECTION<TV> *rigid_body_collection = new RIGID_BODY_COLLECTION<TV>(0, 0);
        RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection = &rigid_body_collection->rigid_geometry_collection;

        // reading geometry
        rigid_body_collection->Read(STREAM_TYPE(T()), basedir, current_frame);

        FetchOutTriangulatedSurfaces(rigid_geometry_collection, materials);
        my_rigid_geometry_collection = rigid_geometry_collection;
        //if a non-zero value was specified, ignore the previous frames
        start_frame = current_frame;        
        // reading all the matrices
        int number_of_geometry_particles = triangulated_surfaces.Size();
        animation_matrices.Resize(number_of_geometry_particles);
        for(int id=1; id<=number_of_geometry_particles; id++) {
            animation_matrices(id).Resize(last_frame - start_frame + 1);
        }
        rigid_body_collection = new RIGID_BODY_COLLECTION<TV>(0, 0);
        rigid_geometry_collection = &rigid_body_collection->rigid_geometry_collection;
        float eps = 1e-3f;
        VECTOR<float,3> scale_vec(1.f - eps, 1.f - eps, 1.f - eps);

        for (unsigned frame = start_frame; frame <= last_frame; frame++) {
            rigid_body_collection->Read(STREAM_TYPE(T()), basedir, frame);

            for(int id=1; id<=number_of_geometry_particles; id++) {
                RIGID_GEOMETRY<TV>& rigid_geometry = rigid_geometry_collection->Rigid_Geometry(id);

                float angle;VECTOR<float,3> axis;
                rigid_geometry.Frame().r.Get_Angle_Axis(angle,axis);
                animation_matrices(id)(frame + 1 - start_frame) = MATRIX<float,4>::Translation_Matrix(rigid_geometry.Frame().t)*MATRIX<float,4>::Rotation_Matrix(axis, angle) *
                                      MATRIX<float,4>::Scale_Matrix(scale_vec);
            }
        }
    }

    void Update(T time) {
        if (!interactive) {
            if (current_frame + 1 > last_frame || world_input->Check_If_Paused()) {
                return;
            }
            current_frame++;

            for(int id=1; id<=triangulated_surfaces.m; id++) {
                triangulated_surfaces(id)->Set_Transform(animation_matrices(id)(current_frame + 1 - start_frame));
                triangulated_surfaces(id)->Update(time);
            }
        } else {
            for(int id=1; id<=triangulated_surfaces.m; id++) {
                triangulated_surfaces(id)->Set_Transform(MATRIX<float,4>::Translation_Matrix(my_selections(id).Get_Object_Location()));
                triangulated_surfaces(id)->Update(time);
            }
        }
    }

    /*Geometry getGeometryObject() {
        return NULL;
    }*/

    void Initialize(OPTIX_RENDER_WORLD<T>* _world_input) {
        world_input = _world_input;
        rt_transform_array.Resize(triangulated_surfaces.m);

        for (int i = 1; i <= triangulated_surfaces.m; i++) {
            triangulated_surfaces(i)->Initialize(world_input);
            rt_transform_array(i) = triangulated_surfaces(i)->getTransform();
        }
        world_input->Add_Key_Listener(this);
        world_input->Add_Mouse_Listener(this);
        world_input->Enable_Fps_Display(enable_fps_display);
        world_input->Set_Image_Output(write_images,image_path);
        world_input->Pause(isPause);
    }

    // Group getGroup() { return rt_group; }
    ARRAY<Transform>& getTransforms() { return rt_transform_array; }

    void KeyPressed(unsigned char key, int x, int y) {
        // std::cout << "key pressed " << key << "\n";
        float offset_constant = 0.005f;
        // static clock_t start_time
        switch(key){
        case 'r':
            current_frame = start_frame;
            break;
        case 's':
            enable_fps_display=!enable_fps_display;
            world_input->Enable_Fps_Display(enable_fps_display);
            break;
        case 'w':
            std::cout<<"Image Writing: "<<(world_input->Toggle_Image_Writing()?"ON":"OFF")<<std::endl;
            break;
        case 'p':
            world_input->Pause(!(world_input->Check_If_Paused()));
            break;
        case 'b':
            if (my_selections.Size() > 0) {
                // rigid_geometry_collection->Rigid_Geometry(1).V() = offset_constant / delta_time;
                // ((OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<float> *)obj)->rigid_geometry_collection_simulation->Rigid_Geometry(real_selection->body_id).V()=object_velocity;
                my_selections(1).Set_Object_Location(my_selections(1).Get_Object_Location() + TV(-offset_constant, 0, 0));}
            break;
        case 'n':
            if (my_selections.Size() > 0)
                my_selections(1).Set_Object_Location(my_selections(1).Get_Object_Location() + TV(offset_constant, 0, 0));
            break;
        case 'j':
            if (my_selections.Size() > 0)
                my_selections(1).Set_Object_Location(my_selections(1).Get_Object_Location() + TV(0, offset_constant, 0));
            break;
        case 'm':
            if (my_selections.Size() > 0)
                my_selections(1).Set_Object_Location(my_selections(1).Get_Object_Location() + TV(0, -offset_constant, 0));
            break;
        default:
            break;
        }
    }

    bool draggingMode;
    int x_old, y_old;
    TV selection_world_point;
    clock_t old_time;

    void Handle_Click(int button, int state, int x, int y, bool ctrl_pressed, bool shift_pressed) {
        // std::cout << "Handle click\n";
        // std::cout << "button == GLUT_LEFT_BUTTON === " << (button == GLUT_LEFT_BUTTON) << "\n";
        if (button == GLUT_LEFT_BUTTON) {
            draggingMode = (state == GLUT_DOWN) && ctrl_pressed;
            // std::cout << "ctrl_pressed = " << ctrl_pressed << "(state == GLUT_DOWN) " << (state == GLUT_DOWN) << "\n";
            // std::cout << "draggingMode = " << draggingMode << "\n";
            x_old = x;
            y_old = y;

            selection_world_point = my_rigid_geometry_collection->Rigid_Geometry(1).X();
            // std::cout << "world_point = " << selection_world_point.x << " " << selection_world_point.y << " " << selection_world_point.z << "\n";
            old_time = clock();
        } else {
            draggingMode = false;
        }
        if (!draggingMode) {
            my_rigid_geometry_collection->Rigid_Geometry(1).V()=TV(0.f, 0.f, 0.f);
        }
        // std::cout << "Handle_Click\n";
    }

    void Handle_Drag(int x, int y) {
        if (draggingMode) {
            // std::cout << "x - x_old " << x - x_old << "\n";
            // std::cout << "y - y_old " << y - y_old << "\n";
            TV delta = world_input->getWorldDeltaFromScreenDelta(x - x_old, y - y_old, selection_world_point);
            selection_world_point += delta;
            if (my_selections.Size() > 0) {
                my_selections(1).Set_Object_Location(my_selections(1).Get_Object_Location() + delta);
            }
            x_old = x;
            y_old = y;
            clock_t current_time = clock();
            if (current_time != old_time)
                my_rigid_geometry_collection->Rigid_Geometry(1).V()=delta / ((current_time - old_time) / (float)CLOCKS_PER_SEC);
            old_time = current_time;
        }
        // std::cout << "Handle_Drag\n";
    }

    OPTIX_SELECTION_RIGID_BODY<TV>* getCurrentSelection() {
        if (my_selections.Size() > 0) {
            return &my_selections(1);
        }
        return NULL;
    }
};
}
#endif
#endif

