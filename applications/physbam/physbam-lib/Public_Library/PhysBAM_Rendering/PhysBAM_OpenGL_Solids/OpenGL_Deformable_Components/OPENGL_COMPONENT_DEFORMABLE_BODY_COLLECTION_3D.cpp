//#####################################################################
// Copyright 2004-2009, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
//#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FREE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Deformable_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D.h>


#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D()
    :OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>(false),
    display_hard_bound_surface_mode(0),display_forces_mode(0),interaction_pair_display_mode(0),
    deformable_body_collection(0,collision_body_list),
    has_embedded_objects(false),has_soft_bindings(false)
{
    deformable_geometry=&deformable_body_collection.deformable_geometry;
    color_map_forces=new OPENGL_COLOR_RAMP<T>();
    color_map_forces->Add_Color((T)-1.,OPENGL_COLOR(0,0,1));
    color_map_forces->Add_Color((T)-.01,OPENGL_COLOR(0,1,1));
    color_map_forces->Add_Color((T).0,OPENGL_COLOR(0,0,0));
    color_map_forces->Add_Color((T).01,OPENGL_COLOR(1,1,0));
    color_map_forces->Add_Color((T)1.,OPENGL_COLOR(1,0,0));
}
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D(const std::string& prefix,const int start_frame)
    :OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>(prefix,start_frame,false),
    display_hard_bound_surface_mode(0),display_forces_mode(0),interaction_pair_display_mode(0),
    deformable_body_collection(0,collision_body_list),
    has_embedded_objects(false),has_soft_bindings(false)
{
    deformable_geometry=&deformable_body_collection.deformable_geometry;
    color_map_forces=new OPENGL_COLOR_RAMP<T>();
    color_map_forces->Add_Color((T)-1.,OPENGL_COLOR(0,0,1));
    color_map_forces->Add_Color((T)-.01,OPENGL_COLOR(0,1,1));
    color_map_forces->Add_Color((T).0,OPENGL_COLOR(0,0,0));
    color_map_forces->Add_Color((T).01,OPENGL_COLOR(1,1,0));
    color_map_forces->Add_Color((T)1.,OPENGL_COLOR(1,0,0));
}
//#####################################################################
// Function ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D()
{
    delete color_map_forces;
}
//#####################################################################
// Function Create_Hard_Bound_Boundary_Surface
//#####################################################################
template<class T,class RW> TRIANGULATED_SURFACE<T>& OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Create_Hard_Bound_Boundary_Surface(TRIANGULATED_SURFACE<T>& boundary_surface)
{
    TRIANGULATED_SURFACE<T>& hard_bound_boundary_surface=*TRIANGULATED_SURFACE<T>::Create(boundary_surface.particles);
    ARRAY<int> particle_map(IDENTITY_ARRAY<>(boundary_surface.particles.array_collection->Size()));
#if 0 // TODO: Fix me
    for(int b=1;b<=deformable_body_collection.soft_bindings.bindings.m;b++){VECTOR<int,2>& binding=deformable_body_collection.soft_bindings.bindings(b);particle_map(binding.x)=binding.y;}
#endif
    for(int t=1;t<=boundary_surface.mesh.elements.m;t++) hard_bound_boundary_surface.mesh.elements.Append(VECTOR<int,3>::Map(particle_map,boundary_surface.mesh.elements(t)));
    return hard_bound_boundary_surface;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Reinitialize(bool force,bool read_geometry)
{
    //std::cout<<"M: Tri Before reinit "<<dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection_simulated->deformable_geometry.structures(1))->Get_Segment_Mesh().elements.m<<std::endl;
    //std::cout<<"M: Before reinit "<<dynamic_cast<LINEAR_SPRINGS<TV>*>(deformable_body_collection_simulated->deformables_forces(2))->segment_mesh.elements.m<<std::endl;
    if(is_interactive){BASE::Reinitialize_From_Simulation(force,read_geometry);
        //std::cout<<"M: Tri After reinit "<<dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection_simulated->deformable_geometry.structures(1))->Get_Segment_Mesh().elements.m<<std::endl;
        //std::cout<<"M: After reinit "<<dynamic_cast<LINEAR_SPRINGS<TV>*>(deformable_body_collection_simulated->deformables_forces(2))->segment_mesh.elements.m<<std::endl;
        return;}    
    if(!(draw && (force || (is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)))) return;
    static bool first_time=true;
    std::string frame_string=STRING_UTILITIES::string_sprintf("%s/%d/",prefix.c_str(),frame);
    std::string static_frame_string=frame_string;
    int static_frame=FILE_UTILITIES::File_Exists(frame_string+"deformable_object_structures")?frame:-1;
    bool read_static_variables=static_frame!=-1 || !deformable_body_collection.deformable_geometry.structures.m;
    deformable_body_collection.Read(STREAM_TYPE(RW()),prefix,prefix,frame,static_frame,read_static_variables,true);

    BASE::Reinitialize(force,false);

    if(FILE_UTILITIES::File_Exists(frame_string+"/interaction_pairs") && interaction_pair_display_mode)
        FILE_UTILITIES::Read_From_File(STREAM_TYPE(RW()),frame_string+"/interaction_pairs",point_triangle_interaction_pairs,edge_edge_interaction_pairs);
    else{
        point_triangle_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();}
    
    std::string filename=frame_string+"/deformable_object_force_data";
    if(FILE_UTILITIES::File_Exists(filename)){
        if(first_time) {std::stringstream ss;ss<<"reading "<<filename<<std::endl;LOG::filecout(ss.str());}
        FILE_UTILITIES::Read_From_File(STREAM_TYPE(RW()),filename,force_data_list);}
    else force_data_list.Remove_All();

    if(read_static_variables){
        int m=deformable_body_collection.deformable_geometry.structures.m;active_list.Resize(m,true,true,true);
        embedded_surface_objects.Delete_Pointers_And_Clean_Memory();embedded_surface_objects.Resize(m);
        boundary_surface_objects.Delete_Pointers_And_Clean_Memory();boundary_surface_objects.Resize(m);
        hard_bound_boundary_surface_objects.Delete_Pointers_And_Clean_Memory();hard_bound_boundary_surface_objects.Resize(m);
        for(int i=1;i<=m;i++){
            STRUCTURE<TV>* structure=deformable_body_collection.deformable_geometry.structures(i);
            if(EMBEDDED_MATERIAL_SURFACE<TV,2>* embedding=dynamic_cast<EMBEDDED_MATERIAL_SURFACE<TV,2>*>(structure)){
                if(first_time) {std::stringstream ss;ss<<"object "<<i<<": embedded triangulated surface\n";LOG::filecout(ss.str());}
                boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(embedding->material_surface,false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Yellow()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Cyan()));
                hard_bound_boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(Create_Hard_Bound_Boundary_Surface(embedding->material_surface),false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Magenta(.5f)));
                embedding->embedded_object.simplicial_object.mesh.Initialize_Segment_Mesh();
                triangulated_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(embedding->embedded_object.simplicial_object,false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()));}
            else if(EMBEDDED_MATERIAL_SURFACE<TV,3>* embedding=dynamic_cast<EMBEDDED_MATERIAL_SURFACE<TV,3>*>(structure)){
                if(first_time) {std::stringstream ss;ss<<"object "<<i<<": embedded tetrahedralized volume\n";LOG::filecout(ss.str());}
                boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(embedding->material_surface,false);
                embedded_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(embedding->embedded_object.embedded_object,false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()));
                hard_bound_boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(Create_Hard_Bound_Boundary_Surface(embedding->material_surface),false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Magenta(.5f)));
                embedding->embedded_object.simplicial_object.mesh.Initialize_Neighbor_Nodes();
                tetrahedralized_volume_objects(i)=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(&embedding->embedded_object.simplicial_object.mesh,
                    &deformable_body_collection.particles,OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()));}
            else if(EMBEDDING<TV>* embedding=dynamic_cast<EMBEDDING<TV>*>(structure)){
                if(first_time) {std::stringstream ss;ss<<"object "<<i<<": embedding\n";LOG::filecout(ss.str());}
                boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(embedding->material_surface,false);
                hard_bound_boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(Create_Hard_Bound_Boundary_Surface(embedding->material_surface),false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Magenta(.5f)));}
            else{if(first_time) {std::stringstream ss;ss<<"object "<<i<<": object unrecognized at body level\n";LOG::filecout(ss.str());}}}}
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++){
#if 0 // TODO: Fix me
        if(deformable_body_collection.soft_bindings.bindings.m) has_soft_bindings=true;
#endif
        if(boundary_surface_objects(i)) has_embedded_objects=true;}
    if(smooth_shading){
        for(int i=1;i<=boundary_surface_objects.m;i++) if(boundary_surface_objects(i))boundary_surface_objects(i)->Initialize_Vertex_Normals();
        for(int i=1;i<=embedded_surface_objects.m;i++) if(embedded_surface_objects(i))embedded_surface_objects(i)->Initialize_Vertex_Normals();}

    first_time=false;
}
//#####################################################################
// Function Set_Display_Modes_For_Geometry_Collection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Set_Display_Modes_For_Geometry_Collection(bool& display_triangulated_surface_objects,bool& display_tetrahedralized_volume_objects,
        bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects) const
{
    if(has_embedded_objects || has_soft_bindings){
        display_triangulated_surface_objects=display_mode!=1 && display_soft_bound_surface_mode!=2;
        display_tetrahedralized_volume_objects=display_mode!=1 && display_mode!=3;// && display_soft_bound_surface_mode;
        display_hexahedralized_volume_objects=display_mode!=1;
        display_free_particles_objects=display_mode!=1;}
    else{
        display_triangulated_surface_objects=display_mode!=1;
        display_tetrahedralized_volume_objects=display_mode!=2;
        display_hexahedralized_volume_objects=display_mode!=2;
        display_free_particles_objects=display_mode!=1;}
}
//#####################################################################
// Function Set_Display_Modes
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Set_Display_Modes(bool& display_triangulated_surface_objects,bool& display_tetrahedralized_volume_objects,
        bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects,bool& display_boundary_surface_objects,
        bool& display_hard_bound_boundary_surface_objects) const
{
    // TODO: fix display modes
    Set_Display_Modes_For_Geometry_Collection(display_triangulated_surface_objects,display_tetrahedralized_volume_objects,
            display_hexahedralized_volume_objects,display_free_particles_objects);

    if(has_embedded_objects || has_soft_bindings){
        display_boundary_surface_objects=display_mode!=2 && display_hard_bound_surface_mode!=2;
        display_hard_bound_boundary_surface_objects=display_mode!=2 && display_hard_bound_surface_mode;}
    else{
        display_boundary_surface_objects=display_mode!=1;
        display_hard_bound_boundary_surface_objects=display_mode!=1 && display_hard_bound_surface_mode;}
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Display(const int in_color) const
{
    BASE::Display(in_color);

    if(!draw || !valid) return;
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}
    bool display_triangulated_surface_objects,display_tetrahedralized_volume_objects,display_hexahedralized_volume_objects,
         display_boundary_surface_objects,display_hard_bound_boundary_surface_objects,display_free_particles_objects;
    Set_Display_Modes(display_triangulated_surface_objects,display_tetrahedralized_volume_objects,
            display_hexahedralized_volume_objects,display_boundary_surface_objects,display_hard_bound_boundary_surface_objects,display_free_particles_objects);

    for(int i=1;i<=boundary_surface_objects.m;i++){
        if(!active_list(i)) continue;
        glPushName(i);
        //if(embedded_surface_objects(i) && display_mode==3){glPushName(3);embedded_surface_objects(i)->Display(in_color);glPopName();}
        if(boundary_surface_objects(i) && display_boundary_surface_objects){
          boundary_surface_objects(i)->wireframe_only=(display_mode==3);glPushName(5);boundary_surface_objects(i)->Display(in_color);glPopName();}
        if(hard_bound_boundary_surface_objects(i) && display_hard_bound_boundary_surface_objects){
            hard_bound_boundary_surface_objects(i)->wireframe_only=(display_mode==3);glPushName(5);hard_bound_boundary_surface_objects(i)->Display(in_color);glPopName();}
        glPopName();}
    if(slice && slice->Is_Slice_Mode()) glPopAttrib();

    if(interaction_pair_display_mode){
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glEnable(GL_BLEND);
        glDisable(GL_LIGHTING);

        // visualize point face interactions
        if(interaction_pair_display_mode==1 || interaction_pair_display_mode==2){
            OpenGL_Begin(GL_LINES);
            for(int k=1;k<=point_triangle_interaction_pairs.Size();k++){
                const POINT_FACE_REPULSION_PAIR<TV>& pair=point_triangle_interaction_pairs(k);
                INDIRECT_ARRAY<const ARRAY_VIEW<TV>,VECTOR<int,4>&> X(deformable_body_collection.particles.X,pair.nodes);
                glColor3f(.5f,1,.5f);
                OpenGL_Line(X(1),TRIANGLE_3D<T>(X(2),X(3),X(4)).Surface(X(1)));
                glColor3f(0,1,1);
                OpenGL_Line(X(1),X(2));OpenGL_Line(X(1),X(3));OpenGL_Line(X(1),X(4));}
            OpenGL_End();
            OpenGL_Begin(GL_TRIANGLES);
            for(int k=1;k<=point_triangle_interaction_pairs.Size();k++){
                const POINT_FACE_REPULSION_PAIR<TV>& pair=point_triangle_interaction_pairs(k);
                INDIRECT_ARRAY<const ARRAY_VIEW<TV>,VECTOR<int,4>&> X(deformable_body_collection.particles.X,pair.nodes);
                glColor4f(0,.6f,.8f,.5f);
                OpenGL_Triangle(X(1),X(3),X(2));
                OpenGL_Triangle(X(1),X(2),X(4));
                OpenGL_Triangle(X(1),X(4),X(3));}
            OpenGL_End();}

        // visualize edge edge interactions
        if(interaction_pair_display_mode==1 || interaction_pair_display_mode==3){
            OpenGL_Begin(GL_LINES);
            for(int k=1;k<=edge_edge_interaction_pairs.Size();k++){
                const EDGE_EDGE_REPULSION_PAIR<TV>& pair=edge_edge_interaction_pairs(k);
                INDIRECT_ARRAY<const ARRAY_VIEW<TV>,VECTOR<int,4>&> X(deformable_body_collection.particles.X,pair.nodes);
                glColor3f(1,1,0);
                VECTOR<T,2> weights;SEGMENT_3D<T>(X(1),X(2)).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(X(3),X(4)),weights);
                OpenGL_Line((1-weights.x)*X(1)+weights.x*X(2),(1-weights.y)*X(3)+weights.y*X(4));
                //OpenGL_Line(X(1),X(3));OpenGL_Line(X(1),X(4));OpenGL_Line(X(2),X(3));OpenGL_Line(X(2),X(4));
                glColor3f(1,.5,0);
                OpenGL_Line(X(1),X(2));OpenGL_Line(X(3),X(4));}
            OpenGL_End();}
        glPopAttrib();}

    //Draw force data
    if(display_forces_mode){
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glEnable(GL_BLEND);
        glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_LINES);
        for(int k=1;k<=force_data_list.Size();k++){
            const FORCE_DATA<TV>& force_data=force_data_list(k);
            if(display_forces_mode==1 && force_data.name!="LINEAR_SPRINGS") continue;
            else if(display_forces_mode==2 && force_data.name!="TRIANGLE_BENDING_SPRINGS") continue;
            else if(display_forces_mode==3 && force_data.name!="LINEAR_ALTITUDE_SPRINGS_3D") continue;
            else if(display_forces_mode==4 && force_data.name!="SCALED_SOLIDS_FORCES") continue;

            OPENGL_COLOR force_color=color_map_forces->Lookup(force_data.state);force_color.Send_To_GL_Pipeline();
            OpenGL_Line(force_data.first_action_point,force_data.second_action_point);}
        OpenGL_End();
        glPopAttrib();}
}
//#####################################################################
// Function Cycle_Hard_Bound_Surface_Display_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Cycle_Hard_Bound_Surface_Display_Mode()
{
    display_hard_bound_surface_mode=(display_hard_bound_surface_mode+1)%3;
    display_soft_bound_surface_mode=!has_embedded_objects && has_soft_bindings?display_hard_bound_surface_mode:1; // TODO: fix names
}
//#####################################################################
// Function Cycle_Forces_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Cycle_Forces_Mode()
{
    display_forces_mode=(display_forces_mode+1)%6;
    std::stringstream ss;
    if(display_forces_mode==0) ss<<"Displaying no forces"<<std::endl;
    else if(display_forces_mode==1) ss<<"Displaying forces for LINEAR_SPRINGS"<<std::endl;
    else if(display_forces_mode==2) ss<<"Displaying forces for TRIANGLE_BENDING_SPRINGS"<<std::endl;
    else if(display_forces_mode==3) ss<<"Displaying forces for LINEAR_ALTITUDE_SPRINGS_3D"<<std::endl;
    else if(display_forces_mode==4) ss<<"Displaying forces for SCALED_SOLIDS_FORCES"<<std::endl;
    else ss<<"Displaying all forces"<<std::endl;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box=BASE::Bounding_Box();
    if(draw && valid && deformable_body_collection.deformable_geometry.structures.m>0){
        for(int i=1;i<=boundary_surface_objects.m;i++) if(boundary_surface_objects(i))box.Enlarge_To_Include_Box(boundary_surface_objects(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Cycle_Interaction_Pair_Display_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Cycle_Interaction_Pair_Display_Mode()
{
    if(!interaction_pair_display_mode && !point_triangle_interaction_pairs.m && !edge_edge_interaction_pairs.m){
        std::string file=STRING_UTILITIES::string_sprintf("%s/%d/interaction_pairs",prefix.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(file))
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(RW()),file,point_triangle_interaction_pairs,edge_edge_interaction_pairs);}
    interaction_pair_display_mode=(interaction_pair_display_mode+1)%4;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* selection=0;
    if(buffer_size>=1){
        if(buffer[1]<=3) return BASE::Get_Selection(buffer,buffer_size);

        selection=new OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>(this);
        selection->body_index=buffer[0];selection->subobject_type=buffer[1];
        switch(selection->subobject_type){
            case 4:selection->subobject=embedded_surface_objects(buffer[0]);break;
            case 5:selection->subobject=boundary_surface_objects(buffer[0]);break;
            default:
                delete selection;
                selection=0;
                return BASE::Get_Selection(buffer,buffer_size);}
        selection->body_selection=selection->subobject->Get_Selection(&buffer[2],buffer_size-2);
        if(!selection->body_selection){delete selection;selection=0;}}
    return selection;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>::
Clear_Highlight()
{
    BASE::Clear_Highlight();
   for(int i=1;i<=boundary_surface_objects.m;i++) if(boundary_surface_objects(i) && active_list(i))boundary_surface_objects(i)->Clear_Highlight();
    for(int i=1;i<=embedded_surface_objects.m;i++) if(embedded_surface_objects(i) && active_list(i)) embedded_surface_objects(i)->Clear_Highlight();
    for(int i=1;i<=hard_bound_boundary_surface_objects.m;i++) if(hard_bound_boundary_surface_objects(i) && active_list(i))hard_bound_boundary_surface_objects(i)->Clear_Highlight();
}
//#####################################################################
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<double,double>;
#endif
