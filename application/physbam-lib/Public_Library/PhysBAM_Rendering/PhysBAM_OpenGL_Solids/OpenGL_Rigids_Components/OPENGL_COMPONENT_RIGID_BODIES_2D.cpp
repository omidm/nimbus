//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h> // TODO: remove once MUSCLE.cpp exists (workaround for windows compiler bug)
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODIES_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
OPENGL_COMPONENT_RIGID_BODIES_2D(const std::string& basedir_input)
    :OPENGL_COMPONENT("Rigid Bodies 2D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_individual_axes(false),draw_node_velocity_vectors(false),draw_articulation_points(false),draw_segmented_curve(true),draw_triangulated_area(false),draw_implicit_curve(false),
    draw_forces_and_torques(false),draw_linear_muscles(false),rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0,0)),articulated_rigid_body(0),
    velocity_field(velocity_vectors,positions,OPENGL_COLOR::Cyan(),0.25,true,true),node_velocity_field(node_velocity_vectors,node_positions,OPENGL_COLOR::Magenta(),0.25,true,true),
    current_selection(0),need_destroy_rigid_body_collection(true)
{
    is_animation=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
OPENGL_COMPONENT_RIGID_BODIES_2D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir_input)
    :OPENGL_COMPONENT("Rigid Bodies 2D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_individual_axes(false),draw_node_velocity_vectors(false),draw_articulation_points(false),draw_segmented_curve(true),draw_triangulated_area(false),draw_implicit_curve(false),
    draw_forces_and_torques(false),draw_linear_muscles(false),rigid_body_collection(rigid_body_collection),articulated_rigid_body(0),
    velocity_field(velocity_vectors,positions,OPENGL_COLOR::Cyan(),0.25,true,true),node_velocity_field(node_velocity_vectors,node_positions,OPENGL_COLOR::Magenta(),0.25,true,true),
    current_selection(0),need_destroy_rigid_body_collection(false)
{
    is_animation=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
~OPENGL_COMPONENT_RIGID_BODIES_2D()
{
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Set_Draw_Mode(const int mode)
{
    draw_segmented_curve=(mode&0x01)!=0;
    draw_triangulated_area=(mode&0x02)!=0;
    draw_implicit_curve=(mode&0x04)!=0;
}
//#####################################################################
// Function Get_Draw_Mode
//#####################################################################
template<class T,class RW> int OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Get_Draw_Mode() const
{
    return 0x04*(int)(draw_implicit_curve)+0x02*(int)(draw_triangulated_area)+0x01*(int)(draw_segmented_curve);
}
//#####################################################################
// Function Read_Hints
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Read_Hints(const std::string& filename)
{
    ARRAY<OPENGL_RIGID_BODY_HINTS,int> opengl_hints;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    Read_Binary<RW>(*input,opengl_hints);delete input;

    for(int i(1);i<=opengl_segmented_curve.Size();i++) if(opengl_segmented_curve(i) && i<=opengl_hints.Size()){
        opengl_segmented_curve(i)->color=opengl_hints(i).material.diffuse;
        use_object_bounding_box(i)=opengl_hints(i).include_bounding_box;}
}
//#####################################################################
// Function Read_Articulated_Information
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Read_Articulated_Information(const std::string& filename)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    // this will need to be changed to reflect multiple articulation points per rigid body
    int numpoints=0;Read_Binary<RW,int>(*input,numpoints);articulation_points.Exact_Resize(numpoints);
    for(int i=1;i<=numpoints;i++) Read_Binary<RW,VECTOR<T,2> >(*input,articulation_points(i));
    delete input;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Reinitialize(const bool force)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_bodies",basedir.c_str(),frame))) return;

        ARRAY<int> needs_init;
        rigid_body_collection.Read(STREAM_TYPE(RW()),basedir,frame,&needs_init); // TODO: avoiding reading triangulated areas

        std::string arb_state_file=STRING_UTILITIES::string_sprintf("%s/%d/arb_state",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(arb_state_file)){
            if(!articulated_rigid_body) articulated_rigid_body=new ARTICULATED_RIGID_BODY<TV>(rigid_body_collection); // TODO: read in the actual particles
            articulated_rigid_body->Read(STREAM_TYPE(RW()),basedir,frame);}
        else{delete articulated_rigid_body;articulated_rigid_body=0;}

        if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame)))
            Read_Articulated_Information(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame));
        int max_number_of_bodies(max(opengl_segmented_curve.Size(),rigid_body_collection.rigid_body_particle.array_collection->Size()));
        // only enlarge array as we read in more geometry to memory
        opengl_segmented_curve.Resize(max_number_of_bodies);
        opengl_triangulated_area.Resize(max_number_of_bodies);
        opengl_levelset.Resize(max_number_of_bodies);
        extra_components.Resize(max_number_of_bodies);
        opengl_axes.Resize(max_number_of_bodies);
        draw_object.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(draw_object,true);
        use_object_bounding_box.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(use_object_bounding_box,true);

        // Initialize bodies which have become active
        for(int i=1;i<=needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_body_collection.Is_Active(id));
            Create_Geometry(id);}

        // Update active bodies / remove inactive bodies
        for(int id(1);id<=rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
            if(rigid_body_collection.Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_body_collection.rigid_body_particle.array_collection->Size()+1;id<=opengl_segmented_curve.Size();id++) Destroy_Geometry(id);

        frame_loaded=frame;
        valid=true;

        std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_forces_and_torques",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename,forces_and_torques);
        else forces_and_torques.Resize(0);}

    Update_Object_Labels();
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Create_Geometry(const int id)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);

    if(!opengl_axes(id)){opengl_axes(id)=new OPENGL_AXES<T>();}
    if(rigid_body.simplicial_object && !opengl_segmented_curve(id)){
        opengl_segmented_curve(id)=new OPENGL_SEGMENTED_CURVE_2D<T>(*rigid_body.simplicial_object);
        opengl_segmented_curve(id)->draw_vertices=true;
        opengl_segmented_curve(id)->Enslave_Transform_To(*opengl_axes(id));}

    // add triangulated area
    TRIANGULATED_AREA<T>* triangulated_area=rigid_body.template Find_Structure<TRIANGULATED_AREA<T>*>();
    if(triangulated_area && !opengl_triangulated_area(id)){
        triangulated_area->mesh.Initialize_Segment_Mesh();
        opengl_triangulated_area(id)=new OPENGL_TRIANGULATED_AREA<T>(*triangulated_area);
        opengl_triangulated_area(id)->Enslave_Transform_To(*opengl_axes(id));}

    if(rigid_body.implicit_object && !opengl_levelset(id)){ // ASSUMES LEVELSET_IMPLICIT_CURVE!
        IMPLICIT_OBJECT<TV>* object_space_implicit_object=rigid_body.implicit_object->object_space_implicit_object;
        if(LEVELSET_IMPLICIT_OBJECT<TV>* levelset_implicit_object=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(object_space_implicit_object)){
            opengl_levelset(id)=new OPENGL_LEVELSET_2D<T>(levelset_implicit_object->levelset);
            // TODO: fix me
            //opengl_levelset(id)->Set_Color_Mode(OPENGL_LEVELSET_2D<T>::COLOR_GRADIENT);
            opengl_levelset(id)->draw_cells=true;opengl_levelset(id)->draw_curve=false;opengl_levelset(id)->draw_area=false;
            opengl_levelset(id)->Enslave_Transform_To(*opengl_axes(id));}}

    // add extra components
    if(opengl_triangulated_area(id)){
        std::string filename_pattern=STRING_UTILITIES::string_sprintf("%s/accumulated_impulses_%d.%%d",basedir.c_str(),id);
        if(FILE_UTILITIES::Frame_File_Exists(filename_pattern,frame)){
            std::stringstream ss;
            ss<<"Adding accumulated impulses to rigid body "<<id<<std::endl;
            LOG::filecout(ss.str());
            OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>* component=
                new OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>(*triangulated_area,filename_pattern);
            component->opengl_vector_field.Enslave_Transform_To(*opengl_axes(id));
            extra_components(id).Append(component);}}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Update_Geometry(const int id)
{
    if(opengl_axes(id)) *opengl_axes(id)->frame=FRAME<VECTOR<float,3> >(Convert_2d_To_3d(rigid_body_collection.Rigid_Body(id).Frame()));
    if(opengl_triangulated_area(id)){
        std::string color_map_filename=STRING_UTILITIES::string_sprintf("%s/%d/stress_map_of_triangulated_area_%d",basedir.c_str(),frame,id);
        if(FILE_UTILITIES::File_Exists(color_map_filename)){
            if(!opengl_triangulated_area(id)->color_map) opengl_triangulated_area(id)->color_map=new ARRAY<OPENGL_COLOR>;
            FILE_UTILITIES::Read_From_File<RW>(color_map_filename,*opengl_triangulated_area(id)->color_map);}
        else if(opengl_triangulated_area(id)->color_map){delete opengl_triangulated_area(id)->color_map;opengl_triangulated_area(id)->color_map=0;}}
    for(int i=1;i<=extra_components(id).m;i++) extra_components(id)(i)->Set_Frame(frame);
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Destroy_Geometry(const int id)
{
    if(opengl_segmented_curve(id)){delete opengl_segmented_curve(id);opengl_segmented_curve(id)=0;}
    if(opengl_triangulated_area(id)){delete opengl_triangulated_area(id)->color_map;delete opengl_triangulated_area(id);opengl_triangulated_area(id)=0;}
    if(opengl_levelset(id)){delete opengl_levelset(id);opengl_levelset(id)=0;}
    if(opengl_axes(id)){delete opengl_axes(id);opengl_axes(id)=0;}
    extra_components(id).Delete_Pointers_And_Clean_Memory();
    draw_object(id)=false;
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Update_Object_Labels()
{
    int number_of_drawn_bodies=draw_object.Count_Matches(true);
    if(draw_velocity_vectors){
        positions.Resize(number_of_drawn_bodies);
        velocity_vectors.Resize(number_of_drawn_bodies);}
    if(draw_node_velocity_vectors){
        node_positions.Resize(number_of_drawn_bodies*4);
        node_velocity_vectors.Resize(number_of_drawn_bodies*4);}

    int idx=0;
    for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        if(draw_object(i)){
            if(draw_velocity_vectors || draw_node_velocity_vectors){idx++;
                if(draw_velocity_vectors){
                    positions(idx)=rigid_body_collection.rigid_body_particle.X(i);
                    velocity_vectors(idx)=rigid_body_collection.rigid_body_particle.V(i);}
                if(draw_node_velocity_vectors){
                    //only valid for squares . . .
                    RIGID_BODY<TV>* rigid_body=&rigid_body_collection.Rigid_Body(i);

                    node_positions((idx-1)*4+1)=rigid_body->World_Space_Point(VECTOR<T,2>(1,1));
                    node_velocity_vectors((idx-1)*4+1)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(1,1)));

                    node_positions((idx-1)*4+2)=rigid_body->World_Space_Point(VECTOR<T,2>(1,-1));
                    node_velocity_vectors((idx-1)*4+2)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(1,-1)));

                    node_positions((idx-1)*4+3)=rigid_body->World_Space_Point(VECTOR<T,2>(-1,1));
                    node_velocity_vectors((idx-1)*4+3)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(-1,1)));

                    node_positions((idx-1)*4+4)=rigid_body->World_Space_Point(VECTOR<T,2>(-1,-1));
                    node_velocity_vectors((idx-1)*4+4)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(-1,-1)));}}
            if(opengl_segmented_curve(i)){
                if(output_positions){
                    rigid_body_collection.Rigid_Body(i).Update_Angular_Velocity();
                    opengl_segmented_curve(i)->Set_Name(STRING_UTILITIES::string_sprintf("%s <%.3f %.3f> [w=%.3f]",rigid_body_collection.Rigid_Body(i).name.c_str(),rigid_body_collection.rigid_body_particle.X(i).x,rigid_body_collection.rigid_body_particle.X(i).y,rigid_body_collection.rigid_body_particle.angular_velocity(i).x));}
                else opengl_segmented_curve(i)->Set_Name(rigid_body_collection.Rigid_Body(i).name);}}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_bodies",basedir.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Display(const int in_color) const
{
    if(draw){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);

        GLint mode;
        glGetIntegerv(GL_RENDER_MODE,&mode);

        if(draw_segmented_curve){
            glPushName(1);
            for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
                glPushName(Value(i));
                if(draw_object(i) && opengl_segmented_curve(i)) opengl_segmented_curve(i)->Display(in_color);
                glPopName();}
            glPopName();}
        if(draw_triangulated_area){
            glPushName(2);
            for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
                glPushName(Value(i));
                if(draw_object(i) && opengl_triangulated_area(i)) opengl_triangulated_area(i)->Display(in_color);
                glPopName();}
            glPopName();}
        if(draw_implicit_curve){
            glPushName(3);
            for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
                glPushName(Value(i));
                if(draw_object(i) && opengl_levelset(i)) opengl_levelset(i)->Display(in_color);
                glPopName();}
            glPopName();}
        if(draw_individual_axes)
            for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++)
                if(draw_object(i) && opengl_axes(i)) opengl_axes(i)->Display(in_color);

        // Articulated rigid bodies
        if(articulated_rigid_body){
            if(draw_articulation_points){
                OPENGL_COLOR articulation_point_color(0.5,0.5,0.5),segment_color(0,0,1);
                if(mode==GL_SELECT){glPushName(4);glPushAttrib(GL_POINT_BIT);glPointSize(OPENGL_PREFERENCES::selection_point_size);}
                for(int i=1;i<=articulation_points.m;i++){
                    glPushName(i);
                    OPENGL_SHAPES::Draw_Dot(articulation_points(i),articulation_point_color,5);
                    glPopName();}
                if(mode!=GL_SELECT) for(int i=1;i<=articulation_points.m;i+=2){
                    OPENGL_SHAPES::Draw_Segment(articulation_points(i),articulation_points(i+1),segment_color,5);}
                if(mode==GL_SELECT){glPopName();glPopAttrib();}
                if(mode!=GL_SELECT && current_selection && current_selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_2D){
                    int joint_id=((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>*)current_selection)->joint_id;
                    OPENGL_SELECTION::Draw_Highlighted_Vertex(articulation_points(joint_id));}}

            if(draw_linear_muscles){
                glPushAttrib(GL_LINE_BIT|GL_ENABLE_BIT|GL_CURRENT_BIT);
                // draw all muscles
                OPENGL_COLOR muscle_color(1,0,0);muscle_color.Send_To_GL_Pipeline();
                glLineWidth(mode==GL_SELECT?OPENGL_PREFERENCES::selection_line_width:5);
                glPushName(5);glPushName(0);
                for(int i=1;i<=articulated_rigid_body->muscle_list->muscles.m;i++){
                    MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(i);
                    glLoadName(i);OpenGL_Begin(GL_LINE_STRIP);
                    OpenGL_Vertex(muscle.attachment_point_1->Embedded_Position());
                    for(int t=1;t<=muscle.via_points.m;t++) OpenGL_Vertex(muscle.via_points(t)->Embedded_Position());
                    OpenGL_Vertex(muscle.attachment_point_2->Embedded_Position());
                    OpenGL_End();}
                glPopName();glPopName();
                // highligh selected one
                if(mode!=GL_SELECT && current_selection && current_selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_MUSCLE_2D){
                    int muscle_id=((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>*)current_selection)->muscle_id;
                    MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(muscle_id);
                    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);OPENGL_PREFERENCES::selection_highlight_color.Send_To_GL_Pipeline();
                    OpenGL_Begin(GL_LINE_STRIP);
                    OpenGL_Vertex(muscle.attachment_point_1->Embedded_Position());
                    for(int t=1;t<=muscle.via_points.m;t++) OpenGL_Vertex(muscle.via_points(t)->Embedded_Position());
                    OpenGL_Vertex(muscle.attachment_point_2->Embedded_Position());
                    OpenGL_End();}
                glPopAttrib();}}

        if(mode!=GL_SELECT){
            if(draw_velocity_vectors) velocity_field.Display(in_color);
            if(draw_node_velocity_vectors) node_velocity_field.Display(in_color);

            if(draw_forces_and_torques && forces_and_torques.Size()==rigid_body_collection.rigid_body_particle.array_collection->Size()){
                T scale=(T)velocity_field.size/24;
                OPENGL_COLOR::Yellow().Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINES);
                for(int i(1);i<=forces_and_torques.Size();i++)
                    OPENGL_SHAPES::Draw_Arrow(rigid_body_collection.rigid_body_particle.X(i),rigid_body_collection.rigid_body_particle.X(i)+scale*forces_and_torques(i).x);
                OpenGL_End();
                for(int i(1);i<=forces_and_torques.Size();i++){
                    std::string label=STRING_UTILITIES::string_sprintf("F=%.3f %.3f, T=%.3f",forces_and_torques(i).x.x,forces_and_torques(i).x.y,forces_and_torques(i).y);
                    OpenGL_String(rigid_body_collection.rigid_body_particle.X(i)+scale*forces_and_torques(i).x,label);}}

            for(int i(1);i<=extra_components.Size();i++)
                for(int j=1;j<=extra_components(i).m;j++)
                    extra_components(i)(j)->Display(in_color);

            if(show_object_names){
                glColor3f(1,1,1);
                for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++)
                    if(draw_object(i) && rigid_body_collection.Rigid_Body(i).name.length())
                        OpenGL_String(rigid_body_collection.rigid_body_particle.X(i),rigid_body_collection.Rigid_Body(i).name);}}
        glPopAttrib();}
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Use_Bounding_Box() const
{
    int num_drawable_objects=0;
    for(int i(1);i<=opengl_segmented_curve.Size();i++)
        if(draw_object(i) && use_object_bounding_box(i))
            num_drawable_objects++;
    return (draw && num_drawable_objects>0);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box;
    if(draw){
        bool first=true;
        for(int i(1);i<=opengl_segmented_curve.Size();i++)
            if(draw_object(i) && use_object_bounding_box(i) && opengl_segmented_curve(i)){
                if(first){
                    box=opengl_segmented_curve(i)->Bounding_Box();
                    first=false;}
                else{box=RANGE<VECTOR<float,3> >::Combine(box,opengl_segmented_curve(i)->Bounding_Box());}}}
    return box;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size>=2){
        if(buffer[0]==1 || buffer[0]==2){
            int body_id(buffer[1]);
            OPENGL_SELECTION* body_selection=0;
            if(buffer[0]==1){ // segmented curve
                PHYSBAM_ASSERT(opengl_segmented_curve(body_id));
                body_selection=opengl_segmented_curve(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
            else{ // triangulated area
                PHYSBAM_ASSERT(opengl_triangulated_area(body_id));
                body_selection=opengl_triangulated_area(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
            if(body_selection) selection=new OPENGL_SELECTION_COMPONENT_RIGID_BODIES_2D<T>(this,body_id,body_selection);}
        else if(buffer[0]==4){ // articulation joints
            selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>(this,buffer[1]);}
        else if(buffer[0]==5){ // muscles
            selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>(this,buffer[1]);}}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if(selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_2D){
        OPENGL_SELECTION_COMPONENT_RIGID_BODIES_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_BODIES_2D<T>*)selection;
        if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D ||
           real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D){ // segmented curve
            if(opengl_segmented_curve(real_selection->body_id)) // might have become inactive
                opengl_segmented_curve(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}
        else if(real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE){ // triangulated area
            if(opengl_triangulated_area(real_selection->body_id)) // might have become inactive
                opengl_triangulated_area(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}}
    else if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_2D){
        current_selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>(this,((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>*)selection)->joint_id);}
    else if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_MUSCLE_2D){
        current_selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>(this,((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>*)selection)->muscle_id);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
    for(int i(1);i<=opengl_segmented_curve.Size();i++)
        if(opengl_segmented_curve(i)) opengl_segmented_curve(i)->Clear_Highlight();
    for(int i(1);i<=opengl_triangulated_area.Size();i++)
        if(opengl_triangulated_area(i)) opengl_triangulated_area(i)->Clear_Highlight();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    if(!selection) return;
    if(selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_2D){
        OPENGL_SELECTION_COMPONENT_RIGID_BODIES_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_BODIES_2D<T>*)selection;

        output_stream<<"Rigid body "<<real_selection->body_id<<std::endl;

        // handle case of body no longer being active
        if(!rigid_body_collection.Is_Active(real_selection->body_id) || !opengl_segmented_curve(real_selection->body_id)){
            output_stream<<"INACTIVE"<<std::endl;
            return;}

        RIGID_BODY<TV> *body=const_cast<RIGID_BODY<TV>*>(&rigid_body_collection.Rigid_Body(real_selection->body_id));
        body->Update_Angular_Velocity();

        if(!body->name.empty()){output_stream<<"Name = "<<body->name<<std::endl;}
        output_stream<<"Mass = "<<body->Mass()<<std::endl;
        output_stream<<"Moment of inertia = "<<body->Inertia_Tensor()<<std::endl;
        output_stream<<"Position = "<<body->X()<<std::endl;
        output_stream<<"Orientation = "<<body->Rotation()<<" (angle "<<body->Rotation().Angle()<<")"<<std::endl;
        output_stream<<"Velocity = "<<body->V()<<std::endl;
        output_stream<<"Angular velocity = "<<body->Angular_Velocity()<<std::endl;
        output_stream<<std::endl;

        MATRIX<T,3> body_transform=body->Frame().Matrix();

        // Have to do this ourselves here because we need to transform particle positions to world space
        if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D ||
           real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D){ // segmented curve
            if(opengl_segmented_curve(real_selection->body_id))
                opengl_segmented_curve(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);

            if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D){
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)real_selection->body_selection)->index;
                output_stream<<"Pointwise object velocity "<<body->Pointwise_Object_Velocity(body->World_Space_Point(body->simplicial_object->particles.X(index)))<<std::endl;}
            else if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D){
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)real_selection->body_selection)->index;
                int node1,node2;body->simplicial_object->mesh.elements(index).Get(node1,node2);
                output_stream<<"Pointwise object velocity "<<body->Pointwise_Object_Velocity(body->World_Space_Point(body->simplicial_object->particles.X(node1)))<<std::endl;
                output_stream<<"Pointwise object velocity "<<body->Pointwise_Object_Velocity(body->World_Space_Point(body->simplicial_object->particles.X(node2)))<<std::endl;}}
        else if(real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE){
            if(opengl_triangulated_area(real_selection->body_id))
                opengl_triangulated_area(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);}}
    else if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_2D){
        PHYSBAM_ASSERT(articulated_rigid_body);
        OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T> *real_selection=(OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>*)selection;
        JOINT_ID joint_id((real_selection->joint_id+1)/2);
        JOINT<TV>* joint=articulated_rigid_body->joint_mesh(joint_id);
        const RIGID_BODY<TV>* parent=articulated_rigid_body->Parent(joint_id),*child=articulated_rigid_body->Child(joint_id);
        FRAME<TV> current_frame=joint->Compute_Current_Joint_Frame(*parent,*child);
        output_stream<<"Joint id "<<joint_id<<" ("<<(!joint->name.empty()?joint->name:"UNNAMED")<<")"<<std::endl;
        output_stream<<"Frame = "<<(real_selection->joint_id%2==0?"child":"parent")<<std::endl;
        //output_stream<<"Articulation point = "<<articulation_points(real_selection->joint_id)<<std::endl;
        VECTOR<T,2> ap1=parent->World_Space_Point(joint->F_pj().t),ap2=child->World_Space_Point(joint->F_cj().t),location=(T).5*(ap1+ap2);
        output_stream<<"Parent = "<<parent->name<<std::endl;
        output_stream<<"Child = "<<child->name<<std::endl;
        output_stream<<"Joint translation = "<<current_frame.t<<std::endl;
        output_stream<<"Joint rotation angle = "<<current_frame.r.Angle()<<std::endl;

        VECTOR<T,2> current_relative_velocity=-RIGID_BODY<TV>::Relative_Velocity(*parent,*child,location); // child w.r.t. parent!
        VECTOR<T,1> current_relative_angular_velocity=-RIGID_BODY<TV>::Relative_Angular_Velocity(*parent,*child); // child w.r.t. parent!
        output_stream<<"Relative velocity at joint = "<<current_relative_velocity<<std::endl;
        output_stream<<"Relative angular velocity = "<<current_relative_angular_velocity<<std::endl;}
    else if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_MUSCLE_2D){
        PHYSBAM_ASSERT(articulated_rigid_body);
        OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T> *real_selection=(OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>*)selection;
        int muscle_id=real_selection->muscle_id;MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(muscle_id);
        ATTACHMENT_POINT<TV>* attachment_point_1=muscle.attachment_point_1;
        ATTACHMENT_POINT<TV>* attachment_point_2=muscle.attachment_point_2;
        output_stream<<"Muscle "<<muscle_id<<" ("<<(!muscle.name.empty()?muscle.name:"UNNAMED")<<")"<<std::endl;
        output_stream<<"Optimal length = "<<muscle.optimal_length<<std::endl;
        output_stream<<"Peak force = "<<muscle.peak_force<<std::endl;
        output_stream<<"Pennation angle = "<<muscle.pennation_angle<<std::endl;
        output_stream<<"Tendon slack length = "<<muscle.tendon_slack_length<<std::endl;
        output_stream<<"Maximum shortening velocity = "<<muscle.max_shortening_velocity<<std::endl;
        output_stream<<"Origin = ("<<attachment_point_1->rigid_body.name<<", "<<attachment_point_1->object_space_position<<")"<<std::endl;
        output_stream<<"Insertion = ("<<attachment_point_2->rigid_body.name<<", "<<attachment_point_2->object_space_position<<")"<<std::endl;
        if(muscle.via_points.m){
            output_stream<<"Via points: ";
            for(int i=1;i<=muscle.via_points.m;i++)
                output_stream<<"("<<muscle.via_points(i)->rigid_body.name<<", "<<muscle.via_points(i)->object_space_position<<") ";
            output_stream<<std::endl;}
        output_stream<<std::endl;
        output_stream<<"Total length = "<<muscle.Total_Length()<<std::endl;
        output_stream<<"Total velocity = "<<muscle.Total_Velocity()<<std::endl;
        if(articulated_rigid_body->muscle_activations.Valid_Index(muscle_id)) output_stream<<"Activation = "<<articulated_rigid_body->muscle_activations(muscle_id)<<std::endl;}
}
//#####################################################################
// Function Set_Draw_Object
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Set_Draw_Object(int i,bool draw_it)
{
    if(i>draw_object.Size()) draw_object.Resize(i);
    draw_object(i)=draw_it;
}
//#####################################################################
// Function Set_Object_Color
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Set_Object_Color(int i,const OPENGL_COLOR &color_input)
{
    if(!opengl_segmented_curve(i)) return;
    opengl_segmented_curve(i)->color=color_input;
    opengl_segmented_curve(i)->color_gray=color_input.Grayscale();
}
//#####################################################################
// Function Get_Draw_Object
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Get_Draw_Object(int i) const
{
    return draw_object(i);
}
//#####################################################################
// Function Set_Use_Object_Bounding_Box
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Set_Use_Object_Bounding_Box(int i,bool use_it)
{
    use_object_bounding_box(i)=use_it;
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Set_Vector_Size(double size)
{
    velocity_field.size=size;
    node_velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Individual_Axes
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Individual_Axes()
{
    draw_individual_axes=!draw_individual_axes;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Output_Positions
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Output_Positions()
{
    output_positions=!output_positions;
    Update_Object_Labels();
}
//#####################################################################
// Function Toggle_Show_Object_Names
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Show_Object_Names()
{
    show_object_names=!show_object_names;
}
//#####################################################################
// Function Toggle_Node_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Node_Velocity_Vectors()
{
    draw_node_velocity_vectors=!draw_node_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Articulation_Points
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Articulation_Points()
{
    draw_articulation_points=!draw_articulation_points;
}
//#####################################################################
// Function Toggle_Articulation_Points
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Linear_Muscles()
{
    draw_linear_muscles=!draw_linear_muscles;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Increase_Vector_Size()
{
    velocity_field.Scale_Vector_Size((T)1.1);
    node_velocity_field.Scale_Vector_Size((T)1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Decrease_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1/(T)1.1);
    node_velocity_field.Scale_Vector_Size(1/(T)1.1);
}
//#####################################################################
// Function Toggle_Forces_And_Torques
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Forces_And_Torques()
{
    draw_forces_and_torques=!draw_forces_and_torques;
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODIES_2D<T,RW>::
Toggle_Draw_Mode()
{
    int mode=Get_Draw_Mode();
    mode=(mode+1)%8;
    Set_Draw_Mode(mode);
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_RIGID_BODIES_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object && body_selection);
    return object->World_Space_Box(body_selection->Bounding_Box());
}
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Empty_Box();
}
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Empty_Box();
}
//#####################################################################
template class OPENGL_COMPONENT_RIGID_BODIES_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RIGID_BODIES_2D<double,double>;
#endif
