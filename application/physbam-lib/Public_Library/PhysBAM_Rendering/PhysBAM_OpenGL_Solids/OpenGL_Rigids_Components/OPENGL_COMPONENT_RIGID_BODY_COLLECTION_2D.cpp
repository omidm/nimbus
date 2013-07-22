//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
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
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(const std::string& basedir_input)
    :OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>(new RIGID_GEOMETRY_PARTICLES<TV>(),basedir_input),draw_articulation_points(false),
    draw_forces_and_torques(false),draw_linear_muscles(false),rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0,0)),articulated_rigid_body(0),
    need_destroy_rigid_body_collection(true)
{
    rigid_geometry_collection=&rigid_body_collection.rigid_geometry_collection;
    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir_input)
    :OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>(rigid_body_collection.rigid_geometry_collection,basedir_input),draw_articulation_points(false),
    draw_forces_and_torques(false),draw_linear_muscles(false),rigid_body_collection(rigid_body_collection),articulated_rigid_body(0),
    need_destroy_rigid_body_collection(false)
{
    rigid_geometry_collection=&rigid_body_collection.rigid_geometry_collection;
    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D()
{
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Read_Articulated_Information
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
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
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame))) return;

        rigid_body_collection.Read(STREAM_TYPE(RW()),basedir,frame,&needs_init,&needs_destroy); // TODO: avoiding reading triangulated areas

        std::string arb_state_file=STRING_UTILITIES::string_sprintf("%s/%d/arb_state",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(arb_state_file)){
            if(!articulated_rigid_body) articulated_rigid_body=new ARTICULATED_RIGID_BODY<TV>(rigid_body_collection); // TODO: read in the actual particles
            articulated_rigid_body->Read(STREAM_TYPE(RW()),basedir,frame);}
        else{delete articulated_rigid_body;articulated_rigid_body=0;}

        if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame)))
            Read_Articulated_Information(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame));

        std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_forces_and_torques",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename,forces_and_torques);
        else forces_and_torques.Resize(0);

        BASE::Reinitialize(force,false);}
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Update_Object_Labels()
{
    BASE::Update_Object_Labels();
    // TODO (zhw): fix this after SIGGRAPH
/*    for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        if(draw_object(i)){
            if(opengl_segmented_curve(i)){
                if(output_positions){
                rigid_body_collection.Rigid_Body(i).Update_Angular_Velocity();}}}}*/
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
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
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size>=2){
        if(buffer[0]<4)
            return BASE::Get_Selection(buffer,buffer_size);
        else if(buffer[0]==4){ // articulation joints
            selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>(this,buffer[1]);}
        else if(buffer[0]==5){ // muscles
            selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>(this,buffer[1]);}}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    BASE::Highlight_Selection(selection);
    if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_2D){
        current_selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>(this,((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_2D<T>*)selection)->joint_id);}
    else if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_MUSCLE_2D){
        current_selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>(this,((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_MUSCLE_2D<T>*)selection)->muscle_id);}
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    if(!selection) return;
    BASE::Print_Selection_Info(output_stream,selection);
    if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_2D){
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
// Function Toggle_Articulation_Points
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Toggle_Articulation_Points()
{
    draw_articulation_points=!draw_articulation_points;
}
//#####################################################################
// Function Toggle_Articulation_Points
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Toggle_Linear_Muscles()
{
    draw_linear_muscles=!draw_linear_muscles;
}
//#####################################################################
// Function Toggle_Forces_And_Torques
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T,RW>::
Toggle_Forces_And_Torques()
{
    draw_forces_and_torques=!draw_forces_and_torques;
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
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<double,double>;
#endif
