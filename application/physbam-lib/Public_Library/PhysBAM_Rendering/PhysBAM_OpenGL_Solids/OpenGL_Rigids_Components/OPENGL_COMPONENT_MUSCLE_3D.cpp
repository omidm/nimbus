//#####################################################################
// Copyright 2006-2007, Kevin Der, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ANALYTIC_SURFACE_MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_MUSCLE_3D.h>
using namespace PhysBAM;
template<class T> OPENGL_COMPONENT_MUSCLE_3D<T>::
OPENGL_COMPONENT_MUSCLE_3D(OPENGL_COLOR_MAP<T>* muscle_color_map_input)
    :OPENGL_COMPONENT("Muscles"),draw_linear_muscles(true),draw_surface_muscles(true),draw_muscle_internal_particles(false),muscle_color_map(muscle_color_map_input),
    current_selection(0),articulated_rigid_body(0)
{
}
template<class T> OPENGL_COMPONENT_MUSCLE_3D<T>::
~OPENGL_COMPONENT_MUSCLE_3D()
{
    opengl_triangulated_surface.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_MUSCLE_3D<T>::
Display(const int in_color) const
{
    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    if(draw_linear_muscles){
        glPushAttrib(GL_LINE_BIT|GL_ENABLE_BIT|GL_CURRENT_BIT);
        // draw all muscles
        OPENGL_COLOR muscle_color(1,0,0);muscle_color.Send_To_GL_Pipeline();
        glLineWidth(mode==GL_SELECT?OPENGL_PREFERENCES::selection_line_width:2);
        glPushName(5);glDisable(GL_LIGHTING);
        for(int i=1;i<=articulated_rigid_body->muscle_list->muscles.m;i++){
            T activation=articulated_rigid_body->muscle_activations.Valid_Index(i)?articulated_rigid_body->muscle_activations(i):0;
            muscle_color_map->Lookup((activation<0)?-log(1-activation):log(1+activation)).Send_To_GL_Pipeline();
            MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(i);
            if(mode==GL_SELECT){glPushName(i);glPushName(0);
                if(!muscle.via_points.m){glLoadName(1);OpenGL_Begin(GL_LINES);OpenGL_Vertex(muscle.attachment_point_1->Embedded_Position());
                    OpenGL_Vertex(muscle.attachment_point_2->Embedded_Position());OpenGL_End();}
                else{glLoadName(1);OpenGL_Begin(GL_LINES);OpenGL_Line(muscle.attachment_point_1->Embedded_Position(),muscle.via_points(1)->Embedded_Position());OpenGL_End();
                for(int t=1;t<muscle.via_points.m;t++){glLoadName(t+1);OpenGL_Begin(GL_LINES);OpenGL_Line(muscle.via_points(t)->Embedded_Position(),muscle.via_points(t+1)->Embedded_Position());OpenGL_End();}
                glLoadName(muscle.via_points.m+1);OpenGL_Begin(GL_LINES);OpenGL_Line(muscle.via_points(muscle.via_points.m)->Embedded_Position(),muscle.attachment_point_2->Embedded_Position());OpenGL_End();}
                glPopName();glPopName();}
            else{OpenGL_Begin(GL_LINES);
                if(!muscle.via_points.m){OpenGL_Line(muscle.attachment_point_1->Embedded_Position(),muscle.attachment_point_2->Embedded_Position());}
                else{OpenGL_Line(muscle.attachment_point_1->Embedded_Position(),muscle.via_points(1)->Embedded_Position());
                for(int t=1;t<muscle.via_points.m;t++){OpenGL_Line(muscle.via_points(t)->Embedded_Position(),muscle.via_points(t+1)->Embedded_Position());}
                OpenGL_Line(muscle.via_points(muscle.via_points.m)->Embedded_Position(),muscle.attachment_point_2->Embedded_Position());}OpenGL_End();}}
        glPopName();
        // highlight selected one
        if(mode!=GL_SELECT && current_selection && current_selection->type==OPENGL_SELECTION::MUSCLE_3D){
            glPushAttrib(GL_DEPTH_BUFFER_BIT);glDepthFunc(GL_LEQUAL);
            OPENGL_SELECTION_MUSCLE_3D<T>* muscle_selection=(OPENGL_SELECTION_MUSCLE_3D<T>*)current_selection;
            int muscle_id=muscle_selection->muscle_id;int segment_id=muscle_selection->segment_id;
            MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(muscle_id);
            glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);OPENGL_PREFERENCES::selection_highlight_color.Send_To_GL_Pipeline();
            OpenGL_Begin(GL_LINES);
            if(!muscle.via_points.m){PHYSBAM_ASSERT(segment_id==1);OpenGL_Line(muscle.attachment_point_1->Embedded_Position(),muscle.attachment_point_2->Embedded_Position());}
            else{PHYSBAM_ASSERT(segment_id<=muscle.via_points.m+1);
                if(segment_id==1){OpenGL_Line(muscle.attachment_point_1->Embedded_Position(),muscle.via_points(1)->Embedded_Position());}
                else if(segment_id==muscle.via_points.m+1){OpenGL_Line(muscle.via_points(muscle.via_points.m)->Embedded_Position(),muscle.attachment_point_2->Embedded_Position());}
                else{OpenGL_Line(muscle.via_points(segment_id-1)->Embedded_Position(),muscle.via_points(segment_id)->Embedded_Position());}}
            OpenGL_End();glPopAttrib();}
        glPopAttrib();}
    if(draw_surface_muscles){
        glPushName(6);
        for(int i=1;i<=opengl_triangulated_surface.m;i++){
            int muscle_index=surface_muscle_indices(i).x;int segment_index=surface_muscle_indices(i).y;
            if(mode==GL_SELECT){glPushName(muscle_index);glPushName(segment_index);glPushName(i);}
            ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>* segment=(ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>*)articulated_rigid_body->muscle_list->muscles(muscle_index)->muscle_segments(segment_index);
            segment->Get_Local_Positions_For_Particles(length_resolution,radial_resolution,opengl_triangulated_surface(i)->surface.particles);
            glPushMatrix();OpenGL_Translate(segment->frame.t);OpenGL_Rotate(segment->frame.r);
            opengl_triangulated_surface(i)->Display();glPopMatrix();if(mode==GL_SELECT){glPopName();glPopName();glPopName();}}
        glPopName();
        if(mode!=GL_SELECT && current_selection && current_selection->type==OPENGL_SELECTION::MUSCLE_SURFACE_3D){
            glPushAttrib(GL_DEPTH_BUFFER_BIT);glDepthFunc(GL_LEQUAL);glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);glDisable(GL_LIGHTING);
            OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>* muscle_selection=(OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>*)current_selection;
            int surface_index=muscle_selection->surface_index;int muscle_index=muscle_selection->muscle_id;int segment_index=muscle_selection->segment_id;
            ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>* segment=(ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>*)articulated_rigid_body->muscle_list->muscles(muscle_index)->muscle_segments(segment_index);
            glPushMatrix();OpenGL_Translate(segment->frame.t);OpenGL_Rotate(segment->frame.r);glLineWidth(2);OPENGL_PREFERENCES::selection_highlight_color.Send_To_GL_Pipeline();
            opengl_triangulated_surface(surface_index)->Draw();glLineWidth(OPENGL_PREFERENCES::line_width);glPopMatrix();glPopAttrib();}}
    if(draw_muscle_internal_particles){
        for(int i=1;i<=muscle_internal_particles.m;i++){glPushName(i);OPENGL_SHAPES::Draw_Dot(muscle_internal_particles(i),OPENGL_COLOR::Red(),5);glPopName();}}
}
//#####################################################################
// Function Read_Muscle_Internal_Particles
//#####################################################################
template<class T> void OPENGL_COMPONENT_MUSCLE_3D<T>::
Read_Muscle_Internal_Particles(const std::string& filename)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    int numpoints=0;Read_Binary<T>(*input,numpoints);muscle_internal_particles.Exact_Resize(numpoints);
    for(int i=1;i<=numpoints;i++) Read_Binary<T>(*input,muscle_internal_particles(i));
    delete input;
}
//#####################################################################
// Function Get_Muscle_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_COMPONENT_MUSCLE_3D<T>::
Get_Muscle_Selection(GLuint *buffer,int buffer_size)
{
    PHYSBAM_ASSERT(buffer_size>=3);
    OPENGL_SELECTION* selection=0;
    if(buffer[0]==5){selection=new OPENGL_SELECTION_MUSCLE_3D<T>(this,buffer[1],buffer[2]);}
    else if(buffer[0]==6){selection=new OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>(this,buffer[1],buffer[2],buffer[3]);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_MUSCLE_3D<T>::
Highlight_Selection(OPENGL_SELECTION* selection)
{
    delete current_selection;current_selection=0;
    if(selection->type==OPENGL_SELECTION::MUSCLE_3D){
        OPENGL_SELECTION_MUSCLE_3D<T>* muscle=(OPENGL_SELECTION_MUSCLE_3D<T>*)selection;
        current_selection=new OPENGL_SELECTION_MUSCLE_3D<T>(this,muscle->muscle_id,muscle->segment_id);}
    else if(selection->type==OPENGL_SELECTION::MUSCLE_SURFACE_3D){
        OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>* muscle=(OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>*)selection;
        current_selection=new OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>(this,muscle->muscle_id,muscle->segment_id,muscle->surface_index);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_MUSCLE_3D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_MUSCLE_3D<T>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const
{
    PHYSBAM_ASSERT(articulated_rigid_body);
    if(selection && selection->object==this){
        int muscle_id=0,segment_id=0;
        if(selection->type==OPENGL_SELECTION::MUSCLE_3D){
            OPENGL_SELECTION_MUSCLE_3D<T> *real_selection=(OPENGL_SELECTION_MUSCLE_3D<T>*)selection;
            muscle_id=real_selection->muscle_id;segment_id=real_selection->segment_id;}
        else if(selection->type==OPENGL_SELECTION::MUSCLE_SURFACE_3D){
            OPENGL_SELECTION_MUSCLE_SURFACE_3D<T> *real_selection=(OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>*)selection;
            muscle_id=real_selection->muscle_id;segment_id=real_selection->segment_id;}
        MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(muscle_id);
        ATTACHMENT_POINT<TV>* attachment_point_1=muscle.attachment_point_1;
        ATTACHMENT_POINT<TV>* attachment_point_2=muscle.attachment_point_2;
        if(selection->type==OPENGL_SELECTION::MUSCLE_SURFACE_3D){
            ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>* surface_segment=(ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>*)muscle.muscle_segments(segment_id);
            output_stream << "Surface type " << surface_segment->curve_type << std::endl;
            output_stream << "Thickness " << surface_segment->curve_thickness << std::endl;
            output_stream << "Offset " << surface_segment->curve_offset_thickness << std::endl; 
            output_stream << "Volume " << surface_segment->Compute_Volume() << std::endl;
                output_stream << "Tendon fraction " << surface_segment->Tendon_Length()/surface_segment->Length() << std::endl << std::endl;}
            output_stream << "Segment " << segment_id << std::endl;
            output_stream << "Type " << muscle.muscle_segments(segment_id)->Name() << std::endl;
            output_stream << "Length " << muscle.muscle_segments(segment_id)->Length() << std::endl << std::endl;
            output_stream << "Muscle " << muscle_id << " (" << (!muscle.name.empty()?muscle.name:"UNNAMED") << ")" << std::endl;
            output_stream << "Optimal length = " << muscle.optimal_length << std::endl;
            output_stream << "Peak force = " << muscle.peak_force << std::endl;
            output_stream << "Pennation angle = " << muscle.pennation_angle << std::endl;
            output_stream << "Tendon slack length = " << muscle.tendon_slack_length << std::endl;
            output_stream << "Maximum shortening velocity = " << muscle.max_shortening_velocity << std::endl;
            output_stream << "Origin = (" << attachment_point_1->rigid_body.name << ", " << attachment_point_1->object_space_position << ")" << std::endl;
            output_stream << "Insertion = (" << attachment_point_2->rigid_body.name << ", " << attachment_point_2->object_space_position << ")" << std::endl;
            if(muscle.via_points.m){
                output_stream << "Via points: ";
                for(int i=1;i<=muscle.via_points.m;i++) 
                    output_stream << "(" << muscle.via_points(i)->rigid_body.name << ", " << muscle.via_points(i)->object_space_position << ") ";
                output_stream << std::endl;}
            output_stream << std::endl;
            output_stream << "Total length = " << muscle.Total_Length() << std::endl;
            output_stream << "Total velocity = " << muscle.Total_Velocity() << std::endl;
            if(articulated_rigid_body->muscle_activations.Valid_Index(muscle_id)) output_stream << "Activation = " << articulated_rigid_body->muscle_activations(muscle_id) << std::endl;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_MUSCLE_3D<T>::
Initialize(ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body_input,std::string basedir,int frame)
{
    if(articulated_rigid_body) return;
    opengl_triangulated_surface.Delete_Pointers_And_Clean_Memory();surface_muscle_indices.Resize(0);
    articulated_rigid_body=articulated_rigid_body_input;

    std::string muscle_info_file=STRING_UTILITIES::string_sprintf("%s/%d/muscle_info",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(muscle_info_file)) Read_Muscle_Internal_Particles(muscle_info_file);

    for(int i=1;i<=articulated_rigid_body->muscle_list->muscles.m;i++)
        for(int j=1;j<=articulated_rigid_body->muscle_list->muscles(i)->muscle_segments.m;j++){
            MUSCLE_SEGMENT<TV>* segment=articulated_rigid_body->muscle_list->muscles(i)->muscle_segments(j);
            if(segment->segment_type==MUSCLE_SEGMENT<TV>::ANALYTIC_SURFACE_SEGMENT){
                TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
                surface->Initialize_Cylinder_Mesh_And_Particles(length_resolution,radial_resolution,20,5);
                OPENGL_MATERIAL material=OPENGL_MATERIAL::Plastic(OPENGL_COLOR(1,0.6f,0.9f));
                OPENGL_TRIANGULATED_SURFACE<T>* opengl_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*surface,false,material,material);
                opengl_triangulated_surface.Append(opengl_surface);
                surface_muscle_indices.Append(PAIR<int,int>(i,j));}}
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_MUSCLE_3D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Empty_Box();
}
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_MUSCLE_SURFACE_3D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Empty_Box();
}
//#####################################################################
template class OPENGL_COMPONENT_MUSCLE_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_MUSCLE_3D<double>;
#endif
