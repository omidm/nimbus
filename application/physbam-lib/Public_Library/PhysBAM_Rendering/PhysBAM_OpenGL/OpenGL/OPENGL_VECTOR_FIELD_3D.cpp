//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_VECTOR_FIELD_3D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box(World_Space_Point(VECTOR<float,3>(vector_locations(1))));
    for (int i=1;i<=vector_locations.m;i++)box.Enlarge_Nonempty_Box_To_Include_Point(World_Space_Point(VECTOR<float,3>(vector_locations(i))));
    return box;
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_VECTOR_FIELD_3D<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    if(is_moving_grid){
        T angle;VECTOR<T,3> axis;rigid_grid_frame.r.Get_Angle_Axis(angle,axis);angle=(T)(angle/pi*180);
        glTranslatef((GLfloat)rigid_grid_frame.t(1),(GLfloat)rigid_grid_frame.t(2),(GLfloat)rigid_grid_frame.t(3));glRotatef((GLfloat)angle,(GLfloat)axis(1),(GLfloat)axis(2),(GLfloat)axis(3));}
    Send_Transform_To_GL_Pipeline();

    VECTOR<T,3> head;
    glMatrixMode(GL_MODELVIEW);

    if(!draw_fancy_arrow){
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT | GL_CURRENT_BIT);
        vector_color.Send_To_GL_Pipeline();
        glLineWidth(1);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
        ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
        for(int i=1;i<=vector_locations.m;i++){
            head=vector_locations(i)+(T)size*vector_field(i);
            OpenGL_Line(vector_locations(i),head,vertices);
            if(draw_arrowhead){
                VECTOR<T,3> orth_vect=vector_field(i).Orthogonal_Vector();
                orth_vect*=(T)(.15*size);
                OpenGL_Line(head,head+orth_vect-(T).15*(T)size*vector_field(i),vertices);
                OpenGL_Line(head,head-orth_vect-(T).15*(T)size*vector_field(i),vertices);
            }
        }
        OpenGL_Draw_Arrays(GL_LINES,3,vertices);
        if(draw_basepoint){
            float old_point_size;
            (OPENGL_COLOR::Yellow(0.6f)).Send_To_GL_Pipeline();
            glGetFloatv(GL_POINT_SIZE,&old_point_size);
            glPointSize(2.0);
            vertices.Resize(0);
            for(int i=1;i<=vector_locations.m;i++)
                if(vector_field(i).Magnitude_Squared()>0) OpenGL_Vertex(vector_locations(i),vertices);
            OpenGL_Draw_Arrays(GL_POINTS,3,vertices);
            vector_color.Send_To_GL_Pipeline();
            glPointSize(old_point_size);}
#ifndef USE_OPENGLES
        if(draw_value){
            glDisable(GL_DEPTH_TEST);
            (vector_color+OPENGL_COLOR(0.8f,0.8f,0.8f)).Send_To_GL_Pipeline();
            for(int i=1;i<=vector_locations.m;i++)
                OpenGL_String(vector_locations(i)+(T)1.1*(T)size*vector_field(i),STRING_UTILITIES::string_sprintf("%.3f %.3f %.3f",vector_field(i).x,vector_field(i).y,vector_field(i).z));
            vector_color.Send_To_GL_Pipeline();
            glEnable(GL_DEPTH_TEST);}
#endif
        glPopAttrib();}
#ifndef USE_OPENGLES
    else{
        double length;
        if(!vector_hat) vector_hat=gluNewQuadric();
        glDisable(GL_CULL_FACE);
        OPENGL_MATERIAL(vector_color).Send_To_GL_Pipeline(); // arrowhead needs material
        for(int i=1;i<=vector_locations.m;i++){
            T len_squared=vector_field(i).Magnitude_Squared();
            if(len_squared>0){
                VECTOR<T,3> orthogonal_vector=vector_field(i).Orthogonal_Vector();
                orthogonal_vector.Normalize();
                head=vector_locations(i)+(T)size*vector_field(i);
                ROTATION<VECTOR<T,3> > vector_hat_orientation=ROTATION<VECTOR<T,3> >::From_Rotated_Vector(VECTOR<T,3>(0,0,1),vector_field(i));
                glPushMatrix();
                glTranslated(head.x,head.y,head.z);
                OpenGL_Rotate(vector_hat_orientation);
                length=sqrt(len_squared);
                T arrow_head_length=(T)(size*length/3);
                T arrow_head_radius=(T)(size*length/4);
                T cylinder_radius=(T)(size*length/16);
                T cylinder_length=(T)(length*size);
                const int arrow_subdivisions=4;
                glTranslated(0,0,-cylinder_length);
                gluCylinder(vector_hat,cylinder_radius,cylinder_radius,cylinder_length-arrow_head_length,arrow_subdivisions,4);
                glTranslated(0,0,cylinder_length-arrow_head_length);
                gluCylinder(vector_hat,size*length/4,0,arrow_head_length,arrow_subdivisions,4);
                gluDisk(vector_hat,0,arrow_head_radius,arrow_subdivisions,1);
                glPopMatrix();}}
         glEnable(GL_CULL_FACE);}
#endif

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> void OPENGL_VECTOR_FIELD_3D<T>::
Scale_Vector_Size(const T scale)
{
    size*=scale;
}
//#####################################################################
// Function Toggle_Arrowhead_Mode
//#####################################################################
template<class T> void OPENGL_VECTOR_FIELD_3D<T>::
Toggle_Arrowhead_Mode()
{
    if(draw_fancy_arrow){draw_fancy_arrow=false;draw_arrowhead=false;}
    else if(draw_arrowhead){draw_fancy_arrow=true;draw_arrowhead=false;}
    else{draw_fancy_arrow=false;draw_arrowhead=true;}
}
//#####################################################################
template class OPENGL_VECTOR_FIELD_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_VECTOR_FIELD_3D<double>;
#endif
