//#####################################################################
// Copyright 2002,2003,Neil Molino,Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace OPENGL_SHAPES
//#####################################################################
#ifndef _OPENGL_SHAPES_h
#define _OPENGL_SHAPES_h

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>

namespace PhysBAM{
namespace OPENGL_SHAPES{

template<class TV>
inline void Draw_Dot(const TV& v,OPENGL_COLOR color = OPENGL_COLOR(1,1,0),const float size=5)
{   glPushAttrib(GL_LIGHTING_BIT | GL_POINT_BIT | GL_TEXTURE_BIT); 
        glDisable(GL_TEXTURE_2D); glDisable(GL_LIGHTING); glPointSize(size);
        ARRAY<typename OPENGL_POLICY<typename TV::SCALAR>::T_GL> vertices; color.Send_To_GL_Pipeline(); OpenGL_Vertex(v,vertices); OpenGL_Draw_Arrays(GL_POINTS, TV::dimension, vertices);
    glPopAttrib();
}

template<class C> 
inline void Draw_Dots(const C& dots,const int n,OPENGL_COLOR color = OPENGL_COLOR(1,1,0),double size=5)
{   glPushAttrib(GL_LIGHTING_BIT | GL_POINT_BIT | GL_TEXTURE_BIT); 
        glDisable(GL_TEXTURE_2D); glDisable(GL_LIGHTING); glPointSize(size);
        ARRAY<typename OPENGL_POLICY<typename C::SCALAR>::T_GL> vertices; color.Send_To_GL_Pipeline(); for(int i=0;i<n;++i) OpenGL_Vertex(dots[i],vertices); OpenGL_Draw_Arrays(GL_POINTS, C::dimension, vertices);
    glPopAttrib();
}

template<class TV>
void Draw_Segment(const TV& v1,const TV& v2,const OPENGL_COLOR& color=OPENGL_COLOR::Yellow(),const double line_width=1)
{   glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT); 
        glDisable(GL_TEXTURE_2D); glDisable(GL_LIGHTING); glLineWidth((GLfloat)line_width); color.Send_To_GL_Pipeline();
        ARRAY<typename OPENGL_POLICY<typename TV::SCALAR>::T_GL> vertices; OpenGL_Line(v1,v2,vertices); OpenGL_Draw_Arrays(GL_LINES, TV::dimension, vertices);
    glPopAttrib();
}

template<class T>
void Draw_Vector(const VECTOR<T,2>& from,const VECTOR<T,2>& v,OPENGL_COLOR color=OPENGL_COLOR(1,1,0),const float size=5);

template<class T>
void Draw_Vector(const VECTOR<T,3>& from,const VECTOR<T,3>& v,OPENGL_COLOR color=OPENGL_COLOR(1,1,0),const float size=5);

inline void Draw_Quad(VECTOR<double,3> corner,VECTOR<double,3> width,VECTOR<double,3> height)
{
    OpenGL_Normal((VECTOR<double,3>::Cross_Product(width,height)).Normalized());
    ARRAY<OPENGL_POLICY<double>::T_GL> vertices;ARRAY<GLfloat> textures;
    OpenGL_Texture(VECTOR<float,2>(0,0),textures);OpenGL_Vertex(corner,vertices);
    OpenGL_Texture(VECTOR<float,2>(1,0),textures);OpenGL_Vertex((corner+width),vertices);
    OpenGL_Texture(VECTOR<float,2>(0,1),textures);OpenGL_Vertex((corner+height),vertices);
    OpenGL_Texture(VECTOR<float,2>(1,1),textures);OpenGL_Vertex((corner+width+height),vertices);
    OpenGL_Draw_Arrays_With_Textures(GL_TRIANGLE_STRIP,3,vertices,textures);
}

inline void Draw_Circle(const double radius,const int slices,const bool fill=true)
{
    ARRAY<OPENGL_POLICY<double>::T_GL> vertices;
        for(int i=0;i<slices;++i){
            double t=2*pi*double(i)/slices; 
            OpenGL_Vertex(VECTOR<double,3>(radius*cos(t),radius*sin(t),0),vertices);}
    OpenGL_Draw_Arrays(fill?GL_TRIANGLE_STRIP:GL_LINE_LOOP,2,vertices);
}

template<class T>
inline void Draw_Circle(const VECTOR<T,2>& center,const double radius,const int slices,const bool fill=true)
{
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
        for(int i=0;i<slices;++i){
            T t=2*(T)pi*i/slices;
            OpenGL_Vertex(VECTOR<T,3>((T)radius*cos(t)+center.x,(T)radius*sin(t)+center.y,0),vertices);}
    OpenGL_Draw_Arrays(fill?GL_TRIANGLE_STRIP:GL_LINE_LOOP,2,vertices);
}

template<class T>
inline void Draw_Box(VECTOR<T,3> corner,VECTOR<T,3> width,VECTOR<T,3> height,VECTOR<T,3> depth)
{   VECTOR<T,3> v000=corner,
        v001=corner+width,
        v010=corner+height,
        v011=corner+width+height,
        v100=corner+depth,
        v101=corner+depth+width,
        v110=corner+depth+height,
        v111=corner+depth+height+width;
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
        OpenGL_Line(v000,v001,vertices);
        OpenGL_Line(v010,v011,vertices);
        OpenGL_Line(v100,v101,vertices);
        OpenGL_Line(v110,v111,vertices);
        
        OpenGL_Line(v000,v010,vertices);
        OpenGL_Line(v001,v011,vertices);
        OpenGL_Line(v100,v110,vertices);
        OpenGL_Line(v101,v111,vertices);

        OpenGL_Line(v000,v100,vertices);
        OpenGL_Line(v001,v101,vertices);
        OpenGL_Line(v010,v110,vertices);
        OpenGL_Line(v011,v111,vertices);
    OpenGL_Draw_Arrays(GL_LINES,3,vertices);
}

template<class T>
inline void Draw_Box(const RANGE<VECTOR<T,3> >& box)
{   VECTOR<T,3> corner(box.min_corner.x,box.min_corner.y,box.min_corner.z),x_width(box.max_corner.x-box.min_corner.x,0,0),
            y_width(0,box.max_corner.y-box.min_corner.y,0),z_width(0,0,box.max_corner.z-box.min_corner.z);
    Draw_Box(corner,x_width,y_width,z_width);
}

inline void Draw_Translucent_Stripe(int x0,int y0,int dx,int dy,const OPENGL_COLOR &color = OPENGL_COLOR::White())
{   glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING); glDisable(GL_TEXTURE_2D); glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        color.Send_To_GL_Pipeline();
        ARRAY<GLshort> vertices;
        OpenGL_Vertex(VECTOR<int,2>(x0,y0),vertices);
        OpenGL_Vertex(VECTOR<int,2>(x0+dx,y0),vertices);
        OpenGL_Vertex(VECTOR<int,2>(x0,y0+dx),vertices);
        OpenGL_Vertex(VECTOR<int,2>(x0+dx,y0+dx),vertices);
        OpenGL_Draw_Arrays(GL_TRIANGLE_STRIP,2,vertices);
    glPopAttrib();
}

//#####################################################################
// Function Draw_Arrow
//#####################################################################
// Assumes OpenGL_Begin(GL_LINES) and attributes are already set
// size is length of arrowhead as fraction of length of arrow line length
// angle is measured from arrow line
#ifndef USE_OPENGLES
template<class T>
inline void Draw_Arrow(const VECTOR<T,2>& startpt,const VECTOR<T,2>& endpt,T arrowhead_size=OPENGL_PREFERENCES::arrowhead_size,
                       T sin_angle=sin(OPENGL_PREFERENCES::arrowhead_angle),T cos_angle=cos(OPENGL_PREFERENCES::arrowhead_angle))
{
    OpenGL_Line(startpt,endpt);
    VECTOR<T,2> direction=endpt-startpt,p=endpt-(cos_angle*arrowhead_size)*direction,dp=(sin_angle*arrowhead_size)*direction.Rotate_Clockwise_90();
    OpenGL_Line(p+dp,endpt);
    OpenGL_Line(p-dp,endpt);
}
#endif

template<class T>
inline void Draw_Arrow(const VECTOR<T,2>& startpt,const VECTOR<T,2>& endpt,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices,T arrowhead_size=OPENGL_PREFERENCES::arrowhead_size,
                       T sin_angle=sin(OPENGL_PREFERENCES::arrowhead_angle),T cos_angle=cos(OPENGL_PREFERENCES::arrowhead_angle))
{
    OpenGL_Line(startpt,endpt,vertices);
    VECTOR<T,2> direction=endpt-startpt,p=endpt-(cos_angle*arrowhead_size)*direction,dp=(sin_angle*arrowhead_size)*direction.Rotate_Clockwise_90();
    OpenGL_Line(p+dp,endpt,vertices);
    OpenGL_Line(p-dp,endpt,vertices);
}
}
}
#endif

