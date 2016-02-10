//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Robert Bridson, Sergey Koltakov, Neil Molino, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace OPENGL_PRIMITIVES
//#####################################################################
#ifndef __OPENGL_PRIMITIVES__
#define __OPENGL_PRIMITIVES__
#include <cstdio>
#ifdef _WIN32
#include <windows.h>
#endif

#ifndef __APPLE__
#ifndef USE_OPENGLES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#else
#include <GLES/gl.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/glues.h>
#endif
#else
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POLICY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_EPS_OUTPUT.h>
namespace PhysBAM{

#define WANT_OPENGL_EPS_OUTPUT // Uncomment this if you want support for dumping eps.
extern OPENGL_EPS_OUTPUT<float>* opengl_eps_output;
#ifdef WANT_OPENGL_EPS_OUTPUT
#define IF_OPENGL_EPS_OUTPUT(x) if(opengl_eps_output){x;}
#else
#define IF_OPENGL_EPS_OUTPUT(x)
#endif

#ifdef USE_OPENGLES
#define glTranslated(x,y,z)     glTranslatef((GLfloat)x,(GLfloat)y,(GLfloat)z)
#define glScaled(x,y,z)         glScalef((GLfloat)x,(GLfloat)y,(GLfloat)z)
#define glRotated(angle,x,y,z)  glRotatef((GLfloat)angle,(GLfloat)x,(GLfloat)y,(GLfloat)z)
#define glNormal3d(x,y,z)       glNormal3f((GLfloat)x,(GLfloat)y,(GLfloat)z)
#define glColor4fv(v)           glColor4f(((GLfloat*)v)[0], ((GLfloat*)v)[1], ((GLfloat*)v)[2], ((GLfloat*)v)[3])
#define glColor3f(r,g,b)        glColor4f(r,g,b,1)
#define GLdouble                GLfloat
#define glGetDoublev            glGetFloatv
#define glPushAttrib(x)
#define glPopAttrib()
#endif

inline void OpenGL_Draw_Arrays_With_Normals(GLenum mode,GLenum type,int dimension,int length,const GLvoid* vertices,const GLvoid* normals)
{glEnableClientState(GL_VERTEX_ARRAY);glEnableClientState(GL_NORMAL_ARRAY);glVertexPointer(dimension,type,0,vertices);glNormalPointer(GL_FLOAT,0,normals);glDrawArrays(mode,0,length);glDisableClientState(GL_VERTEX_ARRAY);glDisableClientState(GL_NORMAL_ARRAY);}

inline void OpenGL_Draw_Arrays_With_Textures(GLenum mode,GLenum type,int dimension,int length,const GLvoid* vertices,const GLvoid* textures)
{glEnableClientState(GL_VERTEX_ARRAY);glEnableClientState(GL_TEXTURE_COORD_ARRAY);glVertexPointer(dimension,type,0,vertices);glTexCoordPointer(2,GL_FLOAT,0,textures);glDrawArrays(mode,0,length);glDisableClientState(GL_VERTEX_ARRAY);glDisableClientState(GL_TEXTURE_COORD_ARRAY);}

inline void OpenGL_Draw_Arrays(GLenum mode,GLenum type,int dimension,int length,const GLvoid* vertices,const GLvoid* colors)
{glEnableClientState(GL_VERTEX_ARRAY);glEnableClientState(GL_COLOR_ARRAY);glVertexPointer(dimension,type,0,vertices);glColorPointer(4,GL_FLOAT,0,colors);glDrawArrays(mode,0,length);glDisableClientState(GL_VERTEX_ARRAY);glDisableClientState(GL_COLOR_ARRAY);}

inline void OpenGL_Draw_Arrays(GLenum mode,GLenum type,int dimension,int length,const GLvoid* vertices)
{glEnableClientState(GL_VERTEX_ARRAY);glVertexPointer(dimension,type,0,vertices);glDrawArrays(mode,0,length);glDisableClientState(GL_VERTEX_ARRAY);}

#ifndef USE_OPENGLES
inline void OpenGL_Draw_Arrays_With_Normals(GLenum mode,int dimension,const ARRAY<GLdouble>& vertices,const ARRAY<GLfloat>& normals)
{assert(vertices.m/dimension==normals.m/3);OpenGL_Draw_Arrays_With_Normals(mode,GL_DOUBLE,dimension,vertices.m/dimension,vertices.base_pointer,normals.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Draw_Arrays_With_Textures(GLenum mode,int dimension,const ARRAY<GLdouble>& vertices,const ARRAY<GLfloat>& textures)
{assert(vertices.m/dimension==textures.m/2);OpenGL_Draw_Arrays_With_Textures(mode,GL_DOUBLE,dimension,vertices.m/dimension,vertices.base_pointer,textures.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Draw_Arrays(GLenum mode,int dimension,const ARRAY<GLdouble>& vertices,const ARRAY<GLfloat>& colors)
{assert(vertices.m/dimension==colors.m/4);OpenGL_Draw_Arrays(mode,GL_DOUBLE,dimension,vertices.m/dimension,vertices.base_pointer,colors.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Draw_Arrays(GLenum mode,int dimension,const ARRAY<GLdouble>& vertices)
{OpenGL_Draw_Arrays(mode,GL_DOUBLE,dimension,vertices.m/dimension,vertices.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}
#endif

inline void OpenGL_Draw_Arrays_With_Normals(GLenum mode,int dimension,const ARRAY<GLfloat>& vertices,const ARRAY<GLfloat>& normals)
{assert(vertices.m/dimension==normals.m/3);OpenGL_Draw_Arrays_With_Normals(mode,GL_FLOAT,dimension,vertices.m/dimension,vertices.base_pointer,normals.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Draw_Arrays_With_Textures(GLenum mode,int dimension,const ARRAY<GLfloat>& vertices,const ARRAY<GLfloat>& textures)
{assert(vertices.m/dimension==textures.m/2);OpenGL_Draw_Arrays_With_Textures(mode,GL_FLOAT,dimension,vertices.m/dimension,vertices.base_pointer,textures.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Draw_Arrays(GLenum mode,int dimension,const ARRAY<GLfloat>& vertices,const ARRAY<GLfloat>& colors)
{assert(vertices.m/dimension==colors.m/4);OpenGL_Draw_Arrays(mode,GL_FLOAT,dimension,vertices.m/dimension,vertices.base_pointer,colors.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Draw_Arrays(GLenum mode,int dimension,const ARRAY<GLfloat>& vertices)
{OpenGL_Draw_Arrays(mode,GL_FLOAT,dimension,vertices.m/dimension,vertices.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Draw_Arrays(GLenum mode,int dimension,const ARRAY<GLshort>& vertices)
{OpenGL_Draw_Arrays(mode,GL_SHORT,dimension,vertices.m/dimension,vertices.base_pointer);
IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex_Array(dimension,vertices));IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

#ifndef USE_OPENGLES
inline void OpenGL_Begin(GLenum mode)
{glBegin(mode);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));}

inline void OpenGL_End()
{glEnd();IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Eps_Emit(const char* str)
{IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Emit(str));}
#endif

inline void OpenGL_Rotate(const ROTATION<VECTOR<float,3> >& r)
{float angle;VECTOR<float,3> axis;r.Get_Angle_Axis(angle,axis);glRotatef(angle*180/(float)pi,axis.x,axis.y,axis.z);}

inline void OpenGL_Rotate(const ROTATION<VECTOR<double,3> >& r)
{double angle;VECTOR<double,3> axis;r.Get_Angle_Axis(angle,axis);glRotated(angle*180/pi,axis.x,axis.y,axis.z);}

inline void OpenGL_Translate(const VECTOR<float,2>& v)
{glTranslatef(v.x,v.y,(float)0);}

inline void OpenGL_Translate(const VECTOR<double,2>& v)
{glTranslated(v.x,v.y,(double)0);}

inline void OpenGL_Translate(const VECTOR<float,3>& v)
{glTranslatef(v.x,v.y,v.z);}

inline void OpenGL_Translate(const VECTOR<double,3>& v)
{glTranslated(v.x,v.y,v.z);}

inline void OpenGL_Scale(const VECTOR<float,3>& v)
{glScalef(v.x,v.y,v.z);}

inline void OpenGL_Scale(const VECTOR<double,3>& v)
{glScaled(v.x,v.y,v.z);}

template<class T>
inline void OpenGL_Transform(const FRAME<VECTOR<T,3> >& frame)
{OpenGL_Translate(frame.t);OpenGL_Rotate(frame.r.Normalized());}

inline void OpenGL_Texture(const VECTOR<float,2>& texture,ARRAY<GLfloat>& textures)
{textures.Append(texture[1]);textures.Append(texture[2]);}

inline void OpenGL_Color(const GLfloat* color,ARRAY<GLfloat>& colors)
{colors.Append(color[0]);colors.Append(color[1]);colors.Append(color[2]);colors.Append(color[3]);}

template<int d>
inline void OpenGL_Vertex(const VECTOR<int,d>& v,ARRAY<GLshort>& vertices)
{for(int i=1;i<=d;i++) vertices.Append(v(i));}

template<int d>
inline void OpenGL_Vertex(const VECTOR<float,d>& v,ARRAY<GLfloat>& vertices)
{for(int i=1;i<=d;i++) vertices.Append(v(i));}

template<int d>
inline void OpenGL_Normal(const VECTOR<float,d>& n, ARRAY<GLfloat>& normals)
{for(int i=1;i<=d;i++) normals.Append(n(i));}

#ifndef USE_OPENGLES
inline void OpenGL_Vertex(const VECTOR<float,3>& v)
{glVertex3f(v.x,v.y,v.z);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<float,2>& v)
{glVertex2f(v.x,v.y);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<float,1>& v)
{glVertex2f(v.x,(float)0);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<float,2>(v.x,0)));}
#endif

inline void OpenGL_Normal(const VECTOR<float,3>& n)
{glNormal3f(n.x,n.y,n.z);}

#ifndef USE_OPENGLES
inline void OpenGL_RasterPos(const VECTOR<float,1>& v)
{glRasterPos2f(v.x,(float)0);}

inline void OpenGL_RasterPos(const VECTOR<float,2>& v)
{glRasterPos2f(v.x,v.y);}

inline void OpenGL_RasterPos(const VECTOR<float,3>& v)
{glRasterPos3f(v.x,v.y,v.z);}
#endif

template<int d>
inline void OpenGL_Vertex(const VECTOR<double,d>& v,ARRAY<OPENGL_POLICY<double>::T_GL>& vertices)
{for(int i=1;i<=d;i++) vertices.Append(v(i));}

template<int d>
inline void OpenGL_Normal(const VECTOR<double,d>& n, ARRAY<GLfloat>& normals)
{for(int i=1;i<=d;i++) normals.Append((GLfloat)n(i));}

#ifndef USE_OPENGLES
inline void OpenGL_Vertex(const VECTOR<double,3>& v)
{glVertex3d(v.x,v.y,v.z);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<double,2>& v)
{glVertex2d(v.x,v.y);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<double,1>& v)
{glVertex2d(v.x,(double)0);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<double,2>(v.x,0)));}
#endif

template<class T,int d>
inline void OpenGL_Line(const VECTOR<T,d>& a, const VECTOR<T,d>& b,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices)
{OpenGL_Vertex(a,vertices);OpenGL_Vertex(b,vertices);}

#ifndef USE_OPENGLES
template<class T,int d>
inline void OpenGL_Line(const VECTOR<T,d>& a, const VECTOR<T,d>& b)
{OpenGL_Vertex(a);OpenGL_Vertex(b);}
#endif

template<class T,int d>
inline void OpenGL_Triangle(const VECTOR<T,d>& a, const VECTOR<T,d>& b, const VECTOR<T,d>& c,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices)
{OpenGL_Vertex(a,vertices);OpenGL_Vertex(b,vertices);OpenGL_Vertex(c,vertices);}

#ifndef USE_OPENGLES
template<class T,int d>
inline void OpenGL_Triangle(const VECTOR<T,d>& a, const VECTOR<T,d>& b, const VECTOR<T,d>& c)
{OpenGL_Vertex(a);OpenGL_Vertex(b);OpenGL_Vertex(c);}
#endif

inline void OpenGL_Normal(const VECTOR<double,3>& n)
{glNormal3d(n.x,n.y,n.z);}

#ifndef USE_OPENGLES
inline void OpenGL_RasterPos(const VECTOR<double,1>& v)
{glRasterPos2f((GLfloat)v.x,(GLfloat)0);}

inline void OpenGL_RasterPos(const VECTOR<double,2>& v)
{glRasterPos2d(v.x,v.y);}

inline void OpenGL_RasterPos(const VECTOR<double,3>& v)
{glRasterPos3d(v.x,v.y,v.z);}
#endif

inline void OpenGL_LookFrom(const FRAME<VECTOR<double,3> >& frame)
{OpenGL_Rotate(frame.r.Inverse().Normalized());OpenGL_Translate(-frame.t);}

#ifndef USE_OPENGLES
template<class TV>
inline void OpenGL_String(const TV& position,const std::string& str,void* font=GLUT_BITMAP_HELVETICA_12)
{OpenGL_RasterPos(position);
for(unsigned int j=0;j<str.length();j++) glutBitmapCharacter(font,str[j]);}
#endif

#ifndef USE_OPENGLES
template<class T>
inline void OpenGL_Quad_2D(const VECTOR<T,2> &bottom_left,const VECTOR<T,2> &top_right)
{OpenGL_Vertex(bottom_left);OpenGL_Vertex(VECTOR<T,2>(bottom_left.x,top_right.y));
 OpenGL_Vertex(top_right);OpenGL_Vertex(VECTOR<T,2>(top_right.x,bottom_left.y));
 IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<T,2>(bottom_left.x,bottom_left.y));opengl_eps_output->Vertex(VECTOR<T,2>(bottom_left.x,top_right.y));
     opengl_eps_output->Vertex(VECTOR<T,2>(top_right.x,top_right.y));opengl_eps_output->Vertex(VECTOR<T,2>(top_right.x,bottom_left.y)));}
#endif

template<class T>
inline void OpenGL_Triangle_Strip_2D(const VECTOR<T,2> &bottom_left,const VECTOR<T,2> &top_right,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices)
{OpenGL_Vertex(bottom_left,vertices);OpenGL_Vertex(VECTOR<T,2>(bottom_left.x,top_right.y),vertices);
 OpenGL_Vertex(VECTOR<T,2>(top_right.x,bottom_left.y),vertices);OpenGL_Vertex(top_right,vertices);}

template<class T>
inline void OpenGL_Quad_2D(const VECTOR<T,2> &bottom_left,const VECTOR<T,2> &top_right,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices)
{OpenGL_Vertex(bottom_left,vertices);OpenGL_Vertex(VECTOR<T,2>(bottom_left.x,top_right.y),vertices);
 OpenGL_Vertex(top_right,vertices);OpenGL_Vertex(VECTOR<T,2>(top_right.x,bottom_left.y),vertices);}

template<class T>
inline void OpenGL_Quad(const VECTOR<T,3> &bottom_left,const VECTOR<T,3> &right,const VECTOR<T,3>& up,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices)
{OpenGL_Vertex(bottom_left,vertices);OpenGL_Vertex(bottom_left+up,vertices);OpenGL_Vertex(bottom_left+up+right,vertices);OpenGL_Vertex(bottom_left+right,vertices);}

#ifndef USE_OPENGLES
template<class T>
inline void OpenGL_Quad(const VECTOR<T,3> &bottom_left,const VECTOR<T,3> &right,const VECTOR<T,3>& up)
{OpenGL_Vertex(bottom_left);OpenGL_Vertex(bottom_left+up);OpenGL_Vertex(bottom_left+up+right);OpenGL_Vertex(bottom_left+right);}
#endif

#ifndef USE_OPENGLES
template<class T>
inline void OpenGL_Clip_Plane(GLenum id,const PLANE<T> &plane)
{GLdouble equation[4]={plane.normal.x,plane.normal.y,plane.normal.z,-VECTOR<T,3>::Dot_Product(plane.normal,plane.x1)};
glClipPlane(id,equation);}
#endif

}
#endif
