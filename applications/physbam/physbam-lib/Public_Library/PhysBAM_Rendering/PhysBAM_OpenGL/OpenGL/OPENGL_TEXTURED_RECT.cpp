//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TEXTURED_RECT.h>
#include <iostream>

using namespace PhysBAM;
using namespace std;

OPENGL_TEXTURED_RECT::OPENGL_TEXTURED_RECT() : texture(0)
{
}

void
OPENGL_TEXTURED_RECT::Set_Texture(OPENGL_TEXTURE *texture_input)
{
    texture = texture_input;
}

void
OPENGL_TEXTURED_RECT::Display(const int in_color) const
{
    if (!texture) return;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    // Force it into fill mode (e.g. it's in wireframe mode)
#ifndef USE_OPENGLES
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
#endif

    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glBindTexture(GL_TEXTURE_2D, texture->id);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);

    ARRAY<OPENGL_POLICY<double>::T_GL> vertices;ARRAY<float> textures;
    OpenGL_Texture(VECTOR<float,2>((float)(texture->min_s),(float)(texture->min_t)),textures);
    OpenGL_Vertex(VECTOR<double,2>(-0.5*width,-0.5*height),vertices);
    OpenGL_Texture(VECTOR<float,2>((float)(texture->min_s),(float)(texture->max_t)),textures);
    OpenGL_Vertex(VECTOR<double,2>(-0.5*width,0.5*height),vertices);
    OpenGL_Texture(VECTOR<float,2>((float)(texture->max_s),(float)(texture->min_t)),textures);
    OpenGL_Vertex(VECTOR<double,2>(0.5*width,-0.5*height),vertices);
    OpenGL_Texture(VECTOR<float,2>((float)(texture->max_s),(float)(texture->max_t)),textures);
    OpenGL_Vertex(VECTOR<double,2>(0.5*width,0.5*height),vertices);
    OpenGL_Draw_Arrays_With_Textures(GL_TRIANGLE_STRIP,2,vertices,textures);

#ifdef USE_OPENGLES
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
#endif

    glPopAttrib();

    glPopMatrix();
}

RANGE<VECTOR<float,3> >
OPENGL_TEXTURED_RECT::Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(-0.5f*(float)width,0.5f*(float)width,-0.5f*(float)height,0.5f*(float)height,0,0));
}
