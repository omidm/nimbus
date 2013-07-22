/* $Id: gl2psTestSimple.c,v 1.1 2006/02/01 19:05:31 irving Exp $ */
/*
 * GL2PS, an OpenGL to PostScript Printing Library
 * Copyright (C) 1999-2005 Christophe Geuzaine <geuz@geuz.org>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of either:
 *
 * a) the GNU Library General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your
 * option) any later version; or
 *
 * b) the GL2PS License as published by Christophe Geuzaine, either
 * version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See either
 * the GNU Library General Public License or the GL2PS License for
 * more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library in the file named "COPYING.LGPL";
 * if not, write to the Free Software Foundation, Inc., 675 Mass Ave,
 * Cambridge, MA 02139, USA.
 *
 * You should have received a copy of the GL2PS License with this
 * library in the file named "COPYING.GL2PS"; if not, I will be glad
 * to provide one.
 *
 * For the latest info about gl2ps, see http://www.geuz.org/gl2ps/.
 * Please report all bugs and problems to <gl2ps@geuz.org>.
 */

/*
  To compile on Linux:
  gcc gl2psTestSimple.c gl2ps.c -lglut -lGL -lGLU -L/usr/X11R6/lib -lX11 -lXi -lXmu -lm

  To compile on MacOSX:
  gcc gl2psTestSimple.c gl2ps.c -framework OpenGL -framework GLUT -framework Cocoa
*/

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include <string.h>
#include "gl2ps.h"

void display(void){
  unsigned int i;
  char *help = "Press 's' to save image or 'q' to quit";  

  glClearColor(0.3, 0.5, 0.8, 0.);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();
  glRotatef(-60., 2., 0., 1.);
  glEnable(GL_LIGHTING);
  glutSolidTorus(0.3, 0.6, 30, 30);
  glDisable(GL_LIGHTING);
  glPopMatrix();
  glColor3f(1.,1.,1.);
  glRasterPos2d(-0.9,-0.9);
  gl2psText(help, "Times-Roman", 24);
  for (i = 0; i < strlen(help); i++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, help[i]);
  glFlush();
}

void keyboard(unsigned char key, int x, int y){
  FILE *fp;
  int state = GL2PS_OVERFLOW, buffsize = 0;

  switch(key){
  case 'q':
    exit(0);
    break;
  case 's':
    fp = fopen("out.eps", "wb");
    while(state == GL2PS_OVERFLOW){
      buffsize += 1024*1024;
      gl2psBeginPage("test", "gl2psTestSimple", NULL, GL2PS_EPS, GL2PS_BSP_SORT, 
		     GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT, 
		     GL_RGBA, 0, NULL, 0, 0, 0,  buffsize, fp, "out.eps");
      display();
      state = gl2psEndPage();
    }
    fclose(fp);
    printf("Done!\n");
    break;
  }
}

int main(int argc, char **argv){
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DEPTH);
  glutInitWindowSize(400, 400);
  glutInitWindowPosition(100, 100);
  glutCreateWindow(argv[0]);

  GLfloat pos[4]={1.,1.,-1.,0.};
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_POSITION, pos);

  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);

  glutMainLoop();
  return 0;
}
