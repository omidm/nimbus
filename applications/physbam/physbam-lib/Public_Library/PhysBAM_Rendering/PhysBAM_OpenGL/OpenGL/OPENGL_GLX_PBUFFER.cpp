//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// pbuffers code taken from pbdemo.c and pbutil.{h,c} in
// http://cvs.sourceforge.net/viewcvs.py/mesa3d/Mesa/xdemos/
// written by Brain Paul (April 1997; Updated on 5 October 2002)
//#####################################################################
#ifdef __linux__
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GLX_PBUFFER.h>
#include <cstdio>
#include <cstring>
#include <iostream>
using namespace PhysBAM;

#ifdef OPENGL_GLX_PBUFFER_SUPPORTED

namespace 
{

    static int XErrorFlag = 0;
    static int HandleXError(Display *dpy, XErrorEvent *event)
    {
        XErrorFlag = 1;
        return 0;
    }

    GLXPbuffer CreatePbuffer(Display *dpy, GLXFBConfig fbConfig, int *pbAttribs)
    {
        int (*oldHandler)(Display *, XErrorEvent *);
        GLXPbuffer pBuffer = None;

        /* Catch X protocol errors with our own error handler */
        oldHandler = XSetErrorHandler(HandleXError);

        XErrorFlag = 0;
        pBuffer = glXCreatePbuffer(dpy, fbConfig, pbAttribs);

        /* Restore original X error handler */
        (void) XSetErrorHandler( oldHandler );

        /* Return pbuffer (may be None) */
        if (!XErrorFlag && pBuffer!=None) return pBuffer;
        else return None;
    }

    /*
     * Print parameters for a GLXFBConfig to stdout.
     * Input:  dpy - the X display
     *         fbConfig - the fbconfig handle
     *         horizFormat - if true, print in horizontal format
     */
    void PrintFBConfigInfo(Display *dpy, GLXFBConfig fbConfig, int width, int height, Bool horizFormat)
    {
        int pbAttribs[] = {
            GLX_PBUFFER_WIDTH, width,
            GLX_PBUFFER_HEIGHT, height,
            GLX_LARGEST_PBUFFER, True,
            GLX_PRESERVED_CONTENTS, False,
            None};
        GLXPbuffer pBuffer;
        int bufferSize, level, doubleBuffer, stereo, auxBuffers;
        int redSize, greenSize, blueSize, alphaSize;
        int depthSize, stencilSize;
        int accumRedSize, accumBlueSize, accumGreenSize, accumAlphaSize;
        int drawableType, renderType, xRenderable, xVisual, id;
        int maxWidth, maxHeight, maxPixels;

        glXGetFBConfigAttrib(dpy, fbConfig, GLX_BUFFER_SIZE, &bufferSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_LEVEL, &level);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_DOUBLEBUFFER, &doubleBuffer);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_STEREO, &stereo);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_AUX_BUFFERS, &auxBuffers);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_RED_SIZE, &redSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_GREEN_SIZE, &greenSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_BLUE_SIZE, &blueSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_ALPHA_SIZE, &alphaSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_DEPTH_SIZE, &depthSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_STENCIL_SIZE, &stencilSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_ACCUM_RED_SIZE, &accumRedSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_ACCUM_GREEN_SIZE, &accumGreenSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_ACCUM_BLUE_SIZE, &accumBlueSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_ACCUM_ALPHA_SIZE, &accumAlphaSize);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_DRAWABLE_TYPE, &drawableType);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_RENDER_TYPE, &renderType);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_X_RENDERABLE, &xRenderable);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_X_VISUAL_TYPE_EXT, &xVisual);
        if (!xRenderable || !(drawableType & GLX_WINDOW_BIT)) xVisual = -1;
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_FBCONFIG_ID, &id);

        glXGetFBConfigAttrib(dpy, fbConfig, GLX_MAX_PBUFFER_WIDTH, &maxWidth);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_MAX_PBUFFER_HEIGHT, &maxHeight);
        glXGetFBConfigAttrib(dpy, fbConfig, GLX_MAX_PBUFFER_PIXELS, &maxPixels);

        pBuffer = CreatePbuffer(dpy, fbConfig, pbAttribs);

        if(horizFormat){
            std::stringstream ss;
            ss<<id<<" ";
            if(xVisual==GLX_STATIC_GRAY) ss<<"StaticGray  ";
            else if(xVisual==GLX_GRAY_SCALE) ss<<"GrayScale   ";
            else if(xVisual==GLX_STATIC_COLOR) ss<<"StaticColor ";
            else if(xVisual==GLX_PSEUDO_COLOR) ss<<"PseudoColor ";
            else if(xVisual==GLX_TRUE_COLOR) ss<<"TrueColor   ";
            else if(xVisual==GLX_DIRECT_COLOR) ss<<"DirectColor ";
            else ss<<"  -none-    ";
            ss<<" "<<bufferSize<<" "<<level<<" "<<(renderType&GLX_RGBA_BIT?'y':'n')<<(renderType&GLX_COLOR_INDEX_BIT?'y':'n')<<(doubleBuffer?"y":"n")<<(stereo?"y":"n")<<"  ";
            ss<<redSize<<" "<<greenSize<<" "<<blueSize<<" "<<alphaSize<<"  ";
            ss<<depthSize<<" "<<stencilSize<<"  ";
            ss<<accumRedSize<<" "<<accumGreenSize<<" "<<accumBlueSize<<" "<<accumAlphaSize;
            ss<<"       "<<(pBuffer?"y":"n")<<std::endl;
            LOG::filecout(ss.str());
        }
        else{
            std::stringstream ss;
            ss<<"Id "<<id<<std::endl;
            ss<<"  Buffer Size: "<<bufferSize<<std::endl;
            ss<<"  Level: "<<level<<std::endl;
            ss<<"  Double Buffer: "<<(doubleBuffer?"yes":"no")<<std::endl;
            ss<<"  Stereo: "<<(stereo?"yes":"no")<<std::endl;
            ss<<"  Aux Buffers: "<<auxBuffers<<std::endl;
            ss<<"  Red Size: "<<redSize<<std::endl;
            ss<<"  Green Size: "<<greenSize<<std::endl;
            ss<<"  Blue Size: "<<blueSize<<std::endl;
            ss<<"  Alpha Size: "<<alphaSize<<std::endl;
            ss<<"  Depth Size: "<<depthSize<<std::endl;
            ss<<"  Stencil Size: "<<stencilSize<<std::endl;
            ss<<"  Accum Red Size: "<<accumRedSize<<std::endl;
            ss<<"  Accum Green Size: "<<accumGreenSize<<std::endl;
            ss<<"  Accum Blue Size: "<<accumBlueSize<<std::endl;
            ss<<"  Accum Alpha Size: "<<accumAlphaSize<<std::endl;
            ss<<"  Drawable Types: ";
            if(drawableType & GLX_WINDOW_BIT) ss<<"Window ";
            if(drawableType & GLX_PIXMAP_BIT) ss<<"Pixmap ";
            if(drawableType & GLX_PBUFFER_BIT) ss<<"PBuffer";
            ss<<std::endl;
            ss<<"  Render Types: ";
            if(renderType & GLX_RGBA_BIT) ss<<"RGBA ";
            if(renderType & GLX_COLOR_INDEX_BIT) ss<<"CI ";
            ss<<std::endl;
            ss<<"  X Renderable: "<<(xRenderable?"yes":"no")<<std::endl;
            ss<<"  Pbuffer: "<<(pBuffer?"yes":"no")<<std::endl;
            LOG::filecout(ss.str());
        }

        if(pBuffer) glXDestroyPbuffer(dpy, pBuffer);
    }

    static GLXPbuffer MakePbuffer(Display *dpy, int screen, int width, int height, bool verbose, GLXFBConfig &fbconfig)
    {
        const int NUM_FB_CONFIGS = 2;
        char fbString[NUM_FB_CONFIGS][100] = {
            "Single Buffered, depth buffer",
            "Double Buffered, depth buffer"
        };
        int fbAttribs[NUM_FB_CONFIGS][100] = {
        {
            /* Single buffered, with depth buffer */
            GLX_RENDER_TYPE, GLX_RGBA_BIT,
            GLX_DRAWABLE_TYPE, GLX_PIXMAP_BIT,
            GLX_RED_SIZE, 1,
            GLX_GREEN_SIZE, 1,
            GLX_BLUE_SIZE, 1,
            GLX_DEPTH_SIZE, 1,
            GLX_DOUBLEBUFFER, 0,
            GLX_STENCIL_SIZE, 0,
            None
        },
        {
            /* Double buffered, with depth buffer */
            GLX_RENDER_TYPE, GLX_RGBA_BIT,
            GLX_DRAWABLE_TYPE, GLX_PIXMAP_BIT,
            GLX_RED_SIZE, 1,
            GLX_GREEN_SIZE, 1,
            GLX_BLUE_SIZE, 1,
            GLX_DEPTH_SIZE, 1,
            GLX_DOUBLEBUFFER, 1,
            GLX_STENCIL_SIZE, 0,
            None
        }
        };
        int pbAttribs[] = {
            GLX_PBUFFER_WIDTH, width,
            GLX_PBUFFER_HEIGHT, height,
            GLX_LARGEST_PBUFFER, True,
            GLX_PRESERVED_CONTENTS, False,
            None};
        GLXFBConfig *fbConfigs;
        GLXPbuffer pBuffer = None;
        int nConfigs;
        int i;
        int attempt;

        for (attempt=0; attempt<NUM_FB_CONFIGS; attempt++) {

            /* Get list of possible frame buffer configurations */
            fbConfigs = glXChooseFBConfig(dpy, screen, fbAttribs[attempt], &nConfigs);
            if (nConfigs==0 || !fbConfigs) {
                std::stringstream ss;
                ss<<"Error: glxChooseFBConfig failed"<<std::endl;
                LOG::filecout(ss.str());
                XCloseDisplay(dpy);
                return 0;
            }

            if (verbose) 
                for (i=0;i<nConfigs;i++) {
                    std::stringstream ss;
                    ss<<"Config "<<i<<std::endl;
                    LOG::filecout(ss.str());
                    PrintFBConfigInfo(dpy, fbConfigs[i], width, height, 0);
                }

            /* Create the pbuffer using first fbConfig in the list that works. */
            for (i=0;i<nConfigs;i++) {
                int alphaSize;
                glXGetFBConfigAttrib(dpy, fbConfigs[i], GLX_ALPHA_SIZE, &alphaSize);
                if(alphaSize<8) continue; // need alpha buffer
                pBuffer = CreatePbuffer(dpy, fbConfigs[i], pbAttribs);
                if (pBuffer) {
                    fbconfig = fbConfigs[i];
                    if (verbose){ 
                        std::stringstream ss;
                        ss<<"Picked config "<<i<<std::endl;
                        LOG::filecout(ss.str());
                    }
                    break;
                }
            }

            if (pBuffer!=None) break;
        }

        if (pBuffer) {
            std::stringstream ss;
            ss<<"Using: "<<fbString[attempt]<<std::endl;
            LOG::filecout(ss.str());
        }

        XFree(fbConfigs);

        return pBuffer;
    }
}

OPENGL_PBUFFER::OPENGL_PBUFFER()
    :verbose(false),display(0),fbconfig(0),pbuffer(0)
{}

OPENGL_PBUFFER::~OPENGL_PBUFFER()
{
    Destroy();
}

bool OPENGL_PBUFFER::
Create(int width, int height)
{
    if (pbuffer) {
        LOG::cerr << "Destroying old pbuffer before creating new one" << std::endl;
        Destroy();
    }

    /* Open the X display */
    display = XOpenDisplay(0);
    if (!display) {
        std::stringstream ss;
        ss<<"Error: couldn't open default X display."<<std::endl; 
        LOG::filecout(ss.str());
        return false; 
    }

    /* Get default screen */
    int screen = DefaultScreen(display);

    /* Test that pbuffers are available */
    char *extensions;
    extensions = (char *)glXQueryServerString(display, screen, GLX_EXTENSIONS);
    if (!strstr(extensions,"GLX_SGIX_fbconfig") || !strstr(extensions,"GLX_SGIX_pbuffer")) {
        std::stringstream ss;
        ss<<"Error: pbuffers not available on this screen"<<std::endl;
        LOG::filecout(ss.str());
        XCloseDisplay(display);
        return false;
    }

    /* Create Pbuffer */
    pbuffer = MakePbuffer(display, screen, width, height, verbose, fbconfig);
    if (pbuffer == None) {
        std::stringstream ss;
        ss<<"Error: couldn't create pbuffer"<<std::endl;
        LOG::filecout(ss.str());
        XCloseDisplay(display);
        return false;
    }

    /* Get corresponding XVisualInfo */
    XVisualInfo *vis_info = glXGetVisualFromFBConfig(display, fbconfig);
    if (!vis_info) {
        std::stringstream ss;
        ss<<"Error: can't get XVisualInfo from FBconfig"<<std::endl;
        LOG::filecout(ss.str());
        XCloseDisplay(display);
        return false;
    }

    /* Create GLX context */
    GLXContext glx_context = glXCreateContext(display, vis_info, 0, True);
    if (!glx_context) {
        /* try indirect */
        glx_context = glXCreateContext(display, vis_info, 0, False);
        std::stringstream ss;
        if (!glx_context) {
            ss<<"Error: Couldn't create GLXContext"<<std::endl;
            LOG::filecout(ss.str());
            XFree(vis_info);
            XCloseDisplay(display);
            return false;
        }
        else {
            ss<<"Warning: using indirect GLXContext"<<std::endl;
            LOG::filecout(ss.str());
        }
    }

    /* Bind context to pbuffer */
    if (!glXMakeCurrent(display, pbuffer, glx_context)) {
        std::stringstream ss;
        ss<<"Error: glXMakeCurrent failed"<<std::endl;
        LOG::filecout(ss.str());
        XFree(vis_info);
        XCloseDisplay(display);
        return false;
    }

    return true;  /* Success!! */
}

void OPENGL_PBUFFER::
Destroy()
{
    if(pbuffer){
        if(verbose) {
            std::stringstream ss;
            ss<<"Destroying pbuffer"<<std::endl;
            LOG::filecout(ss.str());
        }
        glXDestroyPbuffer(display, pbuffer);}
}

#endif
#endif
