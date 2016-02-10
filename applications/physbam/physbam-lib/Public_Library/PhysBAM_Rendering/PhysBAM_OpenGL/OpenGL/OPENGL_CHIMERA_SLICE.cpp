//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_CHIMERA_SLICE
//##################################################################### 
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CHIMERA_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
using namespace PhysBAM;
//#####################################################################
// Function Print_Slice_Info
//#####################################################################
template<class T> void OPENGL_CHIMERA_SLICE<T>::
Print_Slice_Info(std::ostream& output_stream)
{
    for(int i=1;i<=opengl_uniform_slices.m;i++){
        output_stream<<"Grid "<<i<<std::endl;
        opengl_uniform_slices(i)->Print_Slice_Info(output_stream);}
}

//#####################################################################
template class OPENGL_CHIMERA_SLICE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_CHIMERA_SLICE<double>;
#endif
