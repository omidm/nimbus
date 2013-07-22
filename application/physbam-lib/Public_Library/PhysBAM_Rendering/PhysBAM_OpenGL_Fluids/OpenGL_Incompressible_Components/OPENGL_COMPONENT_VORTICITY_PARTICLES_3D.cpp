//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Fluids/OpenGL_Incompressible_Components/OPENGL_COMPONENT_VORTICITY_PARTICLES_3D.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<T,RW>::
OPENGL_COMPONENT_VORTICITY_PARTICLES_3D(const std::string &filename,bool use_ids_input)
    :OPENGL_COMPONENT_PARTICLES_3D<T,VORTICITY_PARTICLES<VECTOR<T,3> >,RW>(filename,"",use_ids_input,false)
{}

template<class T,class RW> void OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<T,RW>::
Reinitialize(bool force)
{
    OPENGL_COMPONENT_PARTICLES_3D<T,VORTICITY_PARTICLES<VECTOR<T,3> >,RW>::Reinitialize(force);
    if(!have_velocities){
        have_velocities=true;
        opengl_vector_field.vector_field.Resize(particles->array_collection->Size());
        int idx=1;
        for(int i=1;i<=particles->array_collection->Size();i++)
            opengl_vector_field.vector_field(idx++)=particles->vorticity(i);}
}
template class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<double,double>;
#endif
