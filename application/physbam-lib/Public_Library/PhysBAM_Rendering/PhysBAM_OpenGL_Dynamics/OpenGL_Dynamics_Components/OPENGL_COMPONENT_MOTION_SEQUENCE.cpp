//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Dynamics/OpenGL_Dynamics_Components/OPENGL_COMPONENT_MOTION_SEQUENCE.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <sstream>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_MOTION_SEQUENCE
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
OPENGL_COMPONENT_MOTION_SEQUENCE(const std::string& filename_input,const bool frame_dependent_data_input,const std::string& segments_filename,const T frame_rate_input,const OPENGL_COLOR& color)
    :OPENGL_COMPONENT("Motion Sequence"),filename(filename_input),segmented_curve(segment_mesh,particles),opengl_segmented_curve(segmented_curve),
     frame_loaded(-1),valid(false),frame_dependent_data(frame_dependent_data_input),frame_rate(frame_rate_input),one_over_frame_rate((T)1/frame_rate)
{
    opengl_segmented_curve.vertex_color=color;
    is_animation=true;opengl_segmented_curve.draw_vertices=true;
    std::string frame_filename=STRING_UTILITIES::string_sprintf(filename.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(frame_filename)) FILE_UTILITIES::Read_From_File<RW>(frame_filename,motion);
    if(FILE_UTILITIES::File_Exists(segments_filename)) FILE_UTILITIES::Read_From_File<RW>(segments_filename,segment_mesh.elements);
    segment_mesh.number_nodes=motion.trajectories.m;
    particles.array_collection->Add_Elements(motion.trajectories.m);
}
//#####################################################################
// Function ~OPENGL_COMPONENT_MOTION_SEQUENCE
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
~OPENGL_COMPONENT_MOTION_SEQUENCE()
{}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
Valid_Frame(int frame_input) const
{
    if(frame_dependent_data) return FILE_UTILITIES::Frame_File_Exists(filename,frame);
    else return motion.time_grid.domain.Lazy_Inside(VECTOR<T,1>(frame_input*one_over_frame_rate));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
Display(const int in_color) const
{
    if(!valid || !draw) return;
    if(slice && slice->Is_Slice_Mode()){glPushAttrib(GL_ENABLE_BIT);slice->Enable_Clip_Planes();}
    opengl_segmented_curve.Display(in_color);
    if(slice && slice->Is_Slice_Mode())glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_segmented_curve.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_MOTION_SEQUENCE<T,RW>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))) return;
    valid=false;
    std::string frame_filename=STRING_UTILITIES::string_sprintf(filename.c_str(),frame);
    if(frame_dependent_data){
        if(FILE_UTILITIES::File_Exists(frame_filename)) FILE_UTILITIES::Read_From_File<RW>(frame_filename,motion);
        else return;}
    if(particles.array_collection->Size() != motion.trajectories.m){particles.array_collection->Delete_All_Elements();particles.array_collection->Add_Elements(motion.trajectories.m);}
    for(int i=1;i<=motion.trajectories.m;i++) particles.X(i)=motion.X(i,one_over_frame_rate*frame);
    frame_loaded=frame;valid=true;
}
//#####################################################################
template class OPENGL_COMPONENT_MOTION_SEQUENCE<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_MOTION_SEQUENCE<double,double>;
#endif
