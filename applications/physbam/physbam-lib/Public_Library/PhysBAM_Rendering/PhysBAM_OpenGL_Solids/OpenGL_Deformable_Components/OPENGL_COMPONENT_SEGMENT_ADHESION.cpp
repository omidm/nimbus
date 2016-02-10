//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_MOTION_SEQUENCE
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_ADHESION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Deformable_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Deformable_Components/OPENGL_COMPONENT_SEGMENT_ADHESION.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_SEGMENT_ADHESION
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
OPENGL_COMPONENT_SEGMENT_ADHESION(const std::string& filename,const OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>& opengl_component_deformable_body,const OPENGL_COLOR& color)
    :OPENGL_COMPONENT("Segment Adhesion"),filename(filename),opengl_component_deformable_body(opengl_component_deformable_body),segmented_curve(segment_mesh,particles),
    opengl_segmented_curve(segmented_curve),frame_loaded(-1),valid(false)
{
    opengl_segmented_curve.color=color;
    is_animation=true;opengl_segmented_curve.draw_vertices=true;
}
//#####################################################################
// Function ~OPENGL_COMPONENT_SEGMENT_ADHESION
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
~OPENGL_COMPONENT_SEGMENT_ADHESION()
{}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
Valid_Frame(const int frame) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
Display(const int in_color) const
{
    opengl_segmented_curve.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_segmented_curve.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SEGMENT_ADHESION<T,RW>::
Reinitialize(bool force)
{
    const PARTICLES<TV>& real_particles=opengl_component_deformable_body.deformable_body_collection.particles;
    (const_cast<OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T,RW>& >(opengl_component_deformable_body)).Reinitialize(force);

    if(!draw || !(force || !valid || (is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0))) return;
    valid=false;
    std::string frame_filename=STRING_UTILITIES::string_sprintf(filename.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(frame_filename)){
        HASHTABLE<VECTOR<int,2>,typename SEGMENT_ADHESION<TV>::SPRING_STATE> springs;
        FILE_UTILITIES::Read_From_File<RW>(frame_filename,springs);
        segment_mesh.elements.Remove_All();
        particles.array_collection->Delete_All_Elements();
        int pairs=0;
        for(HASHTABLE_ITERATOR<VECTOR<int,2>,const typename SEGMENT_ADHESION<TV>::SPRING_STATE> iterator(springs);iterator.Valid();iterator.Next()){
            const typename SEGMENT_ADHESION<TV>::SPRING_STATE& state=iterator.Data();
            int i=particles.array_collection->Add_Element(),j=particles.array_collection->Add_Element();
            particles.X(i)=(1-state.weights[1])*real_particles.X(state.nodes[1])+state.weights[1]*real_particles.X(state.nodes[2]);
            particles.X(j)=(1-state.weights[2])*real_particles.X(state.nodes[3])+state.weights[2]*real_particles.X(state.nodes[4]);
            segment_mesh.elements.Append(VECTOR<int,2>(i,j));
            pairs++;
        }
        std::stringstream ss;
        ss<<"pairs is "<<pairs<<" and number of particles is "<<particles.array_collection->Size()<<std::endl;
        LOG::filecout(ss.str());
    }
}
//#####################################################################
template class OPENGL_COMPONENT_SEGMENT_ADHESION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_SEGMENT_ADHESION<double>;
#endif
