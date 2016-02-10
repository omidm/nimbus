//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h> // TODO: remove once MUSCLE.cpp exists (workaround for windows compiler bug)
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINT_SIMPLICES_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T,RW>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(const std::string& basedir_input)
    :OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>(basedir_input,false),rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0,0)),
    need_destroy_rigid_body_collection(true)
{
    rigid_geometry_collection=&rigid_body_collection.rigid_geometry_collection;
    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T,RW>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir_input)
    :OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>(basedir_input,false),
    rigid_body_collection(rigid_body_collection),
    need_destroy_rigid_body_collection(false)
{
    rigid_geometry_collection=&rigid_body_collection.rigid_geometry_collection;
    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T,RW>::
~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D()
{
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T,RW>::
Reinitialize(const bool force,const bool read_geometry)
{   
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame))) return;
        rigid_body_collection.Read(STREAM_TYPE(RW()),basedir,frame,&needs_init,&needs_destroy);
        BASE::Reinitialize(force,false);}
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T,RW>::
Update_Object_Labels()
{
    BASE::Update_Object_Labels();
    for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        if(draw_object(i)){
            if(opengl_point_simplices(i)){
                if(output_positions){
                    rigid_body_collection.Rigid_Body(i).Update_Angular_Velocity();}}}}
}
//#####################################################################
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<double,double>;
#endif
