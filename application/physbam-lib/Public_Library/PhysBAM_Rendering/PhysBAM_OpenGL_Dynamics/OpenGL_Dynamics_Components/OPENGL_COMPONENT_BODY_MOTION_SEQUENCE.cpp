//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Tessellation/CYLINDER_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Dynamics/OpenGL_Dynamics_Components/OPENGL_COMPONENT_BODY_MOTION_SEQUENCE.h>
#include <sstream>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_BODY_MOTION_SEQUENCE
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
OPENGL_COMPONENT_BODY_MOTION_SEQUENCE(const std::string& filename_input,const bool frame_dependent_data_input, const std::string& rigid_body_file,const T scale_input)
    :OPENGL_COMPONENT("Body Motion Sequence"),filename(filename_input),rigid_filename(rigid_body_file),rigid_body_collection(0,0),opengl_component_rigid_body_collection(rigid_body_collection,rigid_body_file),
    scale(scale_input),ui_scale(1),frame_loaded(-1),valid(false),frame_dependent_data(frame_dependent_data_input)
{
    is_animation=true;opengl_component_rigid_body_collection.selectable=true;
    std::string frame_filename=filename;
    if(FILE_UTILITIES::Is_Animated(filename)){frame_filename=STRING_UTILITIES::string_sprintf(filename.c_str(),frame);}
    if(FILE_UTILITIES::File_Exists(frame_filename)) FILE_UTILITIES::Read_From_File<RW>(frame_filename,body_motion);
    default_length=body_motion.trajectories(1)(frame+1).length;
    frame_rate=(body_motion.time_grid.domain.max_corner.x-body_motion.time_grid.domain.min_corner.x)/(body_motion.time_grid.counts.x-1);one_over_frame_rate=(T)1/frame_rate;
}
//#####################################################################
// Function ~OPENGL_COMPONENT_BODY_MOTION_SEQUENCE
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
~OPENGL_COMPONENT_BODY_MOTION_SEQUENCE()
{}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Valid_Frame(int frame_input) const
{
    if(frame_dependent_data) return FILE_UTILITIES::Frame_File_Exists(filename,frame);
    else return body_motion.time_grid.domain.Lazy_Inside(VECTOR<T,1>(frame_input*one_over_frame_rate));
}
//#####################################################################
// Function Prev_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Prev_Frame()
{
    if(frame<=0) return;
    OPENGL_COMPONENT::Set_Frame(frame-1);
    Reinitialize();
}
//#####################################################################
// Function Next_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Next_Frame()
{
    OPENGL_COMPONENT::Set_Frame(frame+1);
    Reinitialize();
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Save_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Save_Frame()
{
    body_motion.saved_frame=frame;
}
//#####################################################################
// Function Get_Frame
//#####################################################################
template<class T,class RW> int OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Get_Frame()
{
    return frame;
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Display(const int in_color) const
{
    if(!valid || !draw) return;
    if(slice && slice->Is_Slice_Mode()){glPushAttrib(GL_ENABLE_BIT);slice->Enable_Clip_Planes();}
    opengl_component_rigid_body_collection.Display(in_color);
    if(slice && slice->Is_Slice_Mode())glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_component_rigid_body_collection.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T,RW>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))) return;
    valid=false;
    if(frame_dependent_data){
        std::string frame_filename=STRING_UTILITIES::string_sprintf(filename.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(frame_filename)) FILE_UTILITIES::Read_From_File<RW>(frame_filename,body_motion);
        else return;}
    if(body_motion.names.m!=opengl_component_rigid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size()){
        opengl_component_rigid_body_collection.rigid_body_collection.rigid_body_particle.Resize(0);
        opengl_component_rigid_body_collection.rigid_body_collection.rigid_geometry_collection.always_create_structure=true;
        for(int i=1;i<=body_motion.names.m;i++){
            //Assume length is invariant and ui is a scale which must be 1 on input
            int id;
            if(rigid_filename=="cyllink"){
                T radius=(T).18,height=(T)(body_motion.trajectories(i)(1).length*scale*.18);int resolution_radius=16,resolution_height=4;
                RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(opengl_component_rigid_body_collection.rigid_body_collection,true);
                CYLINDER<T> cylinder(TV(-height/2,0,0),TV(height/2,0,0),radius);
                rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(cylinder));
                rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(cylinder,resolution_height,resolution_radius));
                id=opengl_component_rigid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);}
            else id=opengl_component_rigid_body_collection.rigid_body_collection.Add_Rigid_Body(STREAM_TYPE(RW()),rigid_filename,scale*body_motion.trajectories(i)(1).length,true,true,false,true);
            body_motion.ui_position(i).length=(T)1;
            opengl_component_rigid_body_collection.rigid_body_collection.rigid_body_particle.X(id)=TV();
            opengl_component_rigid_body_collection.rigid_body_collection.rigid_body_particle.rotation(id)=ROTATION<TV>();
            opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(id).Update_Bounding_Box();
            id_to_index.Insert(id,i);
            opengl_component_rigid_body_collection.Initialize_One_Body(id,true);}}
    for(int id(1);id<=opengl_component_rigid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
        FRAME<TV> rigid_base_transform_i=rigid_body_base_transform;
        T new_length=body_motion.trajectories(id_to_index.Get(id))(frame+1).length;
        rigid_base_transform_i.t*=new_length/default_length;
        opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(id).Set_Frame(body_motion.trajectories(id_to_index.Get(id))(frame+1).targeted_transform*rigid_base_transform_i);
        T scale=body_motion.ui_position(id_to_index.Get(id)).length;
        if(scale!=1&&scale>0){
            scale/=ui_scale;
            std::stringstream ss;
            ss<<"Scaling by "<<scale<<std::endl;
            LOG::filecout(ss.str());
            ui_scale=body_motion.trajectories(id_to_index.Get(id))(frame+1).length/body_motion.base_position(id_to_index.Get(id)).length;
            for(int i=1;i<=opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(id).structures.m;i++){
                if(rigid_filename=="cyllink"){
                    RIGID_BODY<TV>* rigid_body=&opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(id);
                    if(ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >* analytic_object=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >*>(rigid_body->structures(i))){
                        //T radius=.18
                        T height=analytic_object->analytic.height*scale;int resolution_radius=16,resolution_height=4;
                        LOG::filecout("Rescaling cylinder\n");
                        //CYLINDER<T> cylinder(TV(-height/2,0,0),TV(height/2,0,0),radius);
                        analytic_object->analytic.Set_Endpoints(TV(-height/2,0,0),TV(height/2,0,0));
                        //rigid_body->Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(cylinder));
                        rigid_body->Add_Structure(*TESSELLATION::Generate_Triangles(analytic_object->analytic,resolution_height,resolution_radius));}}
                else opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(id).structures(i)->Rescale(scale);}
            opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(id).Rescale(scale);
            body_motion.ui_position(id_to_index.Get(id)).length=1;}}
    opengl_component_rigid_body_collection.Update_Bodies(false);
    frame_loaded=frame;valid=true;
}
//#####################################################################
template class OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<double>;
#endif
