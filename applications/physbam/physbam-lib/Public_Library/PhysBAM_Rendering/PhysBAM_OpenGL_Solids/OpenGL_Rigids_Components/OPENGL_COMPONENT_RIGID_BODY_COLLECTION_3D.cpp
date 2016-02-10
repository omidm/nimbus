//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/DUALCONTOUR_OCTREE.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/NORMAL_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_MUSCLE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D.h>
#include <sstream>
#include <string>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(bool use_display_lists)
    :OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>(use_display_lists,false),
    rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0,0)),rigid_body_collection_simulated(0),articulated_rigid_body(0),current_selection(0),need_destroy_rigid_body_collection(true)
{
    rigid_geometry_collection=&rigid_body_collection.rigid_geometry_collection;
    Initialize();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(const std::string& basedir_input,bool use_display_lists)
    :OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>(basedir_input,use_display_lists,false),
    rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0,0)),rigid_body_collection_simulated(0),articulated_rigid_body(0),current_selection(0),need_destroy_rigid_body_collection(true)
{
    rigid_geometry_collection=&rigid_body_collection.rigid_geometry_collection;
    Initialize();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir_input,bool use_display_lists)
    :OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>(basedir_input,use_display_lists,false),
    rigid_body_collection(rigid_body_collection),rigid_body_collection_simulated(0),articulated_rigid_body(0),current_selection(0),need_destroy_rigid_body_collection(false)
{
    rigid_geometry_collection=&rigid_body_collection.rigid_geometry_collection;
    Initialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D()
{
    delete opengl_component_muscle_3d;
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Initialize()
{
    draw_articulation_points=false;
    draw_joint_frames=0;
    draw_forces_and_torques=false;
    has_init_destroy_information=true;

    OPENGL_COLOR_RAMP<T>* color_ramp=new OPENGL_COLOR_RAMP<T>();
    color_ramp->Add_Color((T)-log(1.001),OPENGL_COLOR::Blue());
//    color_ramp->Add_Color((T)-log(1.000001),OPENGL_COLOR::Blue());
    color_ramp->Add_Color((T)0,OPENGL_COLOR::Blue((T).5),OPENGL_COLOR::Gray((T).5),OPENGL_COLOR::Red((T).5));
//    color_ramp->Add_Color(log((T)1.00001),OPENGL_COLOR::Red());
    color_ramp->Add_Color(log((T)1.001),OPENGL_COLOR::Red());

    opengl_component_muscle_3d=new OPENGL_COMPONENT_MUSCLE_3D<T>(color_ramp);
}
//#####################################################################
// Function Read_Articulated_Information
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Read_Articulated_Information(const std::string& filename)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    // this will need to be changed to reflect multiple articulation points per rigid body
    int numpoints=0;Read_Binary<RW>(*input,numpoints);articulation_points.Exact_Resize(numpoints);joint_frames.Exact_Resize(numpoints);
    for(int i=1;i<=numpoints;i++) Read_Binary<RW>(*input,articulation_points(i),joint_frames(i));
    Read_Binary<RW>(*input,projected_COM);delete input;
}
//#####################################################################
// Function Update_Articulation_Points
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Update_Articulation_Points()
{
    int num_points=2*articulated_rigid_body->joint_mesh.joints.m;
    articulation_points.Exact_Resize(num_points);
    joint_frames.Exact_Resize(num_points);
    for(int i=1;i<=num_points;i+=2){
        int index=(i+1)/2;
        JOINT<TV>* joint=articulated_rigid_body->joint_mesh.joints(index);
        RIGID_BODY<TV>* parent=articulated_rigid_body->Parent(joint->id_number),*child=articulated_rigid_body->Child(joint->id_number);
        articulation_points(i)=parent->World_Space_Point(joint->F_pj().t);articulation_points(i+1)=child->World_Space_Point(joint->F_cj().t);
        joint_frames(i)=parent->Frame()*joint->F_pj();joint_frames(i+1)=child->Frame()*joint->F_cj();}
}
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Resize_Structures(const int size)
{
    extra_components.Resize(size);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(is_interactive){
        if(draw && (force || frame_loaded!=frame)){
            valid = false;
            int max_number_of_bodies=max(extra_components.Size(),rigid_body_collection.rigid_body_particle.array_collection->Size());
            Resize_Structures(max_number_of_bodies);}
        BASE::Reinitialize_From_Simulation(force);}
    else if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame))) return;

        // TODO: currently reads in all structures, should only read in certain kinds based on read_triangulated_surface,read_implicit_surface,read_tetrahedralized_volume
        rigid_body_collection.Read(STREAM_TYPE(RW()),basedir,frame,&needs_init,&needs_destroy);

        std::string arb_state_file=STRING_UTILITIES::string_sprintf("%s/%d/arb_state",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(arb_state_file)){
            if(!articulated_rigid_body) articulated_rigid_body=new ARTICULATED_RIGID_BODY<TV>(rigid_body_collection); // TODO: read in the actual particles
            articulated_rigid_body->Read(STREAM_TYPE(RW()),basedir,frame);
            Initialize();
            if(opengl_component_muscle_3d) opengl_component_muscle_3d->Initialize(articulated_rigid_body,basedir,frame);}
        else{delete articulated_rigid_body;articulated_rigid_body=0;delete opengl_component_muscle_3d;opengl_component_muscle_3d=0;}

        if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame))){
            Read_Articulated_Information(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame));}

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(extra_components.Size(),rigid_body_collection.rigid_body_particle.array_collection->Size());
        Resize_Structures(max_number_of_bodies);

        std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_forces_and_torques",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename,forces_and_torques);
        else forces_and_torques.Resize(0);

        BASE::Reinitialize(force,false);}
}
//#####################################################################
// Function Reinitialize_Without_Files
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Reinitialize_Without_Files(const bool force)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        LOG::filecout("Reinit being called\n");
        valid=false;

        ARRAY<int> needs_init;
        for(int i(1);i<=articulated_rigid_body->rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(articulated_rigid_body->rigid_body_collection.Is_Active(i)) needs_init.Append(i);

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(extra_components.Size(),rigid_body_collection.rigid_body_particle.array_collection->Size());
        Resize_Structures(max_number_of_bodies);
        BASE::Resize_Structures(max_number_of_bodies);

        // Initialize bodies which have become active
        for(int i=1;i<=needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_body_collection.Is_Active(id));
            Create_Geometry(id);}
        
        BASE::Reinitialize_Without_Files(force);}
}
//#####################################################################
// Function Initialize_One_Body
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Initialize_One_Body(const int body_id,const bool force)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        LOG::filecout("Init for one body being called\n");
        valid=false;

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(extra_components.Size(),rigid_body_collection.rigid_body_particle.array_collection->Size());
        Resize_Structures(max_number_of_bodies);

        BASE::Initialize_One_Body(body_id,force);}
}
//#####################################################################
// Function Update_Bodies
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Update_Bodies(const bool update_arb_points)
{
    if(articulated_rigid_body && update_arb_points) Update_Articulation_Points();
    BASE::Update_Bodies(update_arb_points);
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Create_Geometry(const int id)
{
    BASE::Create_Geometry(id);

    // add extra components
    if(opengl_tetrahedralized_volume(id)){
        RIGID_GEOMETRY<TV>& rigid_geometry=rigid_body_collection.rigid_geometry_collection.Rigid_Geometry(id);
        TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=rigid_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>();
        std::string filename_pattern=STRING_UTILITIES::string_sprintf("%s/accumulated_impulses_%d.%%d",basedir.c_str(),id);
        if(FILE_UTILITIES::Frame_File_Exists(filename_pattern,frame)){
            std::stringstream ss;
            ss<<"Adding accumulated impulses to rigid body "<<id<<std::endl;
            LOG::filecout(ss.str());
            OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T,RW>* component=
                new OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T,RW>(*tetrahedralized_volume,filename_pattern);
            component->opengl_vector_field.Enslave_Transform_To(*opengl_axes(id));
            extra_components(id).Append(component);}}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Update_Geometry(const int id)
{
    BASE::Update_Geometry(id);
    for(int i=1;i<=extra_components(id).m;i++) extra_components(id)(i)->Set_Frame(frame);
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Destroy_Geometry(const int id)
{
    if(draw_object.Size()<id)return; // it's possible that we try to delete geometry that we've never added because of substepping in the simulator
    BASE::Destroy_Geometry(id);
    extra_components(id).Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Display(const int in_color) const
{
    BASE::Display(in_color);

    if(!draw) return;
    GLint mode;glGetIntegerv(GL_RENDER_MODE,&mode);
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}

    // Articulated rigid bodies
    if(articulated_rigid_body){
        if(draw_articulation_points){
            OPENGL_COLOR articulation_point_color(0.8f,0.8f,0.2f),segment_color(0,0,1),com_color(1,0,0);
            if(mode==GL_SELECT){glPushName(4);glPushAttrib(GL_POINT_BIT);glPointSize(OPENGL_PREFERENCES::selection_point_size);}
            for(int i=1;i<=articulation_points.m;i++){
                glPushName(i);
                OPENGL_SHAPES::Draw_Dot(articulation_points(i),articulation_point_color,5);
                glPopName();}
            OPENGL_SHAPES::Draw_Dot(projected_COM,com_color,10);
            if(mode!=GL_SELECT) for(int i=1;i<=articulation_points.m;i+=2){
                OPENGL_SHAPES::Draw_Segment(articulation_points(i),articulation_points(i+1),segment_color,5);}
            if(mode==GL_SELECT){glPopName();glPopAttrib();}
            if(mode!=GL_SELECT && current_selection && current_selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_3D){
                int joint_id=((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>*)current_selection)->joint_id;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(articulation_points(joint_id));}}
        if(opengl_component_muscle_3d) opengl_component_muscle_3d->Display(in_color);}

    RANGE<VECTOR<T,3> > axes_box(0,2,0,2,0,2);
    //RANGE<VECTOR<T,3> > axes_box(0,velocity_field.size,0,velocity_field.size,0,velocity_field.size);
    if(draw_joint_frames==1) for(int i=1;i<=joint_frames.m;i++)(OPENGL_AXES<T>(joint_frames(i),axes_box)).Display();
    else if(draw_joint_frames==2) for(int i=2;i<=joint_frames.m;i+=2)(OPENGL_AXES<T>(joint_frames(i),axes_box)).Display();
    else if(draw_joint_frames==3) for(int i=1;i<=joint_frames.m;i+=2)(OPENGL_AXES<T>(joint_frames(i),axes_box)).Display();

    if(draw_forces_and_torques && forces_and_torques.Size()==rigid_body_collection.rigid_body_particle.array_collection->Size()){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        T scale=(T)velocity_field.size/24;
        OPENGL_COLOR::Yellow().Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int i(1);i<=forces_and_torques.Size();i++) if(rigid_body_collection.Is_Active(i)){
            OpenGL_Line(rigid_body_collection.rigid_body_particle.X(i),rigid_body_collection.rigid_body_particle.X(i)+scale*forces_and_torques(i).x);}
        OpenGL_End();
        for(int i(1);i<=forces_and_torques.Size();i++) if(rigid_body_collection.Is_Active(i)){
            std::string label=STRING_UTILITIES::string_sprintf("F=%.3f %.3f %.3f, T=%.3f %.3f %.3f",forces_and_torques(i).x.x,forces_and_torques(i).x.y,forces_and_torques(i).x.z,forces_and_torques(i).y.x,forces_and_torques(i).y.y,forces_and_torques(i).y.z);
            OpenGL_String(rigid_body_collection.rigid_body_particle.X(i)+scale*forces_and_torques(i).x,label);}
        glPopAttrib();}

    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size>=2){
        if(buffer[0]<=3) return BASE::Get_Selection(buffer,buffer_size);
        else if(buffer[0]==4){ // articulation joints
            selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>(this,buffer[1]);}
        else if(opengl_component_muscle_3d) selection=opengl_component_muscle_3d->Get_Muscle_Selection(buffer,buffer_size);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    BASE::Highlight_Selection(selection);

    delete current_selection;current_selection=0;
    if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_3D){
        current_selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>(this,((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>*)selection)->joint_id);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Clear_Highlight()
{
    BASE::Clear_Highlight();

    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    if(!selection || selection->object!=this) return;

    if(selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>*)selection;

        output_stream<<"Rigid body "<<real_selection->body_id<<std::endl;

        // handle case of body no longer being active
        if(!rigid_body_collection.Is_Active(real_selection->body_id) || !opengl_triangulated_surface(real_selection->body_id)){
            output_stream<<"INACTIVE"<<std::endl;
            return;}

        RIGID_BODY<TV> *body=const_cast<RIGID_BODY<TV>*>(&rigid_body_collection.Rigid_Body(real_selection->body_id));
        body->Update_Angular_Velocity();

        if(!body->name.empty()){output_stream<<"Name ="<<body->name<<std::endl;}
        output_stream<<"Mass = "<<body->Mass()<<std::endl;
        output_stream<<"Inertia tensor = "<<body->Inertia_Tensor()<<std::endl;
        output_stream<<"Angular momentum = "<<body->Angular_Momentum()<<std::endl;
        output_stream<<"Kinetic energy = "<<body->Kinetic_Energy()<<std::endl;
        output_stream<<std::endl;
        BASE::Print_Selection_Info(output_stream,selection);}
    else if(selection->type==OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_3D){
        OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T> *real_selection=(OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>*)selection;
        int articulation_id=real_selection->joint_id;int joint_index=(articulation_id+1/2);
        JOINT<TV>& joint=*articulated_rigid_body->joint_mesh.joints(joint_index);JOINT_ID joint_id=joint.id_number;
        const RIGID_BODY<TV> *parent=articulated_rigid_body->Parent(joint_id),*child=articulated_rigid_body->Child(joint_id);
        output_stream<<"Joint id "<<joint_id<<" ("<<(!joint.name.empty()?joint.name:"UNNAMED")<<")"<<std::endl;
        const char* joint_type_string;
        const std::type_info& type=typeid(joint);
        if(type==typeid(POINT_JOINT<TV>)) joint_type_string="POINT_JOINT";
        else if(type==typeid(RIGID_JOINT<TV>)) joint_type_string="RIGID_JOINT";
        else if(type==typeid(ANGLE_JOINT<TV>)) joint_type_string="ANGLE_JOINT";
        else if(type==typeid(PRISMATIC_TWIST_JOINT<TV>)) joint_type_string="PRISMATIC_TWIST_JOINT";
        else if(type==typeid(NORMAL_JOINT<TV>)) joint_type_string="NORMAL_JOINT";
        else joint_type_string="unknown";
        output_stream<<"Type = "<<joint_type_string<<std::endl;
        output_stream<<"Frame: "<<(articulation_id%2==0?"child":"parent")<<std::endl;
        output_stream<<"Parent = "<<parent->name<<std::endl;
        output_stream<<"Child = "<<child->name<<std::endl;
        output_stream<<"Articulation point "<<articulation_points(articulation_id)<<std::endl;

        FRAME<TV> parent_frame=joint_frames(2*Value(joint_index)-1),child_frame=joint_frames(2*Value(joint_index));
        FRAME<TV> joint_frame=parent_frame.Inverse()*child_frame;
        output_stream<<"Joint translation: "<<joint_frame.t<<std::endl;
        output_stream<<"Joint rotation vector: "<<joint_frame.r.Rotation_Vector()<<std::endl;
        T twist,phi,theta;joint_frame.r.Euler_Angles(twist,phi,theta);
        output_stream<<"Joint Euler angles: twist="<<twist<<", phi="<<phi<<", theta="<<theta<<std::endl;
        VECTOR<T,3> ap1=parent->World_Space_Point(joint.F_pj().t),ap2=child->World_Space_Point(joint.F_cj().t),location=(T).5*(ap1+ap2);

        VECTOR<T,3> current_relative_velocity=-RIGID_BODY<TV>::Relative_Velocity(*parent,*child,location); // child w.r.t. parent!
        VECTOR<T,3> current_relative_angular_velocity=-RIGID_BODY<TV>::Relative_Angular_Velocity(*parent,*child); // child w.r.t. parent!
        output_stream<<"Relative velocity at joint = "<<current_relative_velocity<<std::endl;
        output_stream<<"Relative angular velocity = "<<current_relative_angular_velocity<<std::endl;}
    else if(selection->type==OPENGL_SELECTION::MUSCLE_3D){
        if(opengl_component_muscle_3d) opengl_component_muscle_3d->Print_Selection_Info(output_stream,selection);}
}
//#####################################################################
// Function Toggle_Articulation_Points
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Toggle_Articulation_Points()
{
    draw_articulation_points=!draw_articulation_points;
}
//#####################################################################
// Function Toggle_Joint_Frames
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Toggle_Joint_Frames()
{
    draw_joint_frames++;
    draw_joint_frames%=4;
}
//#####################################################################
// Function Toggle_Forces_And_Torques
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T,RW>::
Toggle_Forces_And_Torques()
{
    draw_forces_and_torques=!draw_forces_and_torques;
}
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<double,double>;
#endif
