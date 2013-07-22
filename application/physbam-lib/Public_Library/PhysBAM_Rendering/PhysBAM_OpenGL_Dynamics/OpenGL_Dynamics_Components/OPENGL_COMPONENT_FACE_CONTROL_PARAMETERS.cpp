//#####################################################################
// Copyright 2005-2007, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Dynamics/OpenGL_Dynamics_Components/OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS.h>
#include <PhysBAM_Dynamics/Motion/ACTIVATION_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/ATTACHMENT_FRAME_CONTROL_SET.h>
#include <sstream>
#include <string>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<T,RW>::
OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS(const std::string& basedir_input,const std::string& filename_input)
    :basedir(basedir_input),filename(filename_input),valid(false),frame_loaded(-1)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(basedir+"/face_control_set_types");
    if(input){
        int control_sets;
        Read_Binary<RW>(*input,control_sets);
        for(int i=1;i<=control_sets;i++){
            int control_set_type;
            Read_Binary<RW>(*input,control_set_type);
            if(control_set_type==FACE_CONTROL_SET<T>::ACTIVATION){
                ACTIVATION_CONTROL_SET<T>* control_set=new ACTIVATION_CONTROL_SET<T>;
                ARRAY<std::string> muscle_names;
                FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/face_control_set_%d.muscle_names",basedir.c_str(),i),muscle_names);
                std::stringstream ss;
                for(int j=1;j<=muscle_names.m;j++){ss<<"Muscle Name "<<j<<" is "<<muscle_names(j)<<std::endl;control_set->Add_Activation(muscle_names(j));}
                LOG::filecout(ss.str());
                face_control_parameters.list.Append(control_set);}
            else if(control_set_type==FACE_CONTROL_SET<T>::ATTACHMENT_FRAME){
                ARRAY<VECTOR<T,3> > X_dummy;ARRAY<ARRAY<int> > nodes_dummy;
                ATTACHMENT_FRAME_CONTROL_SET<T>* control_set=new ATTACHMENT_FRAME_CONTROL_SET<T>(X_dummy,nodes_dummy,0);
                control_set->Read_Jaw_Joint_From_File(STREAM_TYPE(RW()),STRING_UTILITIES::string_sprintf("%s/face_control_set_%d.jaw_joint_parameters",basedir.c_str(),i));
                face_control_parameters.list.Append(control_set);}}
        delete input;}
    is_animation=FILE_UTILITIES::Is_Animated(filename);
}
//#####################################################################
// Function ~OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<T,RW>::
~OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS()
{
    for(int i=1;i<=face_control_parameters.list.m;i++) delete face_control_parameters.list(i);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<T,RW>::
Display(const int in_color) const
{
    if(!valid || !draw) return;

    glMatrixMode(GL_MODELVIEW);glPushMatrix();glLoadIdentity();
    glMatrixMode(GL_PROJECTION);glPushMatrix();glLoadIdentity();
    glPushAttrib(GL_ENABLE_BIT);glDisable(GL_DEPTH_TEST);glDisable(GL_LIGHTING);
    int gl_width=OPENGL_WORLD::Singleton()->window->Width();
    int gl_height=OPENGL_WORLD::Singleton()->window->Height();
    int horizontal_offset=gl_width-540,vertical_offset=gl_height-30;
    gluOrtho2D(0,gl_width,0,gl_height);
    // get strings from controls and convert to list array of strings
    ARRAY<std::string> strings;
    std::ostringstream output_stream;
    face_control_parameters.Print_Diagnostics(output_stream);std::string text=output_stream.str();
    for(std::string::size_type start=0;start<text.length();){
        std::string::size_type end=text.find('\n',start);
        strings.Append(text.substr(start,end-start));
        if(end==std::string::npos) break;start=end+1;}
    int vspace=13;
    // make box
    OPENGL_WORLD::Draw_Transparent_Text_Box(strings,VECTOR<int,2>(horizontal_offset,vertical_offset),vspace,GLUT_BITMAP_8_BY_13,OPENGL_COLOR::Gray(0,0.5));
    // print strings
    OPENGL_COLOR::White().Send_To_GL_Pipeline();
    T height=(T)vertical_offset;
    for(int i=1; i<= strings.m; i++){OpenGL_String(VECTOR<float,2>((float)horizontal_offset,(float)height),strings(i),GLUT_BITMAP_8_BY_13);height-=vspace;}
    // reet opengl state
    glPopAttrib();
    glPopMatrix();glMatrixMode(GL_MODELVIEW);glPopMatrix();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<T,RW>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))) return;
    
    std::string frame_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(frame_filename);
    frame_loaded=frame;valid=true;
    face_control_parameters.template Read<RW>(*input);
    delete input;
}
template class OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_FACE_CONTROL_PARAMETERS<double,double>;
#endif
