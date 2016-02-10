//#####################################################################
// Copyright 2006-2007, Kevin Der, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Geometry/Registry/REGISTRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ANALYTIC_SURFACE_MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_SEGMENT.h>
using namespace PhysBAM;
//#####################################################################
template<class TV> const int MUSCLE_SEGMENT<TV>::activation_memory;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MUSCLE_SEGMENT<TV>::
MUSCLE_SEGMENT()
    :activations(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MUSCLE_SEGMENT<TV>::
MUSCLE_SEGMENT(ATTACHMENT_POINT<TV>* point_1_input,ATTACHMENT_POINT<TV>* point_2_input)
    :segment_type(LINEAR_SEGMENT),point_1(point_1_input),point_2(point_2_input),activations(new QUEUE<T>(activation_memory))
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MUSCLE_SEGMENT<TV>::
~MUSCLE_SEGMENT()
{
    delete activations;
}
//#####################################################################
// Function Update_Frame
//#####################################################################
template<class TV> void MUSCLE_SEGMENT<TV>::
Update_Frame()
{
    TV point_1_embedded_position=point_1->Embedded_Position();
    // TODO: use y direction for non-axisymmetric surfaces
    ROTATION<TV> rotation=ROTATION<TV>::From_Rotated_Vector(TV::Axis_Vector(1),point_2->Embedded_Position()-point_1_embedded_position);
    Set_Frame(FRAME<TV>(point_1_embedded_position,rotation));
}
//#####################################################################
// Class MUSCLE_SEGMENT_REGISTRY
//#####################################################################
namespace {

template<class TV> class MUSCLE_SEGMENT_REGISTRY;

template<class TV>
class MUSCLE_SEGMENT_REGISTRY:public REGISTRY<MUSCLE_SEGMENT<TV>,std::string,MUSCLE_SEGMENT_REGISTRY<TV> >
{};

void Register_Muscle_Segment()
{
    static bool done=false;if(done) return;done=true;
    MUSCLE_SEGMENT_REGISTRY<VECTOR<float,2> >::Register<MUSCLE_SEGMENT<VECTOR<float,2> > >();
    MUSCLE_SEGMENT_REGISTRY<VECTOR<float,3> >::Register<MUSCLE_SEGMENT<VECTOR<float,3> > >();
    MUSCLE_SEGMENT_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_SURFACE_MUSCLE_SEGMENT<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    MUSCLE_SEGMENT_REGISTRY<VECTOR<double,2> >::Register<MUSCLE_SEGMENT<VECTOR<double,2> > >();
    MUSCLE_SEGMENT_REGISTRY<VECTOR<double,3> >::Register<MUSCLE_SEGMENT<VECTOR<double,3> > >();
    MUSCLE_SEGMENT_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_SURFACE_MUSCLE_SEGMENT<double> >();
#endif
}
}
//#####################################################################
// Function Create_From_Input
//#####################################################################
template<class TV> MUSCLE_SEGMENT<TV>* MUSCLE_SEGMENT<TV>::
Create_From_Input(TYPED_ISTREAM& input)
{
    Register_Muscle_Segment();
    std::string name;Read_Binary(input,name);
    MUSCLE_SEGMENT* muscle_segment=MUSCLE_SEGMENT_REGISTRY<TV>::Name_To_Factory(name)->Create();
    muscle_segment->Read_And_Set_Parameters(input);return muscle_segment;
}
//#####################################################################
// Function Read_And_Set_Parameters
//#####################################################################
template<class TV> void MUSCLE_SEGMENT<TV>::
Read_And_Set_Parameters(TYPED_ISTREAM& input)
{
    Read_Binary(input,frame,segment_type);
}
//#####################################################################
// Function Write_Parameters
//#####################################################################
template<class TV> void MUSCLE_SEGMENT<TV>::
Write_Parameters(TYPED_OSTREAM& output) const
{
    Write_Binary(output,frame,segment_type);
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string MUSCLE_SEGMENT<TV>::
Static_Name()
{
    return "muscle_segment";
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV> std::string MUSCLE_SEGMENT<TV>::
Static_Extension()
{
    return "";
}
//#####################################################################
// Function Maximum_Radius
//#####################################################################
template<class TV> typename TV::SCALAR MUSCLE_SEGMENT<TV>::
Maximum_Radius() const
{
    return 0;
}
//#####################################################################
// Function World_Space_Position
//#####################################################################
template<class TV> TV MUSCLE_SEGMENT<TV>::
World_Space_Position(const TV& normalized_local_space_position)
{
    assert(normalized_local_space_position.x>=0 && normalized_local_space_position.x<=1);
    TV local_space_position(normalized_local_space_position);local_space_position.x*=Length();
    return frame*local_space_position;
}
//#####################################################################
// Function Get_Current_World_Space_Position
//#####################################################################
template<class TV> TV MUSCLE_SEGMENT<TV>::
Get_Current_World_Space_Position(const TV& normalized_local_space_point)
{
    return TV();
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> MUSCLE_SEGMENT<TV>* MUSCLE_SEGMENT<TV>::
Create()
{
    return new MUSCLE_SEGMENT();
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void MUSCLE_SEGMENT<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,Name());
    Write_Parameters(output);
}
//#####################################################################
// Function Update_Parameters
//#####################################################################
template<class TV> void MUSCLE_SEGMENT<TV>::
Update_Parameters()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MUSCLE_SEGMENT<TV>::
Initialize()
{
    Update_Frame();
}
//#####################################################################
// Function Set_Current_Activation
//#####################################################################
template<class TV> void MUSCLE_SEGMENT<TV>::
Set_Current_Activation(const T activation)
{
    if(activations->Full()) activations->Dequeue();
    activations->Enqueue(activation);
}
//#####################################################################
// Function Length
//#####################################################################
template<class TV> typename TV::SCALAR MUSCLE_SEGMENT<TV>::
Length() const
{
    return (point_1->Embedded_Position()-point_2->Embedded_Position()).Magnitude();
}
//#####################################################################
template class MUSCLE_SEGMENT<VECTOR<float,1> >;
template class MUSCLE_SEGMENT<VECTOR<float,2> >;
template class MUSCLE_SEGMENT<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MUSCLE_SEGMENT<VECTOR<double,1> >;
template class MUSCLE_SEGMENT<VECTOR<double,2> >;
template class MUSCLE_SEGMENT<VECTOR<double,3> >;
#endif
