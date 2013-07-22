//#####################################################################
// Copyright 2006-2007, Kevin Der, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MUSCLE_SEGMENT
//#####################################################################
#ifndef __MUSCLE_SEGMENT__
#define __MUSCLE_SEGMENT__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <string>
namespace PhysBAM{

template<class TV> class ATTACHMENT_POINT;
template<class TV> class GEOMETRY_PARTICLES;
template<class T> class QUEUE;

template<class TV>
class MUSCLE_SEGMENT:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    static const int activation_memory=5;
    typedef TV VECTOR_T;

    enum MUSCLE_SEGMENT_TYPE {LINEAR_SEGMENT=1,ANALYTIC_SURFACE_SEGMENT,NUM_MUSCLE_SEGMENT_TYPE};
    MUSCLE_SEGMENT_TYPE segment_type;
    ATTACHMENT_POINT<TV> *point_1,*point_2;
    FRAME<TV> frame;
    QUEUE<T>* activations; // don't need to save since it only affects the parameter computation, which is saved

    MUSCLE_SEGMENT();
    MUSCLE_SEGMENT(ATTACHMENT_POINT<TV>* point_1_input,ATTACHMENT_POINT<TV>* point_2_input);

    virtual ~MUSCLE_SEGMENT();

    void Set_Frame(const FRAME<TV>& frame_input)
    {frame=frame_input;}

//#####################################################################
    T Length() const;
    virtual void Set_Current_Activation(const T activation);
    virtual void Update_Parameters();
    virtual void Initialize();
    virtual void Read_And_Set_Parameters(TYPED_ISTREAM& input);
    virtual void Write_Parameters(TYPED_OSTREAM& output) const;
    virtual std::string Name() const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return Static_Name();}
    virtual std::string Extension() const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return Static_Extension();}
    static std::string Static_Name();
    static std::string Static_Extension();
    virtual T Maximum_Radius() const;
    virtual TV World_Space_Position(const TV& normalized_local_space_position);
    virtual TV Get_Current_World_Space_Position(const TV& normalized_local_space_point);
    static MUSCLE_SEGMENT* Create();
    void Write(TYPED_OSTREAM& output) const;
    void Update_Frame();
    static MUSCLE_SEGMENT* Create_From_Input(TYPED_ISTREAM& input);
//#####################################################################
};
}
#endif
