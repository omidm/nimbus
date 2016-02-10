//#####################################################################
// Copyright 2005-2007, Kevin Der, Eftychios Sifakis, Jonathan Su, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MUSCLE_LIST
//#####################################################################
#ifndef __MUSCLE_LIST__
#define __MUSCLE_LIST__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
namespace PhysBAM{

template<class TV> class ATTACHMENT_POINT;
template<class TV> class RIGID_BODY_COLLECTION;

template<class TV>
class MUSCLE_LIST:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    MUSCLE_FORCE_CURVE<T> muscle_force_curve;
    ARRAY<MUSCLE<TV>*> muscles;
    ARRAY<T> muscle_activations;
    ARRAY<ARRAY<TRIPLE<int,ATTACHMENT_POINT<TV>*,ATTACHMENT_POINT<TV>*> >,int> muscle_attachments_on_rigid_body;

    MUSCLE_LIST(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
        :rigid_body_collection(rigid_body_collection_input)
    {}

    void Clean_Memory()
    {muscles.Delete_Pointers_And_Clean_Memory();muscle_activations.Clean_Memory();}

    int Add_Muscle(MUSCLE<TV>* new_muscle)
    {muscles.Append(new_muscle);new_muscle->id=muscles.m;muscle_activations.Append(1);return muscles.m;}

    void Set_Muscle_Activation(const int muscle_id,const T activation)
    {muscle_activations(muscle_id)=activation;}

//#####################################################################
    void Initialize_Muscle_Attachments_On_Rigid_Body();
    void Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame);
    void Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const;
//##################################################################### 
};
}
#endif
