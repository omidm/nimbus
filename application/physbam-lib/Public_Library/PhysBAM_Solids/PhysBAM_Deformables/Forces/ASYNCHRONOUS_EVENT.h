//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ASYNCHRONOUS_EVENT
//##################################################################### 
#ifndef __ASYNCHRONOUS_EVENT__
#define __ASYNCHRONOUS_EVENT__

namespace PhysBAM{

template<class T>
class ASYNCHRONOUS_EVENT
{
public:
    virtual ~ASYNCHRONOUS_EVENT(){}
    virtual T Initialize(const T time)=0; // returns next_time
    virtual T Process(const T next_time)=0; // returns next_time
    virtual T Recompute_Next_Time(const T next_time){return next_time;}
//#####################################################################
};
}
#endif
