//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTIFICIAL_VISCOSITY 
//##################################################################### 
#ifndef __ARTIFICIAL_VISCOSITY__
#define __ARTIFICIAL_VISCOSITY__    

class ARTIFICIAL_VISCOSITY
{
protected:
    int limiter;

    ARTIFICIAL_VISCOSITY()
    {
        Use_Limiter();
    }

public:
    void Use_Limiter()
    {limiter=1;}
    
    void Do_Not_Use_Limiter()
    {limiter=0;}

//#####################################################################
};   
#endif

