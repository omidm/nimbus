//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_SKIP_COLLISION_CHECK
//#####################################################################
#ifndef __RIGID_BODY_SKIP_COLLISION_CHECK__
#define __RIGID_BODY_SKIP_COLLISION_CHECK__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class RIGID_BODY_SKIP_COLLISION_CHECK
{
private:
    unsigned int skip_counter;
    //HASHTABLE<VECTOR<int,2>,unsigned int> pair_last_checked; // occur only at even times
    ARRAY<HASHTABLE<int,unsigned int> > pair_last_checked; // occur only at even times
    ARRAY<unsigned int,int> rigid_body_last_moved; // occur only at odd times

public:
    RIGID_BODY_SKIP_COLLISION_CHECK(const int number_of_rigid_bodies=0)
    {
        Initialize(number_of_rigid_bodies,true);
    }

    void Reset()
    {skip_counter=1;
    for(int i=1;i<=pair_last_checked.m;i++) pair_last_checked(i).Remove_All();
    ARRAYS_COMPUTATIONS::Fill(rigid_body_last_moved,(unsigned int)1);}

    void Initialize(const int number_of_rigid_bodies,const bool reset)
    {if(reset){
        pair_last_checked.Resize(number_of_rigid_bodies,false,false);
        rigid_body_last_moved.Resize(number_of_rigid_bodies,false,false);
        Reset();}
    else{
        for(int i=1;i<=pair_last_checked.m;i++) pair_last_checked(i).Remove_All();
        Resize(number_of_rigid_bodies);}}

    void Resize(const int number_of_rigid_bodies)
    {pair_last_checked.Resize(number_of_rigid_bodies,false,true);
    int old_size=rigid_body_last_moved.m;
    rigid_body_last_moved.Resize(number_of_rigid_bodies,false,true);
    for(int i=old_size+1;i<=number_of_rigid_bodies;i++) rigid_body_last_moved(i)=1;}

    void Set_Last_Checked(const int id_1,const int id_2)
    {if(skip_counter&1) skip_counter++;
    if(id_1>pair_last_checked.m) pair_last_checked.Resize(Value(id_1));
    pair_last_checked(Value(id_1)).Set(Value(id_2),skip_counter);
    if(id_2>pair_last_checked.m) pair_last_checked.Resize(Value(id_2));
    pair_last_checked(Value(id_2)).Set(Value(id_1),skip_counter);}
    //pair_last_checked.Set(VECTOR<int,2>(Value(id_1),Value(id_2)),skip_counter);pair_last_checked.Set(VECTOR<int,2>(Value(id_2),Value(id_1)),skip_counter);}

    void Set_Last_Moved(const int id)
    {if(!(skip_counter&1)) skip_counter++;rigid_body_last_moved(id)=skip_counter;}

    bool Skip_Pair(const int id_1,const int id_2) const
    {return pair_last_checked(Value(id_1)).Get_Default(Value(id_2),0)>rigid_body_last_moved(id_1) && pair_last_checked(Value(id_1)).Get_Default(Value(id_2),0)>rigid_body_last_moved(id_2);}
    //{return pair_last_checked.Get_Default(VECTOR<int,2>(Value(id_1),Value(id_2)),0)>rigid_body_last_moved(id_1) && pair_last_checked.Get_Default(VECTOR<int,2>(Value(id_1),Value(id_2)),0)>rigid_body_last_moved(id_2);}

//#####################################################################
};
}
#endif
