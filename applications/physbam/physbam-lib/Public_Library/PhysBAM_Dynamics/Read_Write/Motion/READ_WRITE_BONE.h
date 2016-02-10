//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BONE
//#####################################################################
#ifndef __READ_WRITE_BONE__
#define __READ_WRITE_BONE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Dynamics/Motion/BONE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<BONE<T>,RW>
{
public:
    static void Read(std::istream& input,BONE<T>& object)
    {Read_Binary<RW>(input,object.transform,object.targeted_transform,object.length);}

    static void Write(std::ostream& output,const BONE<T>& object)
    {Write_Binary<RW>(output,object.transform,object.targeted_transform,object.length);}
};
template<class T> inline std::istream& operator>>(std::istream& input,BONE<T>& bone)
{PHYSBAM_NOT_IMPLEMENTED();}

template<class T> inline std::ostream& operator<<(std::ostream& output,const BONE<T>& bone)
{output<<bone.targeted_transform<<" "<<bone.transform<<" "<<bone.length;return output;}
}
#endif
