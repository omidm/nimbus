//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RIGID_BODY_MASS
//#####################################################################
#ifndef __READ_WRITE_RIGID_BODY_MASS__
#define __READ_WRITE_RIGID_BODY_MASS__

#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_0X0.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_1X1.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
namespace PhysBAM{

template<class RW,class TV,bool world_space>
class Read_Write<RIGID_BODY_MASS<TV,world_space>,RW>
{
public:
    static void Read(std::istream& input,RIGID_BODY_MASS<TV,world_space>& object)
    {Read_Binary<RW>(input,object.mass,object.inertia_tensor);PHYSBAM_ASSERT(object.Valid());}

    static void Write(std::ostream& output,const RIGID_BODY_MASS<TV,world_space>& object)
    {Write_Binary<RW>(output,object.mass,object.inertia_tensor);}
};
template<class TV> inline std::istream& operator>>(std::istream& input,RIGID_BODY_MASS<TV>& mass)
{input>>mass.mass>>mass.inertia_tensor;return input;}

template<class TV> inline std::ostream& operator<<(std::ostream& output,const RIGID_BODY_MASS<TV>& mass)
{output<<mass.mass<<" "<<mass.inertia_tensor;return output;}
}
#endif
