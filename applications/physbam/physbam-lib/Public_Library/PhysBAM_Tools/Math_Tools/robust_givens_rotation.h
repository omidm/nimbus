//#####################################################################
// Copyright 2012, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function robust_givens_rotation
//#####################################################################
//
// Computes the parameters for a Givens-type rotation robustly 
// See LAPACK Working Note 150 by Anderson
//               
//#####################################################################
#ifndef __robust_givens_rotation__
#define __robust_givens_rotation__

#include <cmath>

namespace PhysBAM{

template<class T>
inline void robust_givens_rotation(T& c, T&s, T& r, const T f, const T g)
{if(g==0){c=sign(f);s=0;r=abs(f);}
 else if(f==0){c=0;s=sign(g);r=abs(g);}
 else if(abs(f)>abs(g)){T t=g/f; T u=sign(f)*sqrt(1+t*t);c=1/u;s=t*c;r=f*u;}
 else {T t=f/g;T u=sign(g)*sqrt(1+t*t);s=1/u;c=t*s;r=g*u;}
 return;}
}
#endif

