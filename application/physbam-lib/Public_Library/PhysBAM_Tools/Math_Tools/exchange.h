//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function exchange  
//#####################################################################
//
// exchanges the values of a and b  
//               
//#####################################################################
#ifndef __exchange__
#define __exchange__

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class ENABLE=void> struct EXCHANGE_HELPER
{
    STATIC_ASSERT(!IS_ARRAY<T>::value);

    static void Exchange(T& a,T& b)
    {T c=a;a=b;b=c;}
};

template<class T> void exchange(T& a,T& b)
{EXCHANGE_HELPER<T>::Exchange(a,b);}

template<class T, class ENABLE=void> struct NIMBUS_COPY_HELPER
{
    STATIC_ASSERT(!IS_ARRAY<T>::value);

    static void Nimbus_Copy(T& dest, T& src)
    {dest=src;}
};

template<class T> void nimbus_copy(T& dest,T& src)
{NIMBUS_COPY_HELPER<T>::Nimbus_Copy(dest,src);}

}
#endif
