#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_LOOKUP_FIRE_RLE
//#####################################################################
#ifndef __FACE_LOOKUP_FIRE_RLE__
#define __FACE_LOOKUP_FIRE_RLE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T> class PROJECTION_RLE;

template<class T_GRID>
class FACE_LOOKUP_FIRE_RLE
{
    typedef typename T_GRID::SCALAR T;
public:
    typedef T ELEMENT;

    FACE_LOOKUP_FIRE_RLE(const ARRAY<T>&,const PROJECTION_RLE<T_GRID>&,const ARRAY<T>*){PHYSBAM_NOT_IMPLEMENTED();}
    const ARRAY<T>& Raw_Data() const {PHYSBAM_NOT_IMPLEMENTED();}

    struct LOOKUP
    {
        typedef T ELEMENT;
        LOOKUP(const FACE_LOOKUP_FIRE_RLE<T_GRID>&,const T){PHYSBAM_NOT_IMPLEMENTED();}
        T operator()(const int,const int) const {PHYSBAM_NOT_IMPLEMENTED();}
    };

    template<class T_FACE> LOOKUP Starting_Point_Face(const T_FACE&) const {PHYSBAM_NOT_IMPLEMENTED();}
    LOOKUP Starting_Point_Cell(const int) const {PHYSBAM_NOT_IMPLEMENTED();}
    LOOKUP Region(const T) const {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
};
}
#endif
#endif 
