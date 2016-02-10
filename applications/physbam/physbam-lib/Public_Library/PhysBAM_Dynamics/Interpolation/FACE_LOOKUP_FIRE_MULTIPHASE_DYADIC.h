#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC__
#define __FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYADIC.h>
namespace PhysBAM{

template<class T_GRID>
class FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC
{
    typedef typename T_GRID::SCALAR T;
public:
    typedef T ELEMENT;

    const ARRAY<T>& V_face;

    FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC(const ARRAY<T>&,const PROJECTION_DYADIC<T_GRID>&,const ARRAY<T>*){PHYSBAM_NOT_IMPLEMENTED();}
    const ARRAY<T>& Raw_Data() const {PHYSBAM_NOT_IMPLEMENTED();}

    struct LOOKUP
    {
        typedef T ELEMENT;
        LOOKUP(const FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC<T_GRID>&,const T){PHYSBAM_NOT_IMPLEMENTED();}
        T operator()(const int) const {PHYSBAM_NOT_IMPLEMENTED();}
    };

    template<class T_FACE> LOOKUP Starting_Point_Face(const T_FACE&) const {PHYSBAM_NOT_IMPLEMENTED();}
    LOOKUP Starting_Point_Cell(const int) const {PHYSBAM_NOT_IMPLEMENTED();}
    LOOKUP Region(const T) const {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
};
}
#endif
#endif
