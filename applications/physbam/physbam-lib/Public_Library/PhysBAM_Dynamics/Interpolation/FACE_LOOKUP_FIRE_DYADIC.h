#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2005, Ron Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FACE_LOOKUP_FIRE_DYADIC__
#define __FACE_LOOKUP_FIRE_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_FORWARD.h>
namespace PhysBAM{

template<class T_GRID>
class FACE_LOOKUP_FIRE_DYADIC
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T TV;typedef PROJECTION_DYADIC<T_GRID> T_PROJECTION;
public:
    typedef T ELEMENT;
        
    const ARRAY<T>& V_face;
protected:
    const T_PROJECTION& projection;
    const ARRAY<T>* phi;
    const ARRAY<T>* phi_face;
public:

    FACE_LOOKUP_FIRE_DYADIC(const ARRAY<T>& V_face_input,const T_PROJECTION& projection_input,const ARRAY<T>* phi_input=0,const ARRAY<T>* phi_face_input=0)
        :V_face(V_face_input),projection(projection_input),phi(phi_input),phi_face(phi_face_input)
    {}

    const ARRAY<T>& Raw_Data() const
    {return V_face;}

    class LOOKUP
    {
    private:
        const FACE_LOOKUP_FIRE_DYADIC<T_GRID>& face_lookup;
        const T reference_region;
    public:
        typedef T ELEMENT;

        LOOKUP(const FACE_LOOKUP_FIRE_DYADIC<T_GRID>& face_lookup_input,const T reference_region_input)
            :face_lookup(face_lookup_input),reference_region(reference_region_input)
        {}

        T operator()(const int face) const
        {return face_lookup.projection.Face_Velocity_With_Ghost_Value(face_lookup.V_face,face,reference_region);}
    };

    LOOKUP Starting_Point_Face(const int face) const
    {assert(phi_face);return LOOKUP(*this,(*phi_face)(face));}
    
    LOOKUP Starting_Point_Cell(const int cell) const
    {assert(phi);return LOOKUP(*this,(*phi)(cell));}
    
    LOOKUP Region(const T phi) const
    {return LOOKUP(*this,(*phi));}

//#####################################################################
};
}
#endif
#endif
