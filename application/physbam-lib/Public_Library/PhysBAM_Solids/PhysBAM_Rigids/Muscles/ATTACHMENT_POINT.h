//#####################################################################
// Copyright 2005-2007, Kevin Der, Ron Fedkiw, Eran Guendelman, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ATTACHMENT_POINT
//#####################################################################
#ifndef __ATTACHMENT_POINT__
#define __ATTACHMENT_POINT__

namespace PhysBAM{

template<class TV> class RIGID_BODY;

template<class TV>
class ATTACHMENT_POINT
{
public:
    RIGID_BODY<TV>& rigid_body;
    TV object_space_position;

    ATTACHMENT_POINT(RIGID_BODY<TV>& rigid_body_input,TV object_space_position_input);
    TV Embedded_Position() const;
    TV Embedded_Velocity() const;
    void Apply_Impulse(const TV& impulse);
};
}
#endif
