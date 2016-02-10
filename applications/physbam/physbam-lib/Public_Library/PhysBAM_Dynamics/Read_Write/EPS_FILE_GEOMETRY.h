//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EPS_FILE_GEOMETRY__
#define __EPS_FILE_GEOMETRY__

#include <PhysBAM_Tools/Images/EPS_FILE.h>
namespace PhysBAM{

template<class T> class TRIANGLE_2D;
template<class TV> class SPHERE;

template<class T>
class EPS_FILE_GEOMETRY:public EPS_FILE<T>
{
    typedef VECTOR<T,2> TV;
    typedef EPS_FILE<T> BASE;
public:
    using BASE::Line_Color;using BASE::stream;

    EPS_FILE_GEOMETRY(const std::string& filename,const RANGE<TV>& box=RANGE<TV>(TV(),TV(500,500)))
        :EPS_FILE<T>(filename,box)
    {}

    virtual ~EPS_FILE_GEOMETRY()
    {}

    template<class T_OBJECT>
    void Draw_Object_Colored(const T_OBJECT& object,const VECTOR<T,3>& color)
    {(*stream)<<"gsave ";Line_Color(color);Draw_Object(object);(*stream)<<"grestore"<<std::endl;}

    void Draw_Object(const TRIANGLE_2D<T>& tri)
    {(*stream)<<tri.X[1]<<" moveto "<<tri.X[2]<<" lineto "<<tri.X[3]<<" lineto closepath stroke"<<std::endl;Bound(tri.X);}

    void Draw_Object(const SPHERE<TV>& circle)
    {(*stream)<<circle.center<<" "<<circle.radius<<" 0 360 arc stroke"<<std::endl;Bound(circle.center-circle.radius);Bound(circle.center+circle.radius);}


//#####################################################################
};
}
#endif
