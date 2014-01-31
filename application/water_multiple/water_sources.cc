#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/water_example.h"
#include "application/water_alternate_fine/water_sources.h"
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

namespace {
    typedef application::T T;
} // namespace

WaterSources::WaterSources() {};

void WaterSources::Add_Source(WATER_EXAMPLE<VECTOR<T,1> > *example)
{
    PHYSBAM_FATAL_ERROR();
}

void WaterSources::Add_Source(WATER_EXAMPLE<VECTOR<T,2> > *example)
{
    typedef VECTOR<T,2> TV;
    TV point1,point2;BOX<TV> source;
    point1=TV::All_Ones_Vector()*(T).5;point1(1)=.95;point1(2)=.6;point2=TV::All_Ones_Vector()*(T).65;point2(1)=1;point2(2)=.75;
    source.min_corner=point1;source.max_corner=point2;
    example->sources.Append(new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(source));
}

void WaterSources::Add_Source(WATER_EXAMPLE<VECTOR<T,3> > *example)
{
    typedef VECTOR<T,3> TV;
    TV point1,point2;CYLINDER<T> source;
    point1=TV::All_Ones_Vector()*(T).8;point1(1)=.4;point1(3)=.95;point2=TV::All_Ones_Vector()*(T).8;point2(1)=.4;point2(3)=1;
    source.Set_Endpoints(point1,point2);source.radius=.1;
    IMPLICIT_OBJECT<TV>* analytic=new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(source);
    example->sources.Append(analytic);
}
