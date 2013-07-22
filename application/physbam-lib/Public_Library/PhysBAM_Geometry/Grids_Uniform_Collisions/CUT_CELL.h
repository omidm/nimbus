//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUT_CELL
//#####################################################################
#ifndef __CUT_CELL__
#define __CUT_CELL__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
namespace PhysBAM {
template<class T,int d>
class CUT_CELLS
{
  public:
    CUT_CELLS() : dominant_element(0),geometry(0),visibility(0),visibility_nodes(0)
    {}

    int dominant_element;
    ARRAY<POLYGON<VECTOR<T,d> > > geometry;
    ARRAY<ARRAY<VECTOR<int,d> > > visibility;
    ARRAY<ARRAY<VECTOR<T,d> > > visibility_nodes;

    static void Print_Debug_Information(const POLYGON<VECTOR<T,d> >& geometry)
    {std::stringstream ss;ss<<"GEOMETRY DEBUG INFO: "<<geometry.X.Size()<<std::endl;
        for(int j=1;j<=geometry.X.Size();++j) ss<<"GEOMETRY DEBUG INFO: "<<geometry.X(j)<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<geometry.X(1)<<std::endl;
        LOG::filecout(ss.str());}

    static void Print_Debug_Information(const ARRAY<POLYGON<VECTOR<T,d> > >& geometry)
    {std::stringstream ss;for(int i=1;i<=geometry.Size();++i){
        Print_Debug_Information(geometry(i));
        ss<<"GEOMETRY DEBUG INFO: "<<std::endl;}
        LOG::filecout(ss.str());}

    static void Print_Debug_Information(const ARRAY<POLYGON<VECTOR<T,1> > >& geometry,const RANGE<VECTOR<T,1> >& bounding_box){
        std::stringstream ss;
        ss<<"GEOMETRY DEBUG INFO: "<<(bounding_box.min_corner.x+bounding_box.max_corner.x)*(T).5<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.max_corner.x<<std::endl;
        LOG::filecout(ss.str());
        Print_Debug_Information(geometry);}

    static void Print_Debug_Information(const ARRAY<POLYGON<VECTOR<T,2> > >& geometry,const RANGE<VECTOR<T,2> >& bounding_box){
        std::stringstream ss;
        ss<<"GEOMETRY DEBUG INFO: "<<(bounding_box.min_corner.x+bounding_box.max_corner.x)*(T).5<<"\t"
            <<(bounding_box.min_corner.y+bounding_box.max_corner.y)*(T).5<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.min_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.max_corner.x<<"\t"<<bounding_box.min_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.max_corner.x<<"\t"<<bounding_box.max_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.max_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.min_corner.y<<std::endl;
        LOG::filecout(ss.str());
        Print_Debug_Information(geometry);}

    static void Print_Debug_Information(const ARRAY<POLYGON<VECTOR<T,3> > >& geometry,const RANGE<VECTOR<T,3> >& bounding_box){
        std::stringstream ss;
        ss<<"GEOMETRY DEBUG INFO: "<<(bounding_box.min_corner.x+bounding_box.max_corner.x)*(T).5<<"\t"
            <<(bounding_box.min_corner.y+bounding_box.max_corner.y)*(T).5<<"\t"<<(bounding_box.min_corner.z+bounding_box.max_corner.z)*(T).5<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.min_corner.y<<"\t"<<bounding_box.min_corner.z<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.max_corner.x<<"\t"<<bounding_box.min_corner.y<<"\t"<<bounding_box.min_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.max_corner.x<<"\t"<<bounding_box.max_corner.y<<"\t"<<bounding_box.min_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.max_corner.y<<"\t"<<bounding_box.min_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.max_corner.y<<"\t"<<bounding_box.max_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.max_corner.x<<"\t"<<bounding_box.max_corner.y<<"\t"<<bounding_box.max_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.max_corner.x<<"\t"<<bounding_box.min_corner.y<<"\t"<<bounding_box.max_corner.y<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.min_corner.y<<"\t"<<bounding_box.max_corner.z<<std::endl;
        ss<<"GEOMETRY DEBUG INFO: "<<bounding_box.min_corner.x<<"\t"<<bounding_box.min_corner.y<<"\t"<<bounding_box.min_corner.z<<std::endl;
        LOG::filecout(ss.str());
        Print_Debug_Information(geometry);}
};
}
#endif
