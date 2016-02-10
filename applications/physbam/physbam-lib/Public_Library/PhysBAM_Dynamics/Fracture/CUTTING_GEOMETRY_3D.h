//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CUTTING_GEOMETRY_3D__
#define __CUTTING_GEOMETRY_3D__

#include <PhysBAM_Dynamics/Fracture/CUTTING_GEOMETRY.h>
namespace PhysBAM{

template<class TV,int d_input>
class CUTTING_GEOMETRY_3D:public CUTTING_GEOMETRY<TV,d_input>
{
    typedef typename TV::SCALAR T;
    typedef CUTTING_GEOMETRY<TV,d_input> BASE;
public:

    using BASE::intersection_registry;using BASE::cutting_simplices;using BASE::polygon_mesh;using BASE::current_embedding;using BASE::polygons_per_element;using BASE::current_cutting_polygons;
    using BASE::cutting_polygons_per_cutting_simplex;using BASE::simplices_per_current_embedding_simplex;
    using BASE::verbose;
    using BASE::Get_Particles_On_Simplices;using BASE::Get_Polygon_Edges;

    typedef typename BASE::T_CUTTING_SIMPLICES T_CUTTING_SIMPLICES;
    typedef typename BASE::T_CUTTING_SIMPLEX T_CUTTING_SIMPLEX;

    CUTTING_GEOMETRY_3D(const bool verbose)
        :CUTTING_GEOMETRY<TV,d_input>(verbose)
    {}

//#####################################################################
protected:
    void Get_Simplex_Weights_For_Edge_Triangle_Intersection(const VECTOR<int,3>& simplices,const int triangle_array_index,const VECTOR<int,2>& shared_edge,
        VECTOR<VECTOR<T,2>,3>& all_weights);
    void Intersect_Simplex_With_Old_Simplices_In_Embedding(const int embedding_simplex,const int new_simplex) PHYSBAM_OVERRIDE;
    void Split_Existing_Polygon_Edges() PHYSBAM_OVERRIDE;
    void Split_Existing_Polygons() PHYSBAM_OVERRIDE;
    void Get_Unoriented_Segments_On_Two_Simplices(const int simplex_1,const int simplex_2,ARRAY<VECTOR<int,2> >& new_unoriented_segments_on_simplex) const;
    void Divide_Polygon_Particles_With_New_Segments(ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<VECTOR<int,2> >& possible_segments_to_add,const ARRAY<ARRAY<int> >& polygon_particles,
        const int cutting_simplex,const bool flipped,ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const;
    bool Potential_Segment_Should_Be_Added_To_Polygon(const ARRAY<ARRAY<int > >& particles_for_polygon,const bool flipped,const int simplex,const VECTOR<int,2>& nodes) const;
    bool Point_Is_Inside_Unoriented_Polygon(const ARRAY<int>& polygon_particles,const int simplex_owner,const int point) const;
    void Two_Dimensional_Region_Finding_On_Cutting_Simplex(const int cutting_simplex,const bool flipped,const ARRAY<VECTOR<int,2> >& segments,
        ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions,const bool survive_incomplete_regions=false) const;
    void Inside_Outside_Determination_For_Unconnected_Polygonal_Regions(const int simplex,const bool flipped,const ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions,
        ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const;
    void Draw_Polygon(const int simplex,const bool flipped,const ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions) const;
//#####################################################################
};
}
#endif



