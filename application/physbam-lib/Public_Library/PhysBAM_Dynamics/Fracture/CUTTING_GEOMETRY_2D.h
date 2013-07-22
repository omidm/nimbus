//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CUTTING_GEOMETRY_2D__
#define __CUTTING_GEOMETRY_2D__

#include <PhysBAM_Dynamics/Fracture/CUTTING_GEOMETRY.h>
namespace PhysBAM{

template<class TV,int d_input>
class CUTTING_GEOMETRY_2D:public CUTTING_GEOMETRY<TV,d_input>
{
    typedef typename TV::SCALAR T;
    typedef CUTTING_GEOMETRY<TV,d_input> BASE;
public:

    using BASE::intersection_registry;using BASE::cutting_simplices;using BASE::polygon_mesh;using BASE::current_embedding;using BASE::polygons_per_element;using BASE::current_cutting_polygons;
    using BASE::cutting_polygons_per_cutting_simplex;using BASE::simplices_per_current_embedding_simplex;
    using BASE::verbose;
    using BASE::Get_Particles_On_Simplices;using BASE::Get_Polygon_Edges;

    CUTTING_GEOMETRY_2D(const bool verbose)
        :CUTTING_GEOMETRY<TV,d_input>(verbose)
    {}

//#####################################################################
    void Intersect_Simplex_With_Old_Simplices_In_Embedding(const int embedding_simplex,const int new_simplex) PHYSBAM_OVERRIDE;
    void Split_Existing_Polygon_Edges() PHYSBAM_OVERRIDE;
    void Split_Existing_Polygons() PHYSBAM_OVERRIDE;
    void Divide_Polygon_Particles_With_New_Segments(ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<int>& possible_particles_to_add,const ARRAY<ARRAY<int> >& polygon_particles,
        const int cutting_simplex,const bool flipped,ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const;
//#####################################################################
};
}
#endif
