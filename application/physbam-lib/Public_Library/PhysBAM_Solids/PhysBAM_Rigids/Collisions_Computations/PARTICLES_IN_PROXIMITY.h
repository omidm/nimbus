//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace PARTICLES_IN_PROXIMITY
//##################################################################### 
#ifndef __PARTICLES_IN_PROXIMITY__
#define __PARTICLES_IN_PROXIMITY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/PARTICLE_LEVELSET_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM{

namespace PARTICLES_IN_PROXIMITY
{
//#####################################################################
// Function Normal_Is_Admissible
//#####################################################################
template<class T>
class CONVEX_HULL_COMPARATOR
{
public:
    bool operator()(const PAIR<T,int>& p0,const PAIR<T,int>& p1)
    {
        return p0.x<p1.x;
    };
};
template<class T>
ARRAY<int> Convex_Hull(ARRAY<VECTOR<T,2> >& points)
{
    typedef VECTOR<T,2> TV;

    ARRAY<int> hull;

    int corner=1;
    for(int i=2;i<=points.m;i++)
        if(points(i)(2)<points(corner)(2) || (points(i)(2)==points(corner)(2) && points(i)(1)<points(corner)(1)))
            corner=i;

    ARRAY<PAIR<T,int> > angles(points.m);
    angles(1).x=(T)-pi;
    angles(1).y=1;

    /*for(int i=1;i<=angles.m;i++)
        {std::stringstream ss;ss << "angle " << i << " " << angles(i).x << " " << angles(i).y << std::endl;LOG::filecout(ss.str());}*/

    for(int i=2;i<=points.m;i++)
    {
        TV direction=points(i)-points(corner);
        T angle=atan2(direction(2),direction(1));
        angles(i)=PAIR<T,int>(angle,i);
    }
    CONVEX_HULL_COMPARATOR<T> comparator;
    Sort(angles,comparator);

    /*for(int i=1;i<=angles.m;i++)
      {std::stringstream ss;ss << "angle " << i << " " << angles(i).x << " " << angles(i).y << std::endl;LOG::filecout(ss.str());}*/

    hull.Append(corner);
    int i=2;
    while(i<points.m)
    {
        //{std::stringstream ss;ss << i << " hull " << hull << std::endl;LOG::filecout(ss.str());}
        if(hull.m<2)
        {
            hull.Append(angles(i).y);
            i++;
        }
        else
        {
            TV last_direction_perpendicular=(points(hull(hull.m))-points(hull(hull.m-1))).Perpendicular();
            TV direction=points(angles(i).y)-points(hull(hull.m-1));
            if(TV::Dot_Product(last_direction_perpendicular,direction)>0)
            {
                hull.Append(angles(i).y);
                i++;
            }
            else
                hull.Pop();
        }
    }
    
    return hull;
}
//#####################################################################
// Function Normal_Is_Admissible
//#####################################################################
template<class T>
bool Normal_Is_Admissible(RIGID_BODY<VECTOR<T,1> >& body,int particle,VECTOR<T,1>& normal,T tolerance)
{
    return false;
}
template<class T>
bool Normal_Is_Admissible(RIGID_BODY<VECTOR<T,2> >& body,int particle,VECTOR<T,2>& normal,T tolerance)
{
    return false;
}
template<class T>
bool Normal_Is_Admissible(RIGID_BODY<VECTOR<T,3> >& body,int particle,VECTOR<T,3>& normal,T tolerance)
{
    //TODO:add support for non manifold/closed meshes.
    //TODO:this does not actually work in general, need to find something that can actually determine if a vector is a positive combination of incident face normals
    //     currently rejects acceptable normals

    typedef VECTOR<T,3> TV;

    TRIANGULATED_SURFACE<T>* surface=body.simplicial_object;
    TRIANGLE_MESH& mesh=surface->mesh;

    if(!mesh.topologically_sorted_incident_elements)
        mesh.Initialize_Topologically_Sorted_Incident_Elements();

    ARRAY<int>& incident_elements=mesh.topologically_sorted_incident_elements->operator()(particle);

    TV normal_residual=body.Rotation().Inverse_Rotate(normal);

    for(int i=1;i<=incident_elements.m;i++)
    {
        TV element_normal=surface->Face_Normal(incident_elements(i));
        T dot=TV::Dot_Product(element_normal,normal_residual);
        
        if(dot<0)
            normal_residual-=element_normal*dot;
    }

    return normal_residual.Magnitude()<tolerance;
}
//#####################################################################
// Function Particles_In_Proximity
//#####################################################################
template<class TV>
void Particles_In_Proximity(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<TV>& locations,ARRAY<TV>& normals,ARRAY<typename TV::SCALAR>& distances,typename TV::SCALAR proximity_distance)
{
    typedef typename TV::SCALAR T;
    const int d=TV::dimension;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;

    FRAME<TV> frame=levelset_body.Frame().Inverse_Times(particle_body.Frame());
    MATRIX<T,d> rotation=frame.r.Rotation_Matrix();
    VECTOR<T,d> translation=frame.t;

    {std::stringstream ss;ss<<"M: Comparing intersections between "<<particle_body.particle_index<<" and "<<levelset_body.particle_index<<std::endl;LOG::filecout(ss.str());}
    ARRAY_VIEW<TV>& particles_X=particle_body.simplicial_object->particles.X;
    IMPLICIT_OBJECT<VECTOR<T,d> >& object_space_implicit_object=*levelset_body.implicit_object->object_space_implicit_object;
    for(int p=1;p<=particles_X.Size();p++){T distance;
        TV body_2_location=rotation*particles_X(p)+translation;
        if(object_space_implicit_object.Lazy_Inside_Extended_Levelset_And_Value(body_2_location,distance,proximity_distance)){
            TV location=particle_body.World_Space_Point(particles_X(p));
            TV normal=levelset_body.Implicit_Geometry_Normal(location);
            locations.Append(location);
            normals.Append(normal);
            distances.Append(distance);}}
}
//#####################################################################
// Function Get_Maximum_Edge_Length
//#####################################################################
template<class T> T Get_Maximum_Edge_Length(POINT_SIMPLICES_1D<T>& simplicial_object)
{
    return 0;
}
template<class T> T Get_Maximum_Edge_Length(SEGMENTED_CURVE_2D<T>& simplicial_object)
{
    return 0;
}
template<class T> T Get_Maximum_Edge_Length(TRIANGULATED_SURFACE<T>& simplicial_object)
{
    return simplicial_object.Maximum_Edge_Length();
}
//#####################################################################
// Function Build_Particle_Partition
//#####################################################################
template<class TV>
PARTICLE_PARTITION<TV> Build_Particle_Partition(ARRAY<TV>& locations,typename TV::SCALAR region_threshold)
{
    typedef typename TV::SCALAR T;

    int n=locations.m;

    RANGE<TV> box;
    box.Reset_Bounds(locations(1));
    for(int i=1;i<=n;i++)
        box.Enlarge_To_Include_Point(locations(i));
    box=box.Thickened(region_threshold/(T)2.0);
    
    T volume=box.Robust_Size();

    T cell_edge=cbrt(volume/n);
    VECTOR<int,TV::dimension> counts=VECTOR<int,TV::dimension>(box.Edge_Lengths()*((T)1.0/cell_edge))+VECTOR<int,TV::dimension>::All_Ones_Vector();
    
    PARTICLE_PARTITION<TV> partition(box,counts,GEOMETRY_PARTICLES<TV>(),false);
    
    for(int i=1;i<=n;i++)
        partition.Add_To_Partition(locations(i),i);

    return partition;
}
//#####################################################################
// Function Aggregate_And_Smooth_Convex_Regions
//#####################################################################
template<class TV>
void Aggregate_And_Stagger_Convex_Regions(ARRAY<TV>& locations, ARRAY<TV>& normals,ARRAY<typename TV::SCALAR>& distances,typename TV::SCALAR region_threshold)
{
    typedef typename TV::SCALAR T;

    assert(region_threshold>0);

    int n=locations.m;

    for(int i=1;i<=n;i++)
        if(distances(i)<0)
            distances(i)=0;

    /*PARTICLE_PARTITION<TV> partition=Build_Particle_Partition(locations,region_threshold);
    UNION_FIND<int> regions_union(n);

    for(int i=1;i<=n;i++)
    {
        ARRAY<int> candidates;
        partition.Proximity_List(locations(i),region_threshold,candidates);
        if(candidates.m)
        {
            T min_distance=FLT_MAX;
            int min_index=0;
            for(int j=1;j<=candidates.m;j++)
            {
                if(candidates(j)!=i)
                {
                    T distance=(locations(i)-locations(candidates(j))).Magnitude();
                    if(distances(candidates(j))<=distances(i) && distance<region_threshold && distance<min_distance)
                    {
                        min_distance=distance;
                        min_index=candidates(j);
                    }
                }
            }
            if(min_index>0)
            {
                //{std::stringstream ss;ss << "found lower point nearby " << min_distance << std::endl;LOG::filecout(ss.str());}
                regions_union.Union(i,min_index);
            }
        }
    }
    
    ARRAY<T> regions_min_distance(n);
    for(int i=1;i<=n;i++)
        regions_min_distance(i)=distances(i);

    for(int i=1;i<=n;i++)
        if(regions_min_distance(regions_union.Find(i))>distances(i))
            regions_min_distance(regions_union.Find(i))=distances(i);

    for(int i=1;i<=n;i++)
        distances(i)=distances(i)-regions_min_distance(regions_union.Find(i));*/
}
//#####################################################################
// Function Smooth_Nearby_Outlying_Normals
//#####################################################################
/*template<class T>
VECTOR<VECTOR<T,1>,0> Tangent_Space(VECTOR<T,1>& n)
{
    return VECTOR<VECTOR<T,1>,0>();
};
template<class T>
VECTOR<VECTOR<T,2>,1> Tangent_Space(VECTOR<T,2>& n)
{
    return VECTOR<VECTOR<T,2>,1>(n.Unit_Orthogonal_Vector());
};*/
template<class T>
VECTOR<VECTOR<T,3>,2> Tangent_Space(VECTOR<T,3>& n)
{
    VECTOR<T,3> tangent0=n.Unit_Orthogonal_Vector();
    VECTOR<T,3> tangent1=VECTOR<T,3>::Cross_Product(n,tangent0);
    return VECTOR<VECTOR<T,3>,2>(tangent0,tangent1);
};
template<class T>
ARRAY<int> Eliminate_Redundant_Contact_Points(ARRAY<VECTOR<T,1> >& locations,ARRAY<VECTOR<T,1> >& normals,ARRAY<T>& distances)
{
    return ARRAY<int>(IDENTITY_ARRAY<int>(locations.m));
}
template<class T>
ARRAY<int> Eliminate_Redundant_Contact_Points(ARRAY<VECTOR<T,2> >& locations,ARRAY<VECTOR<T,2> >& normals,ARRAY<T>& distances)
{
    return ARRAY<int>(IDENTITY_ARRAY<int>(locations.m));
}
template<class T>
ARRAY<int> Eliminate_Redundant_Contact_Points(ARRAY<VECTOR<T,3> >& locations,ARRAY<VECTOR<T,3> >& normals,ARRAY<T>& distances)
{
    typedef VECTOR<T,3> TV;

    int n=locations.m;
    if(!n)
        return ARRAY<int>();

    T normal_threshold=(T)0.9;
    T distance_threshold=sqrt((T)2.0-(T)2.0*normal_threshold);

    //group contact points into nearby convex regions with similar normals
    //apply convex hull algorithm to regions
    //keep deepest points in regions
    PARTICLE_PARTITION<TV> partition=Build_Particle_Partition(normals,distance_threshold);
    UNION_FIND<int> regions_union(n);
    
    for(int i=1;i<=n;i++)
    {
        if(regions_union.Find(i)==i)
        {
            ARRAY<int> candidates;
            partition.Proximity_List(normals(i),distance_threshold,candidates);
            //{std::stringstream ss;ss << "candidates " << i << " " << candidates << std::endl;LOG::filecout(ss.str());}
            for(int j=1;j<=candidates.m;j++)
                if(TV::Dot_Product(normals(i),normals(candidates(j)))>normal_threshold)
                    regions_union.Union(i,candidates(j));
        }
    }

    ARRAY<ARRAY<int> > regions(n);

    ARRAY<TV> regions_normals(n);
    for(int i=1;i<=n;i++)
        regions(regions_union.Find(i)).Append(i);

    //{std::stringstream ss;ss << "regions" << std::endl << regions << std::endl;LOG::filecout(ss.str());}

    ARRAY<int> points;

    for(int i=1;i<=n;i++)
        if(regions(i).m)
        {
            ARRAY<int>& region=regions(i);
            TV normal;
            for(int j=1;j<region.m;j++)
                normal+=normals(region(j));
            normal.Normalize();
            VECTOR<TV,TV::dimension-1> tangents=Tangent_Space(normal);
            ARRAY<VECTOR<T,2> > projected_offsets(region.m);
            for(int k=1;k<=region.m;k++)
            {
                TV offset=locations(region(k))-locations(region(1));
                for(int c=1;c<TV::dimension;c++)
                    projected_offsets(k)(c)=TV::Dot_Product(offset,tangents(c));
            }
            ARRAY<int> convex_hull=Convex_Hull(projected_offsets);
            
            int deepest=1;
            for(int k=2;k<=region.m;k++)
                if(distances(region(k))<distances(region(deepest)))
                    deepest=k;
            bool deepest_in_hull=false;
            for(int k=1;k<=convex_hull.m;k++)
            {
                points.Append(region(convex_hull(k)));
                if(convex_hull(k)==deepest)
                    deepest_in_hull=true;
            }
            if(!deepest_in_hull)
                points.Append(region(deepest));
        }

    return points;
}
//#####################################################################
// Function All_Particles_In_Proximity
//#####################################################################
template<class TV>
void All_Particles_In_Proximity(RIGID_BODY<TV>& body_1,RIGID_BODY<TV>& body_2,ARRAY<TV>& locations,ARRAY<TV>& normals,ARRAY<typename TV::SCALAR>& distances,typename TV::SCALAR proximity_distance,const bool stagger_points=true)
{
    typedef typename TV::SCALAR T;

    ARRAY<TV> locations_1;
    ARRAY<TV> normals_1;
    ARRAY<typename TV::SCALAR> distances_1;
    ARRAY<TV> locations_2;
    ARRAY<TV> normals_2;
    ARRAY<typename TV::SCALAR> distances_2;

    locations.Remove_All();
    normals.Remove_All();
    distances.Remove_All();

    T region_threshold_1=0,region_threshold_2=0;
    if(body_1.simplicial_object && body_2.implicit_object){
        Particles_In_Proximity(body_1,body_2,locations_1,normals_1,distances_1,proximity_distance);
        if(locations_1.m){
            if(stagger_points){
                region_threshold_1=Get_Maximum_Edge_Length(*body_1.simplicial_object);
                Aggregate_And_Stagger_Convex_Regions(locations_1,normals_1,distances_1,region_threshold_1);}
            for(int i=1;i<=normals_1.m;i++) normals_1(i)=-normals_1(i);}}
    if(body_1.implicit_object && body_2.simplicial_object){
        Particles_In_Proximity(body_2,body_1,locations_2,normals_2,distances_2,proximity_distance);
        if(locations_2.m){
            if(stagger_points){
                region_threshold_2=Get_Maximum_Edge_Length(*body_2.simplicial_object);
                Aggregate_And_Stagger_Convex_Regions(locations_2,normals_2,distances_2,region_threshold_2);}}}

    ARRAY<TV> locations_merged;
    ARRAY<TV> normals_merged;
    ARRAY<T> distances_merged;

    locations_merged.Append_Elements(locations_1);
    locations_merged.Append_Elements(locations_2);
    normals_merged.Append_Elements(normals_1);
    normals_merged.Append_Elements(normals_2);
    distances_merged.Append_Elements(distances_1);
    distances_merged.Append_Elements(distances_2);

    ARRAY<int> pruned_set=Eliminate_Redundant_Contact_Points(locations_merged,normals_merged,distances_merged);
    for(int i=1;i<=pruned_set.m;i++){
        locations.Append(locations_merged(pruned_set(i)));
        normals.Append(normals_merged(pruned_set(i)));
        distances.Append(distances_merged(pruned_set(i)));}
}

}
}

#endif
