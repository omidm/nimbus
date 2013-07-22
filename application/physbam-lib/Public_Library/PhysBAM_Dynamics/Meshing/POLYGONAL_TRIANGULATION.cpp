//#####################################################################
// Copyright 2006-2008, Kevin Der, Jeffrey Hellrung, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGONAL_TRIANGULATION
//##################################################################### 
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_SIGNED_VOLUME.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/EXACT_SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Dynamics/Meshing/POLYGONAL_TRIANGULATION.h>
using namespace PhysBAM;
//#####################################################################
// Function Triangulate_Convex_Planar_Polygon
//#####################################################################
// robust to non-convexity due to numerical error
template<class T> void POLYGONAL_TRIANGULATION<T>::
Triangulate_Convex_Planar_Polygon(const ARRAY<VECTOR<T,2> >& positions,ARRAY<VECTOR<int,2> > segments,ARRAY<VECTOR<int,3> >& triangles)
{
    HASHTABLE<VECTOR<int,2> > segment_hash;
    ARRAY<ARRAY<int> > incoming_segments(positions.m);ARRAY<ARRAY<int> > outgoing_segments(positions.m);
    ARRAY<ARRAY<int> > all_incoming_segments(positions.m);ARRAY<ARRAY<int> > all_outgoing_segments(positions.m);
    // initiailize: add existing segments
    ARRAY<bool> used(segments.m);
    for(int i=1;i<=segments.m;i++){segment_hash.Insert(segments(i)); // TODO: probably insert sorted
        outgoing_segments(segments(i)(1)).Append(i);incoming_segments(segments(i)(2)).Append(i);
        all_outgoing_segments(segments(i)(1)).Append(i);all_incoming_segments(segments(i)(2)).Append(i);}
    while(used.Contains(false)){assert(used.m==segments.m);
        // make all possible triangles
        for(int segment_id_1=1;segment_id_1<=segments.m;segment_id_1++) if(!used(segment_id_1)){
            VECTOR<int,2> segment_1=segments(segment_id_1);int node_1=segment_1(1);int node_2=segment_1(2);
            ARRAY<int>& incoming_segs_to_1=incoming_segments(node_1);ARRAY<int>& outgoing_segs_from_2=outgoing_segments(node_2);
            bool stop=false;
            for(int i=1;i<=incoming_segs_to_1.m;i++) if(!stop) for(int j=1;j<=outgoing_segs_from_2.m;j++)
                if(!used(incoming_segs_to_1(i))&&!used(outgoing_segs_from_2(j))&&segments(incoming_segs_to_1(i))(1)==segments(outgoing_segs_from_2(j))(2)){
                    int segment_2_id=outgoing_segs_from_2(j);int segment_3_id=incoming_segs_to_1(i);int node_3=segments(segment_3_id)(1);
                    VECTOR<int,3> nodes(node_1,node_2,node_3);VECTOR<VECTOR<T,2>,3> triangle_X(positions.Subset(nodes));
                    if(TRIANGLE_2D<T>::Signed_Area(triangle_X[1],triangle_X[2],triangle_X[3])<0) continue;
                    // ensure no points lie inside the triangle
                    bool point_inside=false;
                    for(int k=1;k<=positions.m;k++) if(!nodes.Contains(k))
                        if(SIMPLEX_INTERACTIONS<T>::Intersection(triangle_X,positions(k))){
                            point_inside=true;break;}
                        if(!point_inside){
                            triangles.Append(nodes);
                            used(segment_id_1)=used(segment_2_id)=used(segment_3_id)=true;
                            incoming_segments(node_1).Remove_Index(i);outgoing_segments(node_1).Remove_Index(outgoing_segments(node_1).Find(segment_id_1));
                            incoming_segments(node_2).Remove_Index(incoming_segments(node_2).Find(segment_id_1));outgoing_segments(node_2).Remove_Index(j);
                            incoming_segments(node_3).Remove_Index(incoming_segments(node_3).Find(segment_2_id));
                            outgoing_segments(node_3).Remove_Index(outgoing_segments(node_3).Find(segment_3_id));
                            stop=true;break;}}}
        if(!used.Contains(false)) break;
        // add the best segment
        int best_p=-1,best_q=-1;T largest_smallest_relative_angle=0; 
        for(int p=1;p<=positions.m;p++) for(int q=p+1;q<=positions.m;q++) if(!segment_hash.Contains(VECTOR<int,2>(p,q))&&!segment_hash.Contains(VECTOR<int,2>(q,p))){
            bool safe=true;
            for(HASHTABLE_ITERATOR<VECTOR<int,2> > segment_hash_iterator(segment_hash);segment_hash_iterator.Valid();segment_hash_iterator.Next()){
                const VECTOR<int,2>& existing_segment=segment_hash_iterator.Key();
                if(existing_segment(1)==p||existing_segment(2)==p||existing_segment(1)==q||existing_segment(2)==q) continue;
                if(SIMPLEX_INTERACTIONS<T>::Intersection(VECTOR<VECTOR<T,2>,2>(positions.Subset(existing_segment)),VECTOR<VECTOR<T,2>,2>(positions(p),positions(q)))){
                    safe=false;break;}}
            if(!safe) continue;
            // it's safe, now check that it's the best one
            ARRAY<int> segments_to_check(all_incoming_segments(p));
            segments_to_check.Append_Elements(all_outgoing_segments(p));
            segments_to_check.Append_Elements(all_incoming_segments(q));
            segments_to_check.Append_Elements(all_outgoing_segments(q));
            T smallest_angle=Find_Sharpest_Angle(positions,segments,VECTOR<int,2>(p,q),segments_to_check);
            if(smallest_angle>largest_smallest_relative_angle){largest_smallest_relative_angle=smallest_angle;best_p=p;best_q=q;}}
        assert(best_p>0&&best_q>>0);
        // add segment
        segment_hash.Insert(VECTOR<int,2>(best_p,best_q));used.Append(false);used.Append(false);
        int first=segments.Append(VECTOR<int,2>(best_p,best_q));
        all_outgoing_segments(best_p).Append(first);all_incoming_segments(best_q).Append(first);
        outgoing_segments(best_p).Append(first);incoming_segments(best_q).Append(first);
        int second=segments.Append(VECTOR<int,2>(best_q,best_p));
        outgoing_segments(best_q).Append(second);incoming_segments(best_p).Append(second);
        all_outgoing_segments(best_q).Append(second);all_incoming_segments(best_p).Append(second);}
}

//#####################################################################
// Function Triangulate_Nonconvex_Planar_Connected_Polygon
//#####################################################################
// assumes all components are connected and segments are directed
// TODO: abstract out functionality common to these two functions
template<class T> void POLYGONAL_TRIANGULATION<T>::
Triangulate_Nonconvex_Planar_Connected_Polygon(const ARRAY<VECTOR<T,2> >& positions,const ARRAY<ARRAY<VECTOR<int,2> > >& segments_input,ARRAY<VECTOR<int,3> >& triangles)
{
    // ensure no 2 vertices are coincident (in float coordinates)
    typedef VECTOR<float,2> TV_float;

    // TODO: Fix the triangulation (this triangulates as if the polygon was simple and convex)
    const ARRAY<VECTOR<int,2> >& segments=segments_input(1);
    ARRAY<int> vertices;
    for(int s=1;s<=segments.m;s++){
        if(segments(s)[2]!=segments(s%segments.m+1)[1]) PHYSBAM_FATAL_ERROR();
        vertices.Append(segments(s)[1]);}
    if(vertices.m<3) PHYSBAM_FATAL_ERROR();
    for(int s=1;s<=segments.m-2;s++)
        triangles.Append(VECTOR<int,3>(vertices(1),vertices(s+1),vertices(s+2)));

#if 0
    return;
    //HASHTABLE<TV_float> float_vertices;
    //for(int v=1;v<=positions.m;v++){TV_float float_vertex(positions(v));
    //    if(float_vertices.Contains(float_vertex)) PHYSBAM_FATAL_ERROR();
    //    float_vertices.Insert(float_vertex);}
    // ensure no 2 segments are intersecting (in float coordinates)
    //for(int c1=1;c1<=segments_input.m;c1++) for(int s1=1;s1<=segments_input(c1).m;s1++) for(int c2=c1,s2=s1+1;c2<=segments_input.m;c2++) for(;s2<=segments_input(c2).m;s2++){
    //    int i,j,k,l;segments_input(c1)(s1).Get(i,j);segments_input(c2)(s2).Get(k,l);
    //    if(i!=k&&i!=l&&j!=k&&j!=l&&EXACT_SIMPLEX_INTERACTIONS::Exact_Segment_Intersection((TV_float)positions(i),(TV_float)positions(j),(TV_float)positions(k),(TV_float)positions(l)))
    //        PHYSBAM_FATAL_ERROR();}
    // polygon boundaries might still not have proper inclusion or (counter-)clockwise ordering
    ARRAY<VECTOR<int,2> > all_segments;ARRAY<bool> segments_used;
    ARRAY<ARRAY<int> > outgoing_segments(positions.m);ARRAY<ARRAY<int> > incoming_segments(positions.m);
    for(int i=1;i<=segments_input.m;i++) for(int j=1;j<=segments_input(i).m;j++){
        int segment_index=all_segments.Append(segments_input(i)(j));segments_used.Append(false);
        outgoing_segments(segments_input(i)(j)(1)).Append(segment_index);
        incoming_segments(segments_input(i)(j)(2)).Append(segment_index);}
    // make segments to connect holes
    HASHTABLE<VECTOR<int,2> > already_connected_holes;
    for(int hole=2;hole<=segments_input.m;hole++){
        for(int j=1;j<=segments_input(hole).m;j++){
            int first_node_on_hole_segment=segments_input(hole)(j)(1);
            // find other node on a different component
            for(int comp=1;comp<=segments_input.m;comp++) if(hole!=comp&&!already_connected_holes.Contains(VECTOR<int,2>(hole,comp).Sorted())){
                for(int k=1;k<=segments_input(comp).m;k++){int candidate_node=segments_input(comp)(k)(1);
                if(!Segment_Intersects(VECTOR<int,2>(first_node_on_hole_segment,candidate_node),all_segments,positions)){
                    int segment_index=all_segments.Append(VECTOR<int,2>(first_node_on_hole_segment,candidate_node));segments_used.Append(false);
                    outgoing_segments(first_node_on_hole_segment).Append(segment_index);
                    incoming_segments(candidate_node).Append(segment_index);
                    segment_index=all_segments.Append(VECTOR<int,2>(candidate_node,first_node_on_hole_segment));segments_used.Append(false);
                    outgoing_segments(candidate_node).Append(segment_index);
                    incoming_segments(first_node_on_hole_segment).Append(segment_index);
                    already_connected_holes.Insert(VECTOR<int,2>(hole,comp).Sorted());
                    goto go_to_next_hole;}}}}go_to_next_hole:{}}
    // make triangles
    int count=0;
    while(count<all_segments.m){
        int globally_best_first_segment_index=-1;
        int globally_best_second_segment_index=-1;
        T best_angle=(T)0;
        bool found=false;
        for(int i=1;i<=segments_used.m;i++)
            if(!segments_used(i)){
                const VECTOR<int,2> first_segment=all_segments(i);
                ARRAY<int>& continuing_segments=outgoing_segments(first_segment(2));
                T min_value_found=1e2;
                int best_second_segment=-1;
                // check angle of the ear
                for(int j=1;j<=continuing_segments.m;j++)
                    if(!segments_used(continuing_segments(j))){
                        const VECTOR<int,2> second_segment=all_segments(continuing_segments(j));
                        T angle=0;
                        if(first_segment(1)==second_segment(2))
                            angle=99.;
                        else
                            angle=VECTOR<T,2>::Oriented_Angle_Between(positions(second_segment(2))-positions(second_segment(1)),
                                                                      positions(first_segment(2))-positions(first_segment(1)));
                        if(angle<min_value_found){
                            min_value_found=angle;
                            best_second_segment=continuing_segments(j);
                        }
                    }
                assert(best_second_segment>0);
                if(min_value_found>0)
                    continue;
                const VECTOR<int,2> second_segment=all_segments(best_second_segment);
                min_value_found=abs(min_value_found+(T)pi); // converting to the positive angle for the ear
                ARRAY<int> segments_to_check(outgoing_segments(first_segment(1)));
                segments_to_check.Append_Elements(incoming_segments(first_segment(1)));
                segments_to_check.Append_Elements(outgoing_segments(second_segment(2)));
                segments_to_check.Append_Elements(incoming_segments(second_segment(2)));
                min_value_found=min(min_value_found,Find_Sharpest_Angle(positions,all_segments,VECTOR<int,2>(second_segment(2),first_segment(1)),segments_to_check));
                // check no points are inside the triangle to be formed
                bool is_an_ear=true;
                for(int j=1;j<=positions.m;j++)
                    if(!first_segment.Contains(j)&&j!=second_segment(2))
                        if(SIMPLEX_INTERACTIONS<T>::Intersection(VECTOR<VECTOR<T,2>,3>(positions(first_segment(1)),positions(first_segment(2)),positions(second_segment(2))),positions(j))){
                            is_an_ear=false;
                            break;
                        }
                if(is_an_ear)
                    if(min_value_found>best_angle){
                        best_angle=min_value_found;
                        globally_best_first_segment_index=i;
                        globally_best_second_segment_index=best_second_segment;
                        found=true;
                    }
            }
            assert(found);
            VECTOR<int,2> first_segment=all_segments(globally_best_first_segment_index);
            VECTOR<int,2> second_segment=all_segments(globally_best_second_segment_index);
            int third_segment_index=all_segments.Find(VECTOR<int,2>(second_segment(2),first_segment(1)));
            if(!third_segment_index){
                int segment_index=all_segments.Append(VECTOR<int,2>(first_segment(1),second_segment(2)));segments_used.Append(false);
                outgoing_segments(first_segment(1)).Append(segment_index);
                incoming_segments(second_segment(2)).Append(segment_index);
                segment_index=all_segments.Append(VECTOR<int,2>(second_segment(2),first_segment(1)));segments_used.Append(true);
                outgoing_segments(second_segment(2)).Append(segment_index);
                incoming_segments(first_segment(1)).Append(segment_index);
            }
            else{
                segments_used(third_segment_index)=true;
            }
            segments_used(globally_best_first_segment_index)=true;
            segments_used(globally_best_second_segment_index)=true;
            count+=3;
            triangles.Append(VECTOR<int,3>(first_segment(1),first_segment(2),second_segment(2)));
    }
#endif
}
//#####################################################################
// Function Segment_Intersects
//#####################################################################
template<class T> bool POLYGONAL_TRIANGULATION<T>::
Segment_Intersects(const VECTOR<int,2>& candidate_segment,const ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<VECTOR<T,2> >& positions)
{
    bool any_intersects=false;
    for(int i=1;i<=all_segments.m;i++){const VECTOR<int,2>& test_segment=all_segments(i);
    if(test_segment.Contains(candidate_segment(1))||test_segment.Contains(candidate_segment(2))) continue;
    bool intersects=SIMPLEX_INTERACTIONS<T>::Intersection(VECTOR<VECTOR<T,2>,2>(positions.Subset(candidate_segment)),VECTOR<VECTOR<T,2>,2>(positions.Subset(test_segment)));
    if(intersects){any_intersects=true;break;}}
    return any_intersects;
}
//#####################################################################
// Function Find_Sharpest_Angle
//#####################################################################
template<class T> T POLYGONAL_TRIANGULATION<T>::
Find_Sharpest_Angle(const ARRAY<VECTOR<T,2> >& positions,const ARRAY<VECTOR<int,2> >& all_segments,const VECTOR<int,2>& given_segment,const ARRAY<int>& test_segments)
{
    T min_value_found=(T)two_pi;
    for(int j=1;j<=test_segments.m;j++){
        const VECTOR<int,2>& test_segment=all_segments(test_segments(j));
        if(test_segment==given_segment||(test_segment(1)==given_segment(2)&&test_segment(2)&&given_segment(1))) continue;
        T angle=abs(VECTOR<T,2>::Oriented_Angle_Between(positions(given_segment(2))-positions(given_segment(1)),
            positions(test_segment(2))-positions(test_segment(1))));
        if(angle>(T)pi*.5) angle=(T)pi-angle;if(angle<min_value_found) min_value_found=angle;}
    return min_value_found;
}
//#####################################################################
// Function Orientation
//#####################################################################
template<class T> int POLYGONAL_TRIANGULATION<T>::
Orientation(const VECTOR<T,2>& x1,const VECTOR<T,2>& x2,const VECTOR<T,2>& x3)
{
    return Adaptive_Signed_Volume< EXACT_FLOAT<T> >(x1,x2,x3).Sign();
}
//#####################################################################
// Function Point_Triangle_Intersects
//#####################################################################
template<class T> bool POLYGONAL_TRIANGULATION<T>::
Point_Triangle_Intersects(const VECTOR<T,2>& y,const VECTOR<T,2>& x1,const VECTOR<T,2>& x2,const VECTOR<T,2>& x3)
{
    bool is_degenerate;
    return (Intersects< EXACT_FLOAT<T> >(y,VECTOR<VECTOR<T,2>,3>(x1,x2,x3),&is_degenerate)||is_degenerate);
}
//#####################################################################
// Function Segment_Segment_Intersects
//#####################################################################
template<class T> bool POLYGONAL_TRIANGULATION<T>::
Segment_Segment_Intersects(const VECTOR<T,2>& x1,const VECTOR<T,2>& x2,const VECTOR<T,2>& y1,const VECTOR<T,2>& y2)
{
    bool is_degenerate;
    return (Intersects< EXACT_FLOAT<T> >(VECTOR<VECTOR<T,2>,2>(x1,x2),VECTOR<VECTOR<T,2>,2>(y1,y2),&is_degenerate)||is_degenerate);
}
//#####################################################################
// Function Triangulate_Nonconvex_Simple_Polygon
//#####################################################################
template<class T> template<class VECTORT2_ARRAY,class INT_ARRAY> int POLYGONAL_TRIANGULATION<T>::
Triangulate_Nonconvex_Simple_Polygon(const VECTORT2_ARRAY& coordinates,const INT_ARRAY& polygon_input,ARRAY< VECTOR<int,3> >& triangles,bool keep_degenerate_triangles)
{
    // assumption: the polygon is positively oriented consecutive vertices are distinct all segments are either pairwise non-intersecting or share 1 or 2 vertices
#ifndef NDEBUG
    assert(Polygon_Orientation(coordinates.Subset(polygon_input))>=0);
    for(int i=1;i<=polygon_input.m;++i){
        int inext=i%polygon_input.m+1,v=polygon_input(i),vnext=polygon_input(inext);
        assert(v!=vnext);assert(coordinates(v)!=coordinates(vnext));}
    for(int i=1;i<=polygon_input.m-2;++i) for(int j=i+2;j<=polygon_input.m;++j){
        int v1=polygon_input(i),v2=polygon_input(i+1),u1=polygon_input(j),u2=polygon_input(j%polygon_input.m+1);
        if(v1==u1||v1==u2||v2==u1||v2==u2) continue;
        const VECTOR<T,2> &x1=coordinates(v1),&x2=coordinates(v2),&y1=coordinates(u1),&y2=coordinates(u2);
        assert(!Segment_Segment_Intersects(x1,x2,y1,y2));}
#endif
    if(polygon_input.m<3) return 0;
    if(polygon_input.m==3){triangles.Append(VECTOR<int,3>(polygon_input(1),polygon_input(2),polygon_input(3)));return 1;}
    ARRAY<int> polygon(polygon_input);
    int num_triangles_added=0;
    HASHTABLE<VECTOR<int,3>,int> triple_to_orientation;
    VECTOR<int,3> triple;
    int infinite_loop_detector=0;
    for(int i=polygon.m;polygon.m>2;i=(i+polygon.m-2)%polygon.m+1){
        if(infinite_loop_detector>polygon.m){
            LOG::cerr<<std::endl;
            LOG::cerr<<"coordinates (input) = "<<std::endl;
            for(int j=1;j<=polygon_input.m;++j) LOG::cerr<<coordinates(polygon_input(j))<<std::endl;
            LOG::cerr<<"polygon_input = "<<polygon_input<<std::endl;
            for(int j=1;j<=polygon.m;++j) LOG::cerr<<coordinates(polygon(j))<<std::endl;
            LOG::cerr<<"polygon = "<<polygon<<std::endl;
            PHYSBAM_FATAL_ERROR("infinite loop");}
        infinite_loop_detector++;
        // attempt to remove the "ear" centered at i
        int iprev=(i+polygon.m-2)%polygon.m+1,inext=i%polygon.m+1,vprev=polygon(iprev),v=polygon(i),vnext=polygon(inext);
        const VECTOR<T,2> &xprev=coordinates(vprev),&x=coordinates(v),&xnext=coordinates(vnext);
        if(xprev==xnext) continue;
        triple.Set(vprev,v,vnext);
        // first check: triangle must be positively oriented
        int orientation;
        if(!triple_to_orientation.Get(triple,orientation)) triple_to_orientation.Insert(triple,orientation=Orientation(xprev,x,xnext));
        if(orientation<0) continue;
        if(orientation>0){ // second check: triangle must not intersect another edge
            bool intersects=false;
            for(int j=inext%polygon.m+1;j!=iprev&&!intersects;j=j%polygon.m+1){
                int u=polygon(j);
                if(u==vprev||u==vnext) continue;
                const VECTOR<T,2>& y=coordinates(u);
                if(Orientation(xprev,x,y)<0 || Orientation(y,x,xnext)<0) continue;
                // y is in the closed sector <xprev,x,xnext>
                if(u!=v&&(intersects=(Orientation(xprev,y,xnext)>=0))) break; // y is in the closed triangle <xprev,x,xnext> minus x
                int uprev=polygon((j+polygon.m-2)%polygon.m+1),unext=polygon(j%polygon.m+1);
                // this is necessary to catch segments passing through v
                intersects=(uprev==v||unext==v);}
            if(intersects){ // record the fact that this triangle is bad by recording its orientation as -1
                triple_to_orientation.Set(triple,-1);continue;}}
        // ear is removable, so remove it
        if(orientation>0||keep_degenerate_triangles){triangles.Append(triple);++num_triangles_added;}
        polygon.Remove_Index(i);i++;infinite_loop_detector=0;}
    return num_triangles_added;
}
//#####################################################################
// Function Triangulate_Nonconvex_Nonsimple_Polygon
//#####################################################################
template<class T> template<class VECTORT2_ARRAY,class INT_ARRAY_ARRAY> int POLYGONAL_TRIANGULATION<T>::
Triangulate_Nonconvex_Nonsimple_Polygon(const VECTORT2_ARRAY& coordinates,const INT_ARRAY_ARRAY& polygon_input,ARRAY< VECTOR<int,3> >& triangles,bool keep_degenerate_triangles)
{
    // assumption: outer loop is positively oriented, inner loops are negatively oriented, loops are simple and pairwise non-intersecting
    if(polygon_input.m==0) return 0;
    if(polygon_input.m==1) return Triangulate_Nonconvex_Simple_Polygon(coordinates,polygon_input(1),triangles,keep_degenerate_triangles);
    ARRAY< ARRAY<int> > polygon(polygon_input.m);
    for(int ell=1;ell<=polygon.m;++ell) polygon(ell)=polygon_input(ell);
    int infinite_loop_detector=0;
    for(int ell=polygon.m;polygon.m>1;ell=(ell+polygon.m-2)%polygon.m+1){
        if(ell==1) continue;
        if(++infinite_loop_detector>polygon.m) PHYSBAM_FATAL_ERROR();
        // attempt to find a non-intersecting segment between loop ell and loop 1
        ARRAY<int>& loop_1=polygon(1);
        const ARRAY<int>& loop_ell=polygon(ell);
        bool spliced=false;
        for(int i=1;i<=loop_1.m&&!spliced;++i){
            int v1=loop_1(i);
            const VECTOR<T,2>& x1=coordinates(v1);
            for(int j=1;j<=loop_ell.m&&!spliced;++j){
                int v2=loop_ell(j);
                const VECTOR<T,2>& x2=coordinates(v2);
                // trial segment is (x1,x2)
                bool intersects=false;
                // see if trial segment intersects with any other loop m
                int m=ell;
                do{
                    const ARRAY<int>& loop_m=polygon(m);
                    for(int k=1;k<=loop_m.m&&!intersects;++k){
                        int u1=loop_m(k),u2=loop_m(k%loop_m.m+1);assert(u1!=u2);
                        if(u1==v1||u1==v2||u2==v1||u2==v2) continue;
                        const VECTOR<T,2> &y1=coordinates(u1),&y2=coordinates(u2);
                        intersects=Segment_Segment_Intersects(x1,x2,y1,y2);}
                }while((m=m%polygon.m+1)!=ell&&!intersects);
                if(intersects) continue;
                // trial segment is good; splice loop_ell into loop_1
                ARRAY<int> temp_loop_1;
                temp_loop_1.Preallocate(loop_1.m+loop_ell.m+2);
                for(int k=1;k<=i;++k) temp_loop_1.Append(loop_1(k));
                for(int k=j;k<=loop_ell.m;++k) temp_loop_1.Append(loop_ell(k));
                for(int k=1;k<=j;++k) temp_loop_1.Append(loop_ell(k));
                for(int k=i;k<=loop_1.m;++k) temp_loop_1.Append(loop_1(k));
                loop_1.Exchange(temp_loop_1);
                spliced=true;}}
        if(spliced){ // remove loop ell
            polygon(ell).Exchange(polygon.Last());
            polygon.Remove_Index(polygon.m);ell++;infinite_loop_detector=0;}}
    assert(polygon.m==1);
    ARRAY<int>& polygon1=polygon(1);

    // unknot vertices, if necessary
    HASHTABLE<int,int> num_paths_through(polygon1.m);
    HASHTABLE<int,int> knotted_vertices;
    for(int i=1;i<=polygon1.m;++i){
        int v=polygon1(i);int& n=num_paths_through.Get_Or_Insert(v,0);
        if(++n>=3) knotted_vertices.Set(v,n);}
    for(HASHTABLE_ITERATOR<int,int> it(knotted_vertices);it.Valid();it.Next()){
        int v=it.Key();
        int n=it.Data();
        // mark where all the loops through v begin
        ARRAY<PAIR<int,int> > vnexts;
        vnexts.Preallocate(n);
        for(int i=1;i<=polygon1.m;++i) if(polygon1(i)==v){int inext=i%polygon1.m+1;vnexts.Append(Tuple(polygon1(inext),inext));}
        assert(vnexts.m==n);
        // (insertion) sort the loops in counter-clockwise order
        const VECTOR<T,2>& x=coordinates(v);
        for(int j=2;j<=n;++j){
            PAIR<int,int> vnextj_pair=vnexts(j);
            const VECTOR<T,2>& xj=coordinates(vnextj_pair.x);
            int k;
            if(Orientation(coordinates(vnexts(1).x),x,xj)>=0){
                for(k=2;k<j&&Orientation(coordinates(vnexts(k).x),x,xj)>0;++k);}
            else{
                for(k=j-1;k>=2&&Orientation(coordinates(vnexts(k).x),x,xj)<0;--k);
                if(k>1) ++k;}
            for(int m=j;m>k;--m) vnexts(m)=vnexts(m-1);
            vnexts(k)=vnextj_pair;}
        // rearrange the loops in polygon1
        ARRAY<int> new_polygon1;
        new_polygon1.Preallocate(polygon1.m);
        for(int j=1;j<=n;++j){
            new_polygon1.Append(v);
            int i=vnexts(j).y;
            for(int u;(u=polygon1(i))!=v;i=i%polygon1.m+1) new_polygon1.Append(u);}
        assert(new_polygon1.m==polygon1.m);
        polygon1.Exchange(new_polygon1);}

    return Triangulate_Nonconvex_Simple_Polygon(coordinates,polygon1,triangles,keep_degenerate_triangles);
}
//#####################################################################
template class POLYGONAL_TRIANGULATION<float>;
template int POLYGONAL_TRIANGULATION<float>::Triangulate_Nonconvex_Nonsimple_Polygon<ARRAY<VECTOR<float,2>,int>,ARRAY<ARRAY<int,int>,int> >(
    ARRAY<VECTOR<float,2>,int> const&,ARRAY<ARRAY<int,int>,int> const&,ARRAY<VECTOR<int,3>,int>&,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POLYGONAL_TRIANGULATION<double>;
template int POLYGONAL_TRIANGULATION<double>::Triangulate_Nonconvex_Nonsimple_Polygon<ARRAY<VECTOR<double,2>,int>,ARRAY<ARRAY<int,int>,int> >(
    ARRAY<VECTOR<double,2>,int> const&,ARRAY<ARRAY<int,int>,int> const&,ARRAY<VECTOR<int,3>,int>&,bool);
#endif
