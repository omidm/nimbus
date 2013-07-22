//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LINEAR_BLEND_SKINNING__
#define __LINEAR_BLEND_SKINNING__
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Motion/BODY_MOTION_SEQUENCE.h>

namespace PhysBAM{

template<class T>
class LINEAR_BLEND_SKINNING
{
    typedef VECTOR<T,3> TV;
public:
    SPARSE_MATRIX_FLAT_NXN<T> h;
    SPARSE_MATRIX_FLAT_NXN<T> surface_laplacian;
    ARRAY<VECTOR_ND<T> > weights;
    ARRAY<VECTOR_ND<T> > closest_points;
    ARRAY<bool> touched_particles,touched_elements;
    ARRAY<T> areas;
    ARRAY<int> row_lengths;
    T constant,motion_frame_rate;
    int start_particles,number_particles;

    TRIANGULATED_SURFACE<T> &surface;
    BODY_MOTION_SEQUENCE<T> &body_motion;

    LINEAR_BLEND_SKINNING(TRIANGULATED_SURFACE<T> &surface_input,BODY_MOTION_SEQUENCE<T> &body_motion_input,T motion_frame_rate_input,int number_particles_input=0,int start_particles_input=1)
        :constant(0.22),motion_frame_rate(motion_frame_rate_input),start_particles(start_particles_input),surface(surface_input),body_motion(body_motion_input)
    {if(number_particles_input) number_particles=number_particles_input; else number_particles=surface.particles.array_collection->Size()-start_particles+1;}

    int Find_Closest_Bone_And_Distance_Squared(T distance_squared,TV location,int frame)
    {
        int bone=0;
        for(int i=1;i<=body_motion.trajectories.m;i++){
            TV bone_length=TV(body_motion.trajectories(i)(frame).length,0,0);
            TV bone_vector=body_motion.trajectories(i)(frame).targeted_transform*bone_length;
            TV distance_vector=location-body_motion.trajectories(i)(frame).targeted_transform*TV(0,0,0);
            T new_distance_squared=(distance_vector-TV::Dot_Product(distance_vector,bone_vector)).Magnitude_Squared();
            if(bone==0||distance_squared>new_distance_squared){distance_squared=new_distance_squared;bone=i;}}
        return bone;
    }

    void Calculate_Weights_Helper(int element,int frame,int sanity_check)
    {
        assert(element>0);assert(element<=surface.mesh.elements.m);
        if(touched_elements(element)) return;
        touched_elements(element)=true;
        T distance_squared=0;int bone;
        int i,j,k;surface.mesh.elements(element).Get(i,j,k);
        
        if(!touched_particles(i-start_particles+1)){
            bone=Find_Closest_Bone_And_Distance_Squared(distance_squared,surface.particles.X(i),frame);
            closest_points(bone)(i-start_particles+1)=(T)1.;
            h.Set_Element(i-start_particles+1,i-start_particles+1,constant/distance_squared);
            touched_particles(i-start_particles+1)=true;}
        if(!touched_particles(j-start_particles+1)){
            bone=Find_Closest_Bone_And_Distance_Squared(distance_squared,surface.particles.X(j),frame);
            closest_points(bone)(j-start_particles+1)=(T)1.;
            h.Set_Element(j-start_particles+1,j-start_particles+1,constant/distance_squared);
            touched_particles(j-start_particles+1)=true;}
        if(!touched_particles(k-start_particles+1)){
            bone=Find_Closest_Bone_And_Distance_Squared(distance_squared,surface.particles.X(k),frame);
            closest_points(bone)(k-start_particles+1)=(T)1.;
            h.Set_Element(k-start_particles+1,k-start_particles+1,constant/distance_squared);
            touched_particles(k-start_particles+1)=true;}

        std::stringstream ss;
        if(element==1) ss<<"vert ("<<i<<","<<j<<","<<k<<")"<<std::endl;
        for(int index=1;index<=(*surface.mesh.adjacent_elements)(element).m;index++){
            int next_element=(*surface.mesh.adjacent_elements)(element)(index);
            if(element==1||next_element==1) ss<<"Found element: "<<element<<" next element: "<<next_element<<std::endl;
            if(touched_elements(next_element)) continue;
            if(element==1||next_element==1) ss<<"Dealing with"<<std::endl;
            int l,m,n,a,b,c,d;surface.mesh.elements(next_element).Get(l,m,n);
            if(i!=l&&i!=m&&i!=n){a=i;b=j;c=k;} else if(j!=l&&j!=m&&j!=n){a=j;b=i;c=k;} else{a=k;b=i;c=j;}
            if(l!=i&&l!=j&&l!=k) d=l; else if(m!=i&&m!=j&&m!=k) d=m; else d=n;
            if(element==1||next_element==1) ss<<"vert ("<<i<<","<<j<<","<<k<<") ("<<l<<","<<m<<","<<n<<") ("<<a<<","<<b<<","<<c<<","<<d<<")"<<std::endl;
            TV u1=(surface.particles.X(b)-surface.particles.X(a)),v1=(surface.particles.X(c)-surface.particles.X(a));
            TV u2=(surface.particles.X(b)-surface.particles.X(d)),v2=(surface.particles.X(c)-surface.particles.X(d));
            T area=TRIANGLE_3D<T>::Area(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k));
            area+=TRIANGLE_3D<T>::Area(surface.particles.X(l),surface.particles.X(m),surface.particles.X(n));
            T cot_theta1=TV::Dot_Product(u1,v1)/TV::Cross_Product(u1,v1).Magnitude(),cot_theta2=TV::Dot_Product(u2,v2)/TV::Cross_Product(u2,v2).Magnitude();
            areas(b-start_particles+1)+=area;areas(c-start_particles+1)+=area;
            if(sanity_check) row_lengths(b-start_particles+1)+=1;row_lengths(c-start_particles+1)+=1;
            if(element==1||next_element==1) ss<<"Adding to "<<b<<" and "<<c<<std::endl;
            surface_laplacian.Add_Element(b-start_particles+1,b-start_particles+1,cot_theta1+cot_theta2);
            surface_laplacian.Add_Element(c-start_particles+1,c-start_particles+1,cot_theta1+cot_theta2);
            surface_laplacian.Set_Symmetric_Elements(b-start_particles+1,c-start_particles+1,(T)3./area*(cot_theta1+cot_theta2));}
        //PHYSBAM_FATAL_ERROR();
        LOG::filecout(ss.str());
    }

    void Calculate_Weights(int frame=0,int sanity_check=0)
    {
        weights.Remove_All();weights.Resize(body_motion.trajectories.m);
        closest_points.Remove_All();closest_points.Resize(body_motion.trajectories.m);
        for(int i=1;i<=body_motion.trajectories.m;i++){closest_points(i).Resize(number_particles);weights(i).Resize(number_particles);}
        surface.mesh.Initialize_Neighbor_Nodes();surface.mesh.Initialize_Adjacent_Elements();surface.mesh.Initialize_Neighbor_Elements();
        row_lengths.Remove_All();row_lengths.Resize(number_particles,false,false);
        for(int i=1;i<=number_particles;i++) row_lengths(i)=(*surface.mesh.neighbor_nodes)(start_particles+i-1).m+1;
        surface_laplacian.Set_Row_Lengths(row_lengths);
        for(int i=1;i<=number_particles;i++) row_lengths(i)=1;
        h.Set_Row_Lengths(row_lengths);
 
        touched_elements.Resize(surface.mesh.elements.m,false,false);for(int i=1;i<=surface.mesh.elements.m;i++) touched_elements(i)=false;
        touched_particles.Resize(number_particles,false,false);for(int i=1;i<=number_particles;i++) touched_particles(i)=false;
        areas.Resize(number_particles,false,false);for(int i=1;i<=number_particles;i++) areas(i)=(T)0;

        for(int i=1;i<=number_particles;i++){
            surface_laplacian.Set_Element(i,i,0);
            for(int j=1;j<=(*surface.mesh.neighbor_nodes)(i).m;j++){surface_laplacian.Set_Symmetric_Elements(i,(*surface.mesh.neighbor_nodes)(i)(j),0);}}

        assert(surface.mesh.elements.m>1);
        for(int i=1;i<=surface.mesh.elements.m;i++) Calculate_Weights_Helper(i,frame+1,sanity_check);

        if(sanity_check){
            for(int i=1;i<=surface.mesh.elements.m;i++) if(!touched_elements(i)) PHYSBAM_FATAL_ERROR("Missed an element");
            for(int i=1;i<=number_particles;i++) if(!touched_particles(i)) {
                std::stringstream ss;
                ss<<"Particle "<<i<<std::endl;
                ARRAY<int> one_ring;
                ss<<"Elements are: (";
                for(int j=1;j<=surface.mesh.elements.m;j++) {
                    int a,b,c;surface.mesh.elements(j).Get(a,b,c);
                    if(a==i||b==i||c==i){
                        ss<<","<<j;
                        one_ring.Append_Unique(a);
                        one_ring.Append_Unique(b);
                        one_ring.Append_Unique(c);}}
                ss<<")"<<std::endl;
                ss<<"One_Ring is "<<one_ring<<std::endl;
                LOG::filecout(ss.str());
                PHYSBAM_FATAL_ERROR("Missed a particle");}
            for(int i=1;i<=number_particles;i++) if(row_lengths(i)>(*surface.mesh.neighbor_nodes)(start_particles+i-1).m+1){
                ARRAY<int> one_ring;
                std::stringstream ss;
                ss<<"Mismatch: set elements for row "<<i<<" was set "<<row_lengths(i)<<" times but should have "<<(*surface.mesh.neighbor_nodes)(start_particles+i-1).m+1<<" times."<<std::endl;
                ss<<"Elements are: (";
                for(int j=1;j<=surface.mesh.elements.m;j++) {
                    int a,b,c;surface.mesh.elements(j).Get(a,b,c);
                    if(a==i||b==i||c==i){
                        ss<<","<<j;
                        one_ring.Append_Unique(a);
                        one_ring.Append_Unique(b);
                        one_ring.Append_Unique(c);}}
                ss<<")"<<std::endl;
                ss<<"One_Ring is "<<one_ring<<std::endl;
                LOG::filecout(ss.str());
                PHYSBAM_FATAL_ERROR();}
        }
        
        for(int i=1;i<=number_particles;i++) surface_laplacian.Set_Element(i,i,surface_laplacian(i,i)*(T)6./areas(i));

        for(int i=1;i<=body_motion.trajectories.m;i++){
            VECTOR_ND<T> b(number_particles);
            surface_laplacian.Negate();
            for(int j=1;j<=number_particles;j++) surface_laplacian.Add_Element(j,j,h(j,j));
            h.Times(closest_points(i),b);
            surface_laplacian.Gauss_Seidel_Solve(weights(i),b);}
    }

    void Update_Particles(ARRAY_VIEW<TV> X,T time)
    {
        int frame=(int)(time*motion_frame_rate)+1;
        T alpha=time*motion_frame_rate-frame+1;
        for(int i=start_particles;i<=number_particles;i++){
            X(i)=TV();
            for(int j=1;j<=body_motion.trajectories.m;j++){
                FRAME<TV> transform=FRAME<TV>::Interpolation(body_motion.trajectories(j)(frame).targeted_transform,body_motion.trajectories(j)(frame+1).targeted_transform,alpha);
                X(i)+=transform.t*weights(j)(i);}}
    }

//#####################################################################
};
}
#endif
