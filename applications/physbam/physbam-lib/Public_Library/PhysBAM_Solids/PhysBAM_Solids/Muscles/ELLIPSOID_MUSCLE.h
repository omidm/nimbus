//#####################################################################
// Copyright 2006, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELLIPSOID_MUSCLE
//#####################################################################
#ifndef __ELLIPSOID_MUSCLE__
#define __ELLIPSOID_MUSCLE__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
namespace PhysBAM{

template<class TV>
class ELLIPSOID_MUSCLE:public MUSCLE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef MUSCLE<TV> BASE;
    using BASE::attachment_point_1;using BASE::attachment_point_2;using BASE::via_points;

    ARRAY<ELLIPSOID<T>*> ellipsoid_list;
    ARRAY<T> volume_list;
    T shorter_radii;
    T three_fourths_div_pi;

    ELLIPSOID_MUSCLE(const MUSCLE_FORCE_CURVE<T>& force_curve_input,T shorter_radii_input=(T)1):MUSCLE<TV>(force_curve_input),shorter_radii(shorter_radii_input),three_fourths_div_pi(3/(4*pi))
    {}

    virtual ~ELLIPSOID_MUSCLE()
    {
        ellipsoid_list.Delete_Pointers_And_Clean_Memory();
    }

    void Clean_Memory() PHYSBAM_OVERRIDE
    {
        ellipsoid_list.Delete_Pointers_And_Clean_Memory();
    }

    void Update_Positions_And_Radii()
    {
        if(!via_points.m)ellipsoid_list(1)->center=(attachment_point_1->Position()+attachment_point_2->Position())/2.0;
        else{
            ellipsoid_list(1)->center=(attachment_point_1->Position()+via_points(1)->Position())/2.0;
            for(int i=1;i<via_points.m;i++)ellipsoid_list(i+1)->center=(via_points(i)->Position()+via_points(i+1)->Position())/2.0;
            ellipsoid_list(ellipsoid_list.m)->center=(via_points(via_points.m)->Position()+attachment_point_2->Position())/2.0;}

        T radius_x,radius_yz;
        if(!via_points.m){
            radius_x=((attachment_point_1->Position()-attachment_point_2->Position()).Magnitude())/2.0;
            radius_yz=sqrt(three_fourths_div_pi*volume_list(1)*(1/radius_x));
            ellipsoid_list(1)->radii=DIAGONAL_MATRIX<T,3>(radius_x,radius_yz,radius_yz);}
        else{
            radius_x=((attachment_point_1->Position()-via_points(1)->Position()).Magnitude())/2.0;
            radius_yz=sqrt(three_fourths_div_pi*volume_list(1)*(1/radius_x));
            ellipsoid_list(1)->radii=DIAGONAL_MATRIX<T,3>(radius_x,radius_yz,radius_yz);
            for(int i=1;i<via_points.m;i++){
                radius_x=((via_points(i)->Position()-via_points(i+1)->Position()).Magnitude())/2.0;
                radius_yz=sqrt(three_fourths_div_pi*volume_list(i+1)*(1/radius_x));
                ellipsoid_list(i+1)->radii=DIAGONAL_MATRIX<T,3>(radius_x,radius_yz,radius_yz);}
            radius_x=((via_points(via_points.m)->Position()-attachment_point_2->Position()).Magnitude())/2.0;
            radius_yz=sqrt(three_fourths_div_pi*volume_list(via_points.m+1)*(1/radius_x));
            ellipsoid_list(via_points.m+1)->radii=DIAGONAL_MATRIX<T,3>(radius_x,radius_yz,radius_yz);}
    }

    ARRAY<ROTATION<TV> > Update_Orientations()
    {
        ARRAY<ROTATION<TV> > orientations;
        orientations.Resize(ellipsoid_list.m);
        VECTOR<T,3> rotation_vector(1,0,0);
        if(!via_points.m)orientations(1)=ROTATION<TV>::Rotation_Quaternion(rotation_vector,attachment_point_1->Position()-attachment_point_2->Position());
        else{
            orientations(1)=ROTATION<TV>::Rotation_Quaternion(rotation_vector,attachment_point_1->Position()-via_points(1)->Position());
            for(int i=1;i<via_points.m;i++)orientations(i+1)=ROTATION<TV>::Rotation_Quaternion(rotation_vector,via_points(i)->Position()-via_points(i+1)->Position());
            orientations(via_points.m+1)=ROTATION<TV>::Rotation_Quaternion(rotation_vector,via_points(via_points.m)->Position()-attachment_point_2->Position());}
        return orientations;
    }

    ARRAY<int> Find_Touching_Tets(TETRAHEDRALIZED_VOLUME<T>* tet_vol)
    {
        ARRAY<int> touching_particle_indices;
        VECTOR<T,3> curr_position;
        for (int i=1;i<=tet_vol->particles.array_collection->Size();i++){
            curr_position=tet_vol->particles.X(i);
            if(Inside_Ellipsoid_Test(curr_position))touching_particle_indices.Append(i);}

        ARRAY<bool> tet_indices;tet_indices.Resize(tet_vol->tetrahedron_list.m);
        for(int i=1;i<tet_indices.m;i++)tet_indices(i)=false;

        for(int i=1;i<=touching_particle_indices.m;i++)
            for(int j=1;j<=(tet_vol->mesh->incident_elements(i)).m;j++)tet_indices((tet_vol->mesh->incident_elements(i))(j))=true;

        ARRAY<int> touching_tet_indices;
        for(int i=1;i<tet_indices.m;i++)if(tet_indices(i))touching_tet_indices.Append(i);

        return touching_tet_indices;
    }

    bool Inside_Ellipsoid_Test(VECTOR<T,3>& position)
    {
        {std::stringstream ss;ss<<"POS 1: "<<position(1)<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"POS 2: "<<position(2)<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"POS 3: "<<position(3)<<std::endl;LOG::filecout(ss.str());}

        ARRAY<ROTATION<TV> > orientations=Update_Orientations();
        VECTOR<T,3> inv_rot_pos;
        for(int i=1;i<=ellipsoid_list.m;i++){
            inv_rot_pos=orientations(i).Inverse()*position;
            T result=sqr(inv_rot_pos(1)-(ellipsoid_list(i)->center)(1))/sqr((ellipsoid_list(i)->radius)(1))+
                sqr(inv_rot_pos(2)-(ellipsoid_list(i)->center)(2))/sqr((ellipsoid_list(i)->radius)(2))+
                sqr(inv_rot_pos(3)-(ellipsoid_list(i)->center)(3))/sqr((ellipsoid_list(i)->radius)(3));
            if(result<=1)return true;}
        return false;
    }

    void Set_Shorter_Radii(T shorter_radii_input)
    {
        shorter_radii=shorter_radii_input;
    }

    void Setup_Ellipsoid_List(int ellipsoid_list_size)
    {
        ellipsoid_list.Resize(ellipsoid_list_size);
        for(int i=1;i<=ellipsoid_list.m;i++)ellipsoid_list(i)=new ELLIPSOID<T>();

        // calculate volumes in volume list
        volume_list.Resize(ellipsoid_list_size);
        T radius_x,four_thirds_pi_times_shorter_radii_squared=(4/3)*pi*sqr(shorter_radii);

        if (!via_points.m){
            radius_x=((attachment_point_1->Position()-attachment_point_2->Position()).Magnitude())/2.0;
            volume_list(1)=four_thirds_pi_times_shorter_radii_squared*radius_x;}
        else{
            radius_x=((attachment_point_1->Position()-via_points(1)->Position()).Magnitude())/2.0;
            volume_list(1)=four_thirds_pi_times_shorter_radii_squared*radius_x;
            for(int i=1;i<via_points.m;i++){
                radius_x=((via_points(i)->Position()-via_points(i+1)->Position()).Magnitude())/2.0;
                volume_list(i+1)=four_thirds_pi_times_shorter_radii_squared*radius_x;}

            radius_x=((via_points(via_points.m)->Position()-attachment_point_2->Position()).Magnitude())/2.0;
            volume_list(via_points.m+1)=four_thirds_pi_times_shorter_radii_squared*radius_x;}
    }

    /*
    ELLIPSOID<T>* Read_Ellipsoid(TYPED_ISTREAM& input_stream)
    {
        VECTOR<T,3> center,radius;
        Read_Binary(input_stream,center,radius);
        ELLIPSOID<T>* e=new ELLIPSOID<T>(center,radius);
        LOG::cout<<"CENTER: "<<e->center<<std::endl;
        LOG::cout<<"RADIUS: "<<e->radius<<std::endl;
        return e;
    }

    void Write_Ellipsoid(TYPED_OSTREAM& output_stream,const ELLIPSOID<T>* e) const
    {
        Write_Binary(output_stream,e->center,e->radius);
    }
    */

    void Read(TYPED_ISTREAM& input_stream,RIGID_BODY_COLLECTION<TV>& rigid_body_collection) PHYSBAM_OVERRIDE
    {
        MUSCLE<TV>::Read(input_stream,rigid_body_collection);

        Read_Binary(input_stream,shorter_radii);
        Set_Shorter_Radii(shorter_radii);

        int ellipsoid_list_size;
        Read_Binary(input_stream,ellipsoid_list_size);
        Setup_Ellipsoid_List(ellipsoid_list_size);

        // THIS IS THE INFORMATION THAT WAS READ IN
        {std::stringstream ss;ss<<"THE FOLLOWING INFO WAS READ IN\n";LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"Shorter_Radii: "<<shorter_radii<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"ELLIPSOID LIST SIZE: "<<ellipsoid_list.m<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"FIRST ELLIPSOID: "<<ellipsoid_list(1)->center<<" "<<ellipsoid_list(1)->radii<<std::endl;LOG::filecout(ss.str());}
    }

    void Write(TYPED_OSTREAM& output_stream) const PHYSBAM_OVERRIDE
    {
        MUSCLE<TV>::Write(output_stream);

        Write_Binary(output_stream,shorter_radii);
        Write_Binary(output_stream,ellipsoid_list.m);

        // THIS IS THE INFORMATION THAT WAS WRITTEN OUT
        {std::stringstream ss;ss<<"THE FOLLOWING INFO WAS WRITTEN OUT\n";LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"Shorter_Radii: "<<shorter_radii<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"ELLIPSOID LIST SIZE: "<<ellipsoid_list.m<<std::endl;LOG::filecout(ss.str());}
    }
};
}
#endif
