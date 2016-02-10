//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_CAMERA
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_CAMERA__
#define __OPTIX_CAMERA__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>

namespace PhysBAM{

// Note: after each invokation of a function, camera stays consistent
template<class T>
class OPTIX_CAMERA:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    TV position; // camera position
    TV focal_point; // where the image plane is located
    TV look_vector; // points from the position to the focal point - normalized
    TV horizontal_vector; // points to the right on the omage plane - normalized
    TV vertical_vector; // point up in the image plane - normalized

    T proj_width,proj_height,proj_distance;
    unsigned int screen_width,screen_height;

    OPTIX_CAMERA(TV position_input=TV((T)0,(T)0.1,(T)1.0),TV focal_point_input=TV(0,0,0),TV horizontal_vector_input=TV((T)1.0,0,0),TV vertical_vector_input=TV(0,(T)1.0,0))
        :position(position_input),focal_point(focal_point_input),horizontal_vector(horizontal_vector_input),vertical_vector(vertical_vector_input)
    {
        look_vector=(focal_point-position).Normalized();
        horizontal_vector=TV::Cross_Product(look_vector,vertical_vector).Normalized();
        vertical_vector=TV::Cross_Product(horizontal_vector,look_vector).Normalized();
        proj_height=proj_width=(T)1.0; // angle of view is 45 degrees
        proj_distance=(T)1.0;
        screen_width=800;
        screen_height=800;
    }

    void Position_And_Aim_Camera(const TV& position_input,const TV& look_at_point,const TV& pseudo_up_vector)
    {
        T focal_distance=(focal_point-position)*look_vector;
        focal_point=position+focal_distance*look_vector;
        position=position_input;
        look_vector=(look_at_point-position).Normalized();
        horizontal_vector=TV::Cross_Product(look_vector,pseudo_up_vector).Normalized();
        vertical_vector=TV::Cross_Product(horizontal_vector,look_vector).Normalized();
    }

    void Focus_Camera_And_Give_Projection_Metrics(unsigned int width, unsigned int height,const T focal_distance,const T angle_of_view)
    {
        T aspect_ratio=width/(T)height;
        focal_point=position+focal_distance*look_vector;
        proj_width=(T)2*focal_distance*tan((T).5*angle_of_view);proj_height=width/aspect_ratio;
    }

    void rotateAroundRight(T d) 
    {
        MATRIX<T,3> Right_Rotation_Matrix = MATRIX<T,3>::Rotation_Matrix(horizontal_vector,d);
        T focal_distance=TV::Dot_Product((focal_point-position),look_vector);
        look_vector=Right_Rotation_Matrix*look_vector;
        position=focal_point-look_vector*focal_distance;

        vertical_vector=Right_Rotation_Matrix*vertical_vector;
        horizontal_vector=Right_Rotation_Matrix*horizontal_vector;
    }

    void rotateY(T d) 
    {
        MATRIX<T,3> OY_Rotation_Matrix=MATRIX<T,3>::Rotation_Matrix_Y_Axis(d);
        T focal_distance=TV::Dot_Product((focal_point-position),look_vector);
        vertical_vector=OY_Rotation_Matrix*vertical_vector;
        horizontal_vector=OY_Rotation_Matrix*horizontal_vector;
        look_vector=OY_Rotation_Matrix*look_vector;
        position=focal_point-look_vector*focal_distance;
    }

    void moveStraight(T d) 
    {
        position+=d*look_vector;
        if (TV::Dot_Product(focal_point-position,look_vector)<0){
            look_vector=(focal_point-position).Normalized();
            horizontal_vector*=-1;
        }
    }

    VECTOR<T,3> getWorldDeltaFromScreenDelta(int dx,int dy,TV V) 
    {
        V-=position;
        TV v=V;
        float k=proj_distance/TV::Dot_Product(v,look_vector);
        v*=k;
        float frac=V.Magnitude()/v.Magnitude();
        float proj_dx=(dx/(float)screen_width)*proj_width;
        float proj_dy=(dy/(float)screen_height)*proj_height;

        return frac*(proj_dx*horizontal_vector-proj_dy*vertical_vector);
    }

//#####################################################################
};
}
#endif
#endif
