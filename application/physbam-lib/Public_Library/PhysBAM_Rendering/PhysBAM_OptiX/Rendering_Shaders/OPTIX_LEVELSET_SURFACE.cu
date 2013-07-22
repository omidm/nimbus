//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <optix_world.h>
#include "OPTIX_HELPERS.h"
using namespace optix;

rtDeclareVariable(float, scene_epsilon, , );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float3, phi_tex_delta, , );

rtTextureSampler<float, 3, cudaReadModeElementType> phi_tex;

// rtDeclareVariable(float3, back_hit_point, attribute back_hit_point, );
// rtDeclareVariable(float3, front_hit_point, attribute front_hit_point, );

rtDeclareVariable(float3, geometric_normal, attribute geometric_normal, );
rtDeclareVariable(float3, shading_normal, attribute shading_normal, );
__device__ bool isInsideDomain(float3 point) {
    float3 up_corner = make_float3(1.f, 1.f, 1.f);
    float3 low_corner = make_float3(0.f, 0.f, 0.f);

    return !(point.x < low_corner.x || point.x > up_corner.x ||
        point.y < low_corner.y || point.y > up_corner.y ||
        point.z < low_corner.z || point.z > up_corner.z);
}

__device__ float3 operator*(float3 a, uint3 b) {
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

__device__ float getIntensity(float3 point) {
    float3 up_corner = make_float3(1.f, 1.f, 1.f);
    float3 low_corner = make_float3(0.f, 0.f, 0.f);

    if (!isInsideDomain(point)) {
        return -1.f;
    }
    float3 tex_coord = (point - low_corner) / (up_corner - low_corner);
    // tex_coord = tex_coord * 2 - make_float3(1.f, 1.f, 1.f);
    // return tex_coord.x*tex_coord.x + tex_coord.y*tex_coord.y + tex_coord.z*tex_coord.z - 1;
    return tex3D(phi_tex, tex_coord.x, tex_coord.y, tex_coord.z);
}

__device__ float3 boxnormal(float t) {
  float3 up_corner = make_float3(1.f, 1.f, 1.f);
  float3 low_corner = make_float3(0.f, 0.f, 0.f);

  float3 t0 = (low_corner - ray.origin)/ray.direction;
  float3 t1 = (up_corner - ray.origin)/ray.direction;
  float3 neg = make_float3(t==t0.x?1:0, t==t0.y?1:0, t==t0.z?1:0);
  float3 pos = make_float3(t==t1.x?1:0, t==t1.y?1:0, t==t1.z?1:0);
  return pos-neg;
}

__device__ float intersect_box(optix::Ray ray, float &nearS, float &farS) {
    float3 up_corner = make_float3(1.f, 1.f, 1.f);
    float3 low_corner = make_float3(0.f, 0.f, 0.f);

    float3 t1 = (up_corner - ray.origin) / ray.direction;
    float3 t2 = (low_corner - ray.origin) / ray.direction;

    float3 T1 = fminf(t1, t2);
    float3 T2 = fmaxf(t1, t2);

    float Tnear = fmaxf(T1);
    float Tfar = fminf(T2);

    if (Tnear < Tfar && Tfar > 0) {
        if (Tnear < 0) {
            Tnear = scene_epsilon * 10;
        }
        /*
        shading_normal = geometric_normal = boxnormal(Tnear);
        if(rtPotentialIntersection(Tnear)) {
            shading_normal = geometric_normal = boxnormal(Tnear);
            front_hit_point = ray.origin + Tnear * ray.direction + geometric_normal * scene_epsilon;
            back_hit_point = front_hit_point - 2 * scene_epsilon * geometric_normal;
            rtReportIntersection(0);
        }
        */
        nearS = Tnear;
        farS = Tfar;
        return Tnear;
    }
    return -1;
}

RT_PROGRAM void intersect(int) {
    float d, d_far;
    intersect_box(ray, d, d_far);
    float step = fminf(phi_tex_delta);
    if (d > 0) {
        float3 hit_point = ray.origin + (d + scene_epsilon) * ray.direction;
        bool t = !isInsideDomain(ray.origin) || getIntensity(ray.origin) > 0;
        float current_distance;

        while (d < d_far) {
            current_distance = getIntensity(hit_point);
            if ((t && current_distance < 0) || (!t && current_distance > 0)) {
                if(rtPotentialIntersection(d)) {
                    shading_normal = geometric_normal = normalize(make_float3((getIntensity(hit_point + make_float3(phi_tex_delta.x, 0.f, 0.f)) - getIntensity(hit_point - make_float3(phi_tex_delta.x, 0.f, 0.f))) / (2 * phi_tex_delta.x),
                                                                              (getIntensity(hit_point + make_float3(0.f, phi_tex_delta.y, 0.f)) - getIntensity(hit_point - make_float3(0.f, phi_tex_delta.y, 0.f))) / (2 * phi_tex_delta.y),
                                                                              (getIntensity(hit_point + make_float3(0.f, 0.f, phi_tex_delta.z)) - getIntensity(hit_point - make_float3(0.f, 0.f, phi_tex_delta.z))) / (2 * phi_tex_delta.z)));
                    // shading_normal = geometric_normal = make_float3(0.f, 1.f, 0.f);
                    /*
                    front_hit_point = hit_point + scene_epsilon * geometric_normal;
                    back_hit_point = hit_point - scene_epsilon * geometric_normal;
                    */
                    // back_hit_point = front_hit_point = hit_point;
                    rtReportIntersection(0);
                }
                d -= step;
                step *= 0.1;
                if (step / 2 < 1e-4)
                    break;
            }
            d += step;
            hit_point = ray.origin + d * ray.direction;
        }
    }
    // shading_normal = geometric_normal = front_hit_point = back_hit_point = make_float3(0.f, 0.f, 0.f);
}

RT_PROGRAM void box_bounds (int, float result[6]) {
  optix::Aabb* aabb = (optix::Aabb*)result;

  aabb->m_min = make_float3(0, 0, 0);
  aabb->m_max = make_float3(1, 1, 1);
}


