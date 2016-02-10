//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <optix_world.h>
#include "OPTIX_HELPERS.h"
#include "OPTIX_INTERSECTION_REFINEMENT.h"
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtBuffer<float3> vertex_buffer;
rtBuffer<uint3> vindex_buffer;

rtBuffer<float3> normal_buffer;
// rtBuffer<uint3> nindex_buffer;

rtDeclareVariable(float3, geometric_normal, attribute geometric_normal, );
rtDeclareVariable(float3, shading_normal, attribute shading_normal, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float3, front_hit_point, attribute front_hit_point, );
rtDeclareVariable(float3, back_hit_point, attribute back_hit_point, );
rtDeclareVariable(float, scene_epsilon, , );
rtDeclareVariable(int, interpolate_normal, , );
rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );

RT_PROGRAM void intersect(int primIdx) {
  uint3 v_idx = vindex_buffer[primIdx];

  float3 p0 = vertex_buffer[ v_idx.x ];
  float3 p1 = vertex_buffer[ v_idx.y ];
  float3 p2 = vertex_buffer[ v_idx.z ];

  // Intersect ray with triangle
  float3 n;
  float  t, beta, gamma;
  if( intersect_triangle( ray, p0, p1, p2, n, t, beta, gamma ) ) {

    if(  rtPotentialIntersection( t ) ) {

      uint3 n_idx = v_idx;//nindex_buffer[ primIdx ];

      if ( normal_buffer.size() == 0 || n_idx.x < 0 || n_idx.y < 0 || n_idx.z < 0 ) {
        shading_normal = normalize( n );
      } else {
        float3 n0 = normal_buffer[ n_idx.x ];
        float3 n1 = normal_buffer[ n_idx.y ];
        float3 n2 = normal_buffer[ n_idx.z ];
        if (interpolate_normal) {
            shading_normal = normalize( n1*beta + n2*gamma + n0*(1.0f-beta-gamma) );
        } else {
            shading_normal = normalize( n );
        }
      }
      /* shading_normal = */geometric_normal = normalize( n );

      back_hit_point = front_hit_point = ray.origin + t * ray.direction;
      /*
      if (dot(ray.direction, shading_normal) > 0.f) {
        back_hit_point = front_hit_point - shading_normal * scene_epsilon;
        front_hit_point = front_hit_point + shading_normal * scene_epsilon;
      } else {
        back_hit_point = front_hit_point + shading_normal * scene_epsilon;
        front_hit_point = front_hit_point - shading_normal * scene_epsilon;
      }*/

      // refine_and_offset_hitpoint( ray.origin + t*ray.direction, ray.direction, geometric_normal, p0, back_hit_point, front_hit_point );
      // back_hit_point = front_hit_point = ray.origin + t * ray.direction;// + shading_normal * scene_epsilon;

      /*
      geometric_normal = shading_normal = normalize(n);
      int3 t_idx = tindex_buffer[ primIdx ];
      if ( texcoord_buffer.size() == 0 || t_idx.x < 0 || t_idx.y < 0 || t_idx.z < 0 ) {
        texcoord = make_float3( 0.0f, 0.0f, 0.0f );
      } else {
        float2 t0 = texcoord_buffer[ t_idx.x ];
        float2 t1 = texcoord_buffer[ t_idx.y ];
        float2 t2 = texcoord_buffer[ t_idx.z ];
        texcoord = make_float3( t1*beta + t2*gamma + t0*(1.0f-beta-gamma) );
      }

      rtReportIntersection(material_buffer[primIdx]);
      */
      rtReportIntersection(0);
    }
  }
}

RT_PROGRAM void box_bounds (int primIdx, float result[6]) {
  const uint3 v_idx = vindex_buffer[primIdx];

  const float3 v0   = vertex_buffer[ v_idx.x ];
  const float3 v1   = vertex_buffer[ v_idx.y ];
  const float3 v2   = vertex_buffer[ v_idx.z ];

  const float  area = length(cross(v1-v0, v2-v0));
  optix::Aabb* aabb = (optix::Aabb*)result;

  if(area > 0.0f && !isinf(area)) {
    aabb->m_min = fminf( fminf( v0, v1), v2 );
    aabb->m_max = fmaxf( fmaxf( v0, v1), v2 );
  } else {
    aabb->invalidate();
  }
}

