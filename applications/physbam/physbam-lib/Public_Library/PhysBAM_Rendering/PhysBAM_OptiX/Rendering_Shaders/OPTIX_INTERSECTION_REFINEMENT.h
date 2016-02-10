/*
 * Copyright (c) 2008 - 2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and proprietary
 * rights in and to this software, related documentation and any modifications thereto.
 * Any use, reproduction, disclosure or distribution of this software and related
 * documentation without an express license agreement from NVIDIA Corporation is strictly
 * prohibited.
 *
 * TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, THIS SOFTWARE IS PROVIDED *AS IS*
 * AND NVIDIA AND ITS SUPPLIERS DISCLAIM ALL WARRANTIES, EITHER EXPRESS OR IMPLIED,
 * INCLUDING, BUT NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE.  IN NO EVENT SHALL NVIDIA OR ITS SUPPLIERS BE LIABLE FOR ANY
 * SPECIAL, INCIDENTAL, INDIRECT, OR CONSEQUENTIAL DAMAGES WHATSOEVER (INCLUDING, WITHOUT
 * LIMITATION, DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS INTERRUPTION, LOSS OF
 * BUSINESS INFORMATION, OR ANY OTHER PECUNIARY LOSS) ARISING OUT OF THE USE OF OR
 * INABILITY TO USE THIS SOFTWARE, EVEN IF NVIDIA HAS BEEN ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGES
 */

#ifdef __CUDACC__

#include <optix.h>
#include <optixu/optixu_math_namespace.h>


// Plane intersection -- used for refining triangle hit points.  Note
// that this skips zero denom check (for rays perpindicular to plane normal)
// since we know that the ray intersects the plane.
__device__ __inline__ float intersectPlane( const optix::float3& origin,
                                            const optix::float3& direction,
                                            const optix::float3& normal,
                                            const optix::float3& point )
{
  // Skipping checks for non-zero denominator since we know that ray intersects this plane
  return -( optix::dot( normal, origin-point ) ) / optix::dot( normal, direction );

}

// Offset the hit point using integer arithmetic
__device__ __inline__ optix::float3 offset( const optix::float3& hit_point, const optix::float3& normal )
{
  using namespace optix;

  const float epsilon = 1.0e-4f;
  const float offset  = 4096.0f*2.0f;

  float3 offset_point = hit_point;
  if( (__float_as_int( hit_point.x )&0x7fffffff)  < __float_as_int( epsilon ) ) {
    offset_point.x += epsilon * normal.x;
  } else {
    offset_point.x = __int_as_float( __float_as_int( offset_point.x ) + int(copysign( offset, hit_point.x )*normal.x) );
  }

  if( (__float_as_int( hit_point.y )&0x7fffffff) < __float_as_int( epsilon ) ) {
    offset_point.y += epsilon * normal.y;
  } else {
    offset_point.y = __int_as_float( __float_as_int( offset_point.y ) + int(copysign( offset, hit_point.y )*normal.y) );
  }

  if( (__float_as_int( hit_point.z )&0x7fffffff)  < __float_as_int( epsilon ) ) {
    offset_point.z += epsilon * normal.z;
  } else {
    offset_point.z = __int_as_float( __float_as_int( offset_point.z ) + int(copysign( offset, hit_point.z )*normal.z) );
  }

  return offset_point;
}

// Refine the hit point to be more accurate and offset it for reflection and
// refraction ray starting points.
__device__ __inline__ void refine_and_offset_hitpoint( const optix::float3& original_hit_point, const optix::float3& direction,
                                                       const optix::float3& normal, const optix::float3& p,
                                                       optix::float3& back_hit_point,
                                                       optix::float3& front_hit_point )
{
  using namespace optix;

  // Refine hit point
  float  refined_t          = intersectPlane( original_hit_point, direction, normal, p );
  float3 refined_hit_point  = original_hit_point + refined_t*direction;

  // Offset hit point
  if( dot( direction, normal ) > 0.0f ) {
    back_hit_point  = offset( refined_hit_point,  normal );
    front_hit_point = offset( refined_hit_point, -normal );
  } else {
    back_hit_point  = offset( refined_hit_point, -normal );
    front_hit_point = offset( refined_hit_point,  normal );
  }
}
#endif
