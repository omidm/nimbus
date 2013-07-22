//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <optix_world.h>
#include <optixu/optixu_math_namespace.h>
#include "OPTIX_COMMONSTRUCTS.h"
#include "OPTIX_HELPERS.h"
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtDeclareVariable(float,front_box_hit_t,attribute front_box_hit_t,);
rtDeclareVariable(float,back_box_hit_t,attribute back_box_hit_t,);
rtDeclareVariable(optix::Ray,ray,rtCurrentRay,);
rtDeclareVariable(PerRayData_radiance,prd,rtPayload,);
rtDeclareVariable(PerRayData_shadow,prd_shadow,rtPayload,);
rtDeclareVariable(float3,low_corner,,);
rtDeclareVariable(float3,up_corner,,);
rtDeclareVariable(float3,one_over_box_size,,);
rtDeclareVariable(float,t_step,,)=0.016f;
rtDeclareVariable(float,t_step_2,,)=0.016f;
rtDeclareVariable(rtObject,top_opaque_object,,);

rtTextureSampler<float,3,cudaReadModeElementType> soot_texture;
rtTextureSampler<float,3,cudaReadModeElementType> temperature_texture;
rtTextureSampler<float4,1,cudaReadModeElementType> fire_color_texture;

__device__ float4 temperature_color_table(float normalized_temperature)
{
	return tex1D(fire_color_texture,normalized_temperature);
}

RT_PROGRAM void closest_hit_radiance_fire()
{
	PerRayData_radiance bg_prd;
	bg_prd.importance=1.0f;
	bg_prd.depth=0;
	optix::Ray bg_ray_radiance(ray.origin,ray.direction,RAY_TYPE_RADIANCE,RAY_T_MIN,RAY_T_MAX);
	rtTrace(top_opaque_object,bg_ray_radiance,bg_prd);
	float3 color=bg_prd.result;

	float3 start_point=ray.origin+back_box_hit_t*ray.direction;
	int n=int((back_box_hit_t-front_box_hit_t)/t_step);
	float absorption=4.0f*t_step;
	for(int i=0;i<n;i++){
		float3 pos=start_point-(float)i*t_step*ray.direction;
		float3 tex_coord=(pos-low_corner)*one_over_box_size;
		float sampled_density=tex3D(soot_texture,tex_coord.x,tex_coord.y,tex_coord.z);
		float normalized_sampled_temperature=tex3D(temperature_texture,tex_coord.x,tex_coord.y,tex_coord.z)/3000.0f;
		float attenuation=expf(-sampled_density*absorption);
		float3 emitting_color=make_float3(temperature_color_table(normalized_sampled_temperature));
		color=color*attenuation+0.3f*emitting_color*sampled_density;
	}

	prd.result=color;
}

RT_PROGRAM void closest_hit_radiance_smoke()
{
    float3 hit_point=ray.origin+front_box_hit_t*ray.direction;
    int n=int((back_box_hit_t-front_box_hit_t)/t_step);
    float color=0;
	float shadow=0;
    float alpha=0;
    float coef_density_to_alpha=0.8f;
	float coef_density_to_shadow=0.0f;
	float coef_density_to_color=1.0f;
	bool flag_compute_bg_color=true;
    for(int i=0;i<n;i++){
        float3 pos=hit_point+(float)i*t_step*ray.direction;
        float3 tex_coord=(pos-low_corner)*one_over_box_size;
        float sampled_density=tex3D(soot_texture,tex_coord.x,tex_coord.y,tex_coord.z);
		float sampled_alpha=1.0f-expf(-coef_density_to_alpha*sampled_density);
		float sampled_shadow=coef_density_to_shadow*sampled_density;
		float sampled_color=coef_density_to_color*sampled_density;	////not associtated with opacity yet
		alpha+=(1.0f-alpha)*sampled_alpha;
		if(alpha>1.0f){flag_compute_bg_color=false;break;}
		shadow=sampled_shadow+(1.0f-sampled_alpha)*shadow;
		if(shadow>1.0f){flag_compute_bg_color=false;break;}
		color+=(1.0f-alpha)*sampled_color*sampled_alpha*(1.0f-shadow);
    }

	float3 bg_color=make_float3(0,0,0);
	if(flag_compute_bg_color){
		PerRayData_radiance bg_prd;
		bg_prd.importance=1.0f;
		bg_prd.depth=0;
		optix::Ray bg_ray_radiance(ray.origin,ray.direction,RAY_TYPE_RADIANCE,RAY_T_MIN,RAY_T_MAX);
		rtTrace(top_opaque_object,bg_ray_radiance,bg_prd);
		bg_color=bg_prd.result;
	}

    prd.result=make_float3(color,color,color)+(1.0f-alpha)*bg_color;
}

RT_PROGRAM void closest_hit_shadow_smoke()
{
	float3 start_point=ray.origin+front_box_hit_t*ray.direction;
	int n=int((back_box_hit_t-front_box_hit_t)/t_step_2);
	float attenuation=1.0f;
	float coef_density_to_transpancy=0.2f;

	for(int i=0;i<n;i++){
		float3 pos=start_point+(float)i*t_step*ray.direction;
		float3 tex_coord=(pos-low_corner)*one_over_box_size;
		float sampled_density=tex3D(soot_texture,tex_coord.x,tex_coord.y,tex_coord.z);
		float sampled_transparency=expf(-coef_density_to_transpancy*sampled_density);	////more accurate
		attenuation*=sampled_transparency;
	}

	prd_shadow.attenuation*=attenuation;
}