//#####################################################################
// Copyright 2011-
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_CHIMERA_VOXELS
//#####################################################################
#ifndef __RENDERING_CHIMERA_VOXELS__
#define __RENDERING_CHIMERA_VOXELS__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/SMOOTH_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_VOXELS.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

template<class T>
class triple_y_greater_than
{
public:
    bool operator()(TRIPLE<int,T,T> a,TRIPLE<int,T,T> b)
    {return a.y>b.y;}
};

template<class T>
class RENDERING_CHIMERA_VOXELS:public RENDERING_UNIFORM_VOXELS<T>
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    typedef ARRAY<TRIPLE<int,T,T> > POINT_SAMPLES;
public:
    using RENDERING_VOXELS<T>::box;using RENDERING_VOXELS<T>::small_number;using RENDERING_VOXELS<T>::precompute_single_scattering;
    using RENDERING_UNIFORM_VOXELS<T>::data;using RENDERING_UNIFORM_VOXELS<T>::data_scale;using RENDERING_UNIFORM_VOXELS<T>::data_offset;using RENDERING_UNIFORM_VOXELS<T>::data_clamp_low_value;using RENDERING_UNIFORM_VOXELS<T>::data_lowest_value;using RENDERING_UNIFORM_VOXELS<T>::data_clamp_high_value;using RENDERING_UNIFORM_VOXELS<T>::data_highest_value;using RENDERING_UNIFORM_VOXELS<T>::number_of_smoothing_steps;using RENDERING_UNIFORM_VOXELS<T>::volumetric_step;using RENDERING_UNIFORM_VOXELS<T>::map_from_accessor_index_to_my_index;using RENDERING_UNIFORM_VOXELS<T>::default_voxel_source_interpolation;using RENDERING_UNIFORM_VOXELS<T>::default_voxel_light_interpolation;

    ARRAY<GRID<TV>*> grids;
    ARRAY<FRAME<TV>*> grid_frames;
    ARRAY<ARRAY<ARRAY<TV,VECTOR<int,3> >*> > precomputed_lights;
    ARRAY<ARRAY<ARRAY<bool,VECTOR<int,3> >*> > precomputed_light_valids;

    //we inherit this from uniform voxels: ARRAY<ARRAY<T,VECTOR<int,3> >*> data; // defined at each grid point
private:
    const ARRAY<ARRAY<int> >* side_neighbor_ranks;

    bool use_custom_interpolations;
    bool use_cubic_for_density;

    int number_of_ghost_cells;
    
    ARRAY<INTERPOLATION_UNIFORM<GRID<TV>,T>*> voxel_source_interpolations;
    ARRAY<INTERPOLATION_UNIFORM<GRID<TV>,TV>*> voxel_light_interpolations;
    ARRAY<RANGE<TV> > blend_inner_domains;
    ARRAY<RANGE<TV> > blend_outer_domains;
    

public:

    RENDERING_CHIMERA_VOXELS(const ARRAY<GRID<TV>*>& grids_input,const ARRAY<FRAME<TV>*>& grid_frames_input,const ARRAY<ARRAY<T,VECTOR<int,3> >*>& data_input,const T volumetric_step,const ARRAY<ARRAY<int> >* side_neighbor_ranks_input)
        :RENDERING_UNIFORM_VOXELS<T>(*grids_input(1),*data_input(1),volumetric_step),grids(grids_input),grid_frames(grid_frames_input),side_neighbor_ranks(side_neighbor_ranks_input),use_custom_interpolations(false),use_cubic_for_density(false),number_of_ghost_cells(3)
    {
        box.Change_Size(-0.5*grids(1)->dX.Max());
        for(int grid_index=2;grid_index<=grids.m;++grid_index) data.Append(data_input(grid_index)); // grid 1 is already appended in base class
        if(!use_custom_interpolations) for(int grid_index=1;grid_index<=grids.m;++grid_index)
            {voxel_source_interpolations.Append(&default_voxel_source_interpolation);voxel_light_interpolations.Append(&default_voxel_light_interpolation);}
        for(int grid_index=1;grid_index<=grids.m;++grid_index){
            blend_inner_domains.Append(Interpolation_Domain(grid_index));
            LOG::cerr<<"Domain of grid "<<grid_index<<" is "<<blend_inner_domains(blend_inner_domains.m)<<std::endl;
            blend_outer_domains.Append(Interpolation_Domain(grid_index,number_of_ghost_cells));}
    }

    RANGE<TV> Interpolation_Domain(int grid_index,int n_ghost_cells=0) const PHYSBAM_OVERRIDE
    {
        T thicken_cells=T(-.5)+n_ghost_cells;
        RANGE<TV> domain=grids(grid_index)->domain;
        TV dx=grids(grid_index)->dX;
        if(!side_neighbor_ranks) domain.Change_Size(thicken_cells*dx);
        else{
            for(int n=1;n<=GRID<TV>::number_of_neighbors_per_node;n++)
            if((*side_neighbor_ranks)(grid_index)(n)<0){
                TV_INT direction=GRID<TV>::Node_Neighbor(TV_INT(),n);
                int axis=direction.Dominant_Axis();
                if(direction(axis)==1) domain.max_corner(axis)+=thicken_cells*dx(axis);
                else domain.min_corner(axis)-=thicken_cells*dx(axis);}
        }
        return domain;
    }

    void Print_Parameters()
    {
        for(int grid_index=1;grid_index<=grids.m;++grid_index) std::cout<<"grid "<<grid_index<<" min_corner"<<grids(grid_index)->domain.min_corner<<" max_corner"<<grids(grid_index)->domain.max_corner<<std::endl;
        std::stringstream ss;
        for(int i=1;i<=1/*data.Size()*/;i++){
            ss<<"data_scale("<<i<<")="<<data_scale(i)<<" data_offset("<<i<<")="<<data_offset(i)<<std::endl;
            ss<<"data_clamp_low_value("<<i<<")="<<data_clamp_low_value(i)<<" data_lowest_value("<<i<<")="<<data_lowest_value(i)<<std::endl;
            ss<<"data_clamp_high_value("<<i<<")="<<data_clamp_high_value(i)<<" data_highest_value("<<i<<")="<<data_highest_value(i)<<std::endl;
            if(side_neighbor_ranks){
                ss<<"side_neighbor_ranks"<<std::endl;
                for(int i=1;i<=(*side_neighbor_ranks).m;++i)
                    ss<<i<<": "<<(*side_neighbor_ranks)(i)<<std::endl;
            }
            ss<<std::endl;
        }
        LOG::filecout(ss.str());
        //for(int i=127;i<=131;++i)for(int j=33;j<=36;++j)
        //    LOG::cerr<<"g2c["<<i<<", "<<j<<", "<<32<<"]"<<" = "<<(*data(2))(i,j,32)<<std::endl;
    }

    T Volumetric_Integration_Step(const RAY<TV> &ray,const T xi) const PHYSBAM_OVERRIDE
    {return xi*volumetric_step;}

    int Get_Point_Samples(const TV& location) const PHYSBAM_OVERRIDE
    {
        int enclosing_gi=0; T finest_enclosing_grid_cell_size=FLT_MAX;
        for(int grid_index=1;grid_index<=grids.m;++grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(location);
            if(blend_inner_domains(grid_index).Lazy_Inside(object_space_location) && grids(grid_index)->Cell_Size()<finest_enclosing_grid_cell_size){
                enclosing_gi=grid_index;finest_enclosing_grid_cell_size=grids(grid_index)->Cell_Size();
            }
        }
        /*if(enclosing_gi==3 || enclosing_gi==4 || enclosing_gi==9 || enclosing_gi==10) LOG::cerr<<" "<<enclosing_gi<<std::endl;*/
        return enclosing_gi;
    }
    POINT_SAMPLES* Get_Point_Samples_In_Ghost_Reqion(const int enclosing_gi, const TV& location) const PHYSBAM_OVERRIDE
    {
        POINT_SAMPLES *point_samples=new POINT_SAMPLES;
        for(int grid_index=1;grid_index<=grids.m;++grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(location);
            if(grid_index!=enclosing_gi && blend_outer_domains(grid_index).Lazy_Inside(object_space_location) && grids(grid_index)->Cell_Size()<grids(enclosing_gi)->Cell_Size()){
                TV dist_to_min= blend_inner_domains(grid_index).min_corner-object_space_location;
                TV dist_to_max=-blend_inner_domains(grid_index).max_corner+object_space_location;
                dist_to_min/=T(number_of_ghost_cells)*grids(grid_index)->dX;
                dist_to_max/=T(number_of_ghost_cells)*grids(grid_index)->dX;
                T min_dist=max((T)0,dist_to_min.Max(),dist_to_max.Max());
                point_samples->Append(TRIPLE<int,T,T>(grid_index,grids(grid_index)->Cell_Size(),1-min_dist));
                //LOG::cerr<<"Appending: ("<<grid_index<<", "<<grids(grid_index)->dX<<", "<<1-min_dist<<") os_p="<<object_space_location<<" os_ci="<<(object_space_location-(grids(grid_index)->domain.min_corner-(T)0.5*grids(grid_index)->dX))/grids(grid_index)->dX<<std::endl;
            }
        }
        return point_samples;
    }

    T Source_Term(const int source_term_index,const TV& location) const PHYSBAM_OVERRIDE
    {
        T value(0);
        if(int enclosing_gi=Get_Point_Samples(location)){
            TV object_space_location=grid_frames(enclosing_gi)->Inverse_Times(location);
            value=voxel_source_interpolations(enclosing_gi)->Clamped_To_Array(*grids(enclosing_gi),*data(enclosing_gi),object_space_location);
            if(use_cubic_for_density) value=value*value*value;
            POINT_SAMPLES* point_samples=Get_Point_Samples_In_Ghost_Reqion(enclosing_gi,location);
            if(point_samples->m){
                if(point_samples->m==1){
                    int gi=(*point_samples)(1).x;T alpha=(*point_samples)(1).z;
                    T finer_value=voxel_source_interpolations(gi)->Clamped_To_Array(*grids(gi),*data(gi),grid_frames(gi)->Inverse_Times(location));
                    if(use_cubic_for_density) finer_value=finer_value*finer_value*finer_value;
                    //if(finer_value!=0) LOG::cerr<<"gi="<<gi<<" value="<<value<<" finer_value="<<finer_value<<" alpha="<<alpha;
                    value=alpha*finer_value+(1-alpha)*value;
                    //if(finer_value!=0) LOG::cerr<<" afblend value="<<value<<std::endl;
                }
                else{
                    Sort<POINT_SAMPLES,triple_y_greater_than<T> >(*point_samples,triple_y_greater_than<T>());
                    for(int i=1;i<=point_samples->m;++i){
                        int gi=(*point_samples)(i).x;T alpha=(*point_samples)(i).z;
                        T finer_value=voxel_source_interpolations(gi)->Clamped_To_Array(*grids(gi),*data(gi),grid_frames(gi)->Inverse_Times(location));
                        if(use_cubic_for_density) finer_value=finer_value*finer_value*finer_value;
                        value=alpha*finer_value+(1-alpha)*value;
                    }
                }
            }
            delete point_samples;
        }
        if(data_clamp_low_value(1)) value=max(value,data_lowest_value(1));
        if(data_clamp_high_value(1)) value=min(value,data_highest_value(1));
        return value;
    }

    void Get_Node_Locations(ARRAY<TV>& locations) PHYSBAM_OVERRIDE
    {int n_locations=0;
    for(int grid_index=1;grid_index<=grids.m;++grid_index) {n_locations+=grids(grid_index)->counts.Product();std::cout<<"counts.Product()"<<grids(grid_index)->counts.Product()<<std::endl;}
    std::cout<<"n_locations="<<n_locations<<std::endl;
    locations.Resize(n_locations);
    map_from_accessor_index_to_my_index.Resize(locations.m);int index=1;
    for(int grid_index=1;grid_index<=grids.m;++grid_index)for(int i=1;i<=grids(grid_index)->counts.x;i++)for(int j=1;j<=grids(grid_index)->counts.y;j++)for(int ij=1;ij<=grids(grid_index)->counts.z;ij++){
        map_from_accessor_index_to_my_index(index)=VECTOR<int,3>(i,j,ij);
        locations(index)=(*grid_frames(grid_index))*grids(grid_index)->X(i,j,ij);index++;}} // (yuey) world space location is stored here

    bool Use_Precomputed_Light_Data(const TV& location,const int light_index) const PHYSBAM_OVERRIDE
    {if(!precompute_single_scattering)return false;
        T finest_gi=0, finest_enclosing_grid_cell_size=FLT_MAX;
        for(int grid_index=1;grid_index<=grids.m;++grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(location);
            if(grids(grid_index)->domain.Inside(object_space_location,0.5*grids(grid_index)->dX.Max()) && grids(grid_index)->Cell_Size()<finest_enclosing_grid_cell_size){
                finest_gi=grid_index;
                finest_enclosing_grid_cell_size=grids(grid_index)->Cell_Size();
            }
        }
        if(finest_gi){
            TV object_space_location=grid_frames(finest_gi)->Inverse_Times(location);
            VECTOR<int,3> index=INTERPOLATION_UNIFORM<GRID<TV>,TV>::Clamped_Index_End_Minus_One(*grids(finest_gi),*precomputed_lights(finest_gi)(light_index),object_space_location);
            int i=index.x,j=index.y,ij=index.z;
            ARRAY<bool,VECTOR<int,3> >& valid=*precomputed_light_valids(finest_gi)(light_index);
            bool i_j_ij=valid(i,j,ij),i_j_ij1=valid(i,j,ij+1),i_j1_ij=valid(i,j+1,ij),
            i_j1_ij1=valid(i,j+1,ij+1),i1_j_ij=valid(i+1,j,ij),i1_j_ij1=valid(i+1,j,ij+1),
            i1_j1_ij=valid(i+1,j+1,ij),i1_j1_ij1=valid(i+1,j+1,ij+1);
            return i_j_ij&&i_j_ij1&&i_j1_ij&&i_j1_ij1&&i1_j_ij&&i1_j_ij1&&i1_j1_ij&&i1_j1_ij1;}
        else return false;
    }
    
    void Set_Precomputed_Light_Data(const int location_index,const int light_index,const TV& light_value) PHYSBAM_OVERRIDE
    {int n_locations=0,grid_index=1;
    for(;grid_index<=grids.m;++grid_index){
        n_locations+=grids(grid_index)->counts.Product();
        if(location_index<=n_locations) break;}
    assert(grid_index<=grids.m);
    (*precomputed_lights(grid_index)(light_index))(map_from_accessor_index_to_my_index(location_index))=light_value;}

    void Set_Precomputed_Light_Valid(const int location_index,const int light_index,const bool value) PHYSBAM_OVERRIDE
    {int n_locations=0,grid_index=1;
    for(;grid_index<=grids.m;++grid_index){
        n_locations+=grids(grid_index)->counts.Product();
        if(location_index<=n_locations) break;}
    assert(grid_index<=grids.m);
    (*precomputed_light_valids(grid_index)(light_index))(map_from_accessor_index_to_my_index(location_index))=value;}

    TV Precomputed_Light_Data(const TV& location,const int light) const PHYSBAM_OVERRIDE
    {
        TV value=TV();
        T finest_gi=0, finest_enclosing_grid_cell_size=FLT_MAX;
        for(int grid_index=1;grid_index<=grids.m;++grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(location);
            if(grids(grid_index)->domain.Inside(object_space_location,0.5*grids(grid_index)->dX.Max()) && grids(grid_index)->Cell_Size()<finest_enclosing_grid_cell_size){
                finest_gi=grid_index;
                finest_enclosing_grid_cell_size=grids(grid_index)->Cell_Size();
            }
        }
        if(finest_gi){
            TV object_space_location=grid_frames(finest_gi)->Inverse_Times(location);
            if(use_custom_interpolations) value=voxel_light_interpolations(finest_gi)->Clamped_To_Array(*grids(finest_gi),*precomputed_lights(finest_gi)(light),object_space_location);
            else value=default_voxel_light_interpolation.Clamped_To_Array(*grids(finest_gi),*precomputed_lights(finest_gi)(light),object_space_location);}
        return value;
    }

    void Set_Custom_Source_And_Light_Interpolations(const ARRAY<INTERPOLATION_UNIFORM<GRID<TV>,T>*>& source_interpolations,const ARRAY<INTERPOLATION_UNIFORM<GRID<TV>,TV>*>& light_interpolations)
    {use_custom_interpolations=true;
    voxel_source_interpolations=source_interpolations;
    voxel_light_interpolations=light_interpolations;}

    void Set_Default_Source_And_Light_Interpolation()
    {use_custom_interpolations=false;}

    void Set_Cubic_Density(bool cubic=true)
    {use_cubic_for_density=cubic;}

protected:
    void Prepare_For_Precomputation(RENDER_WORLD<T>& world) PHYSBAM_OVERRIDE
    {precomputed_lights.Resize(grids.m);precomputed_light_valids.Resize(grids.m);std::cout<<"grids.m"<<grids.m<<std::endl;
    for(int grid_index=1;grid_index<=grids.m;++grid_index){
        GRID<TV>* local_grid=grids(grid_index);
        precomputed_lights(grid_index).Resize(world.Lights().m);precomputed_light_valids(grid_index).Resize(world.Lights().m);
        for(int i=1;i<=precomputed_lights(grid_index).m;i++)precomputed_lights(grid_index)(i)=new ARRAY<TV,VECTOR<int,3> >(1,local_grid->counts.x,1,local_grid->counts.y,1,local_grid->counts.z);
        for(int i=1;i<=precomputed_light_valids(grid_index).m;i++){precomputed_light_valids(grid_index)(i)=new ARRAY<bool,VECTOR<int,3> >(1,local_grid->counts.x,1,local_grid->counts.y,1,local_grid->counts.z);precomputed_light_valids(grid_index)(i)->Fill(false);}}}
    
    void Postprocess_Light_Field()
    {
        if(number_of_smoothing_steps) for(int grid_index=1;grid_index<=grids.m;++grid_index) for(int light=1;light<=precomputed_lights(grid_index).m;light++){
            std::stringstream ss;
            ss<<"Smoothing light "<<light<<" "<<number_of_smoothing_steps<<" steps"<<std::endl;
            LOG::filecout(ss.str());
            SMOOTH::Smooth<GRID<TV> >(*precomputed_lights(grid_index)(light),number_of_smoothing_steps,0);}
    }

};   
}
#endif

