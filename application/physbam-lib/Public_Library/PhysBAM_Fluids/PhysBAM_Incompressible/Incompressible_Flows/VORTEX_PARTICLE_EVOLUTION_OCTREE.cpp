#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Computations/VORTICITY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EULER_STEP_PARTICLES.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/DELETE_POINTS.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_OCTREE.h>
using namespace PhysBAM;
//#####################################################################
// Function Set_Radius
//#####################################################################
template<class T> inline T VORTEX_PARTICLE_EVOLUTION_OCTREE<T>::
Gaussian_Kernel(const T distance_squared)
{
    const T normalization=1/(T)pow(2*pi,1.5);
    return one_over_radius_cubed*normalization*exp(-(T).5*one_over_radius_squared*distance_squared);
}
//#####################################################################
// Function Set_Radius
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_OCTREE<T>::
Set_Radius(const T radius_input)
{
    radius=radius_input;radius_squared=sqr(radius);one_over_radius_squared=(T)1/sqr(radius);one_over_radius_cubed=(T)1/cube(radius);
}
//#####################################################################
// Function Compute_Body_Force
//#####################################################################
template<class T> struct COMPUTE_BODY_FORCE_HELPER{VORTICITY_PARTICLES<VECTOR<T,3> >* particles;int p;T radius_squared;RANGE<VECTOR<T,3> > box;VECTOR<T,3> missing_vorticity;OPERATION_HASH<>* hash;const ARRAY<VECTOR<T,3> >* V_ghost; 
    ARRAY<VECTOR<T,3> >* force;ARRAY<VECTOR<T,3> >* node_locations;T dt;T time;T force_scaling;VORTEX_PARTICLE_EVOLUTION_OCTREE<T>* evolution;OCTREE_GRID<T>* grid;};
template<class T> static void
Compute_Body_Force_Helper(OCTREE_CELL<T>* cell,COMPUTE_BODY_FORCE_HELPER<T>* helper)
{
    VECTOR<T,3> center=cell->Center(),dx_over_two=cell->DX()*T(0.5),lower=center-dx_over_two,upper=center+dx_over_two;
    if(!RANGE<VECTOR<T,3> >(lower,upper).Intersection(helper->box,(T)1e-5))return; // cull box
    
    if(cell->Has_Children()) for(int k=0;k<8;k++) Compute_Body_Force_Helper(cell->Child(k),helper);
    else for(int k=0;k<8;k++){
        int node_index=cell->Node(k);
        if(!helper->hash->Is_Marked_Current(node_index)){
            VECTOR<T,3> X_minus_Xp=helper->grid->Node_Location(node_index)-helper->particles->X(helper->p);T distance_squared=X_minus_Xp.Magnitude_Squared();
            if(distance_squared<=helper->radius_squared && distance_squared > 1e-5) (*helper->force)(node_index)-=(T)(helper->force_scaling*helper->evolution->Gaussian_Kernel(distance_squared)/sqrt(distance_squared))*VECTOR<T,3>::Cross_Product(X_minus_Xp,helper->missing_vorticity);
            helper->hash->Mark(node_index);}}
}
template<class T> void VORTEX_PARTICLE_EVOLUTION_OCTREE<T>::
Compute_Body_Force(const ARRAY<VECTOR<T,3> >& V_ghost,ARRAY<VECTOR<T,3> >& force,const T dt,const T time)
{
    // TODO: Andy, look at this code like we discussed. -Frank
    grid_vorticity.Resize(grid.number_of_nodes);VORTICITY_DYADIC<T>::Vorticity(grid,V_ghost,grid_vorticity);OPERATION_HASH<> hash(grid.number_of_nodes);
    VECTOR<T,3> radius_vector(radius,radius,radius);
    if(apply_individual_particle_forces){
        for(int p=1;p<=vorticity_particles.array_collection->Size();p++){
            hash.Next_Operation();
            VECTOR<T,3> local_grid_vorticity=LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,grid_vorticity,vorticity_particles.X(p));
            VECTOR<T,3> missing_vorticity=vorticity_particles.vorticity(p)-local_grid_vorticity;
            VECTOR<T,3> sign_check=missing_vorticity*vorticity_particles.vorticity(p);
            for(int a=1;a<=3;a++) if(sign_check[a]<0) missing_vorticity[a]=0;
            VECTOR<T,3> lower_vector=vorticity_particles.X(p)-radius_vector,upper_vector=vorticity_particles.X(p)+radius_vector;
            VECTOR<int,3> lower_index=grid.uniform_grid.Clamp_To_Cell(lower_vector),upper_index=grid.uniform_grid.Clamp_To_Cell(upper_vector);
            COMPUTE_BODY_FORCE_HELPER<T> helper;
            helper.p=p;helper.radius_squared=radius_squared;helper.box=RANGE<VECTOR<T,3> >(lower_vector,upper_vector);helper.missing_vorticity=missing_vorticity;helper.hash=&hash;helper.V_ghost=&V_ghost;helper.force=&force;
            helper.grid=&grid;helper.dt=dt;helper.time=time;helper.particles=&vorticity_particles;helper.force_scaling=force_scaling;helper.evolution=this;
            for(int i=lower_index.x;i<=upper_index.x;i++) for(int j=lower_index.y;j<=upper_index.y;j++) for(int ij=lower_index.z;ij<=upper_index.z;ij++) Compute_Body_Force_Helper<T>(grid.cells(i,j,ij),&helper);}}
    else PHYSBAM_NOT_IMPLEMENTED(); // don't have scattered interpolation in octree
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_OCTREE<T>::
Euler_Step(const ARRAY<VECTOR<T,3> >& V_ghost,const T dt,const T time)
{
    LOG::Time("Advancing vorticity particles");
    // vortex stretching/tilting term  - omega dot grad V
    T minimum_cell_size=grid.Minimum_Edge_Length();T one_over_two_dx=1/(2*minimum_cell_size);VECTOR<T,3> x_offset(minimum_cell_size,0,0),y_offset(0,minimum_cell_size,0),z_offset(0,0,minimum_cell_size);
    for(int p=1;p<=vorticity_particles.array_collection->Size();p++){
        VECTOR<T,3> X=vorticity_particles.X(p);
        VECTOR<T,3> new_vorticity=dt*MATRIX<T,3>(one_over_two_dx*(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+x_offset)-LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-x_offset)),
            one_over_two_dx*(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+y_offset)-LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-y_offset)),
            one_over_two_dx*(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+z_offset)-LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-z_offset)))*vorticity_particles.vorticity(p);
        if(renormalize_vorticity_after_stretching_tilting) vorticity_particles.vorticity(p)=vorticity_particles.vorticity(p).Magnitude()*new_vorticity.Normalized();
        else vorticity_particles.vorticity(p)=new_vorticity;}
    // advance vortex particle position
    EULER_STEP_PARTICLES<OCTREE_GRID<T> >::Euler_Step_Node(vorticity_particles.X,grid,V_ghost,dt);
    
    // delete particles outside grid
    int deleted_count=POINT_CLOUD_COMPUTATIONS::Delete_Points_Outside_Range(vorticity_particles,grid.Domain());
    std::stringstream ss;ss<<"Deleted "<<deleted_count<<" particles"<<std::endl;LOG::filecout(ss.str());
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_OCTREE<T>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
   LOG::Time("Writing Vortex Specific Data");
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/vorticity_particles",output_directory.c_str(),frame),vorticity_particles);
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/grid_vorticity",output_directory.c_str(),frame),grid_vorticity);
   LOG::Stop_Time();
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_OCTREE<T>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& input_directory,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/vorticity_particles",input_directory.c_str(),frame),vorticity_particles);
}
//#####################################################################
template class VORTEX_PARTICLE_EVOLUTION_OCTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTEX_PARTICLE_EVOLUTION_OCTREE<double>;
#endif
#endif
