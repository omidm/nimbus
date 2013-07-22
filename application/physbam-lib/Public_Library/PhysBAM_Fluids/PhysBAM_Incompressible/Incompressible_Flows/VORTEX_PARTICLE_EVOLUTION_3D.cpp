//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EULER_STEP_PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Set_Radius
//#####################################################################
template<class T> inline T VORTEX_PARTICLE_EVOLUTION_3D<T>::
Gaussian_Kernel(const T distance_squared)
{
//    const T normalization=1/(T)pow(2*pi,1.5);
//    return one_over_radius_cubed*normalization*exp(-(T).5*one_over_radius_squared*distance_squared);
    return exp(-distance_squared*4); // returns 1 when distance=0 and almost 0 when distance=1
}
//#####################################################################
// Function Set_Radius
//#####################################################################
/*template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Set_Radius(const T radius_input)
{
    radius=radius_input;radius_squared=sqr(radius);one_over_radius_squared=(T)1/sqr(radius);one_over_radius_cubed=(T)1/cube(radius);scattered_interpolation.Set_Radius_Of_Influence(radius);
}*/
//#####################################################################
// Function Compute_Body_Force
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Compute_Body_Force(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& force,const T dt,const T time)
{
    ARRAY<T,VECTOR<int,3> > grid_vorticity_magnitude(grid.Domain_Indices(2),false);
    VORTICITY_UNIFORM<TV>::Vorticity(grid,FACE_LOOKUP_UNIFORM<GRID<TV> >(face_velocities_ghost),grid_vorticity,grid_vorticity_magnitude);

    if(apply_individual_particle_forces){
        // compute missing vorticity per particle
        VORTICITY_PARTICLES<VECTOR<T,3> > missing_vorticity_particles;missing_vorticity_particles.array_collection->Add_Elements(vorticity_particles.array_collection->Size());
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> > vorticity_interpolation;

        for(int p=1;p<=vorticity_particles.array_collection->Size();p++){
            // find missing vorticity
            VECTOR<T,3> local_grid_vorticity=vorticity_interpolation.Clamped_To_Array(grid,grid_vorticity,vorticity_particles.X(p));
            VECTOR<T,3> missing_vorticity=vorticity_particles.vorticity(p)-local_grid_vorticity;
            VECTOR<T,3> sign_check=missing_vorticity*vorticity_particles.vorticity(p);
            for(int a=1;a<=3;a++) if(sign_check[a]<0) missing_vorticity[a]=0;
            missing_vorticity_particles.X(p)=vorticity_particles.X(p);
            missing_vorticity_particles.vorticity(p)=missing_vorticity;
            missing_vorticity_particles.radius(p)=vorticity_particles.radius(p);}
        
        if(mpi_grid){
            T max_radius = (T)0;
            if(vorticity_particles.array_collection->Size()>=1) max_radius=ARRAYS_COMPUTATIONS::Max(vorticity_particles.radius);
            Exchange_Boundary_Particles_Flat(*mpi_grid,missing_vorticity_particles,max_radius);}

        T small_number=(T)1e-4*grid.min_dX;
        for(int p=1;p<=missing_vorticity_particles.array_collection->Size();p++){
            T radius=missing_vorticity_particles.radius(p);
            VECTOR<int,3> box_radius((int)(radius/grid.dX.x)+1,(int)(radius/grid.dX.y)+1,(int)(radius/grid.dX.z)+1);
            VECTOR<int,3> index=grid.Clamped_Index(missing_vorticity_particles.X(p));
            VECTOR<int,3> lower=clamp_min(index-box_radius,VECTOR<int,3>(1,1,1)),upper=clamp_max(index+box_radius,grid.counts);
            for(int i=lower.x;i<=upper.x;i++) for(int j=lower.y;j<=upper.y;j++) for(int ij=lower.z;ij<=upper.z;ij++){
                VECTOR<T,3> X_minus_Xp=grid.X(i,j,ij)-missing_vorticity_particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
                if(distance_squared>small_number&&distance_squared<=sqr(radius)) {
//                    force(i,j,ij)-=(T)(force_scaling*Gaussian_Kernel(distance_squared)/sqrt(distance_squared))*VECTOR<T,3>::Cross_Product(X_minus_Xp,missing_vorticity_particles.vorticity(p));}}}
                    T distance=sqrt(distance_squared);
                    force(i,j,ij)-=(T)(force_scaling*Gaussian_Kernel(sqr(distance/radius))/distance)*VECTOR<T,3>::Cross_Product(X_minus_Xp,missing_vorticity_particles.vorticity(p));}}}}
    else{
        if(mpi_grid) PHYSBAM_NOT_IMPLEMENTED(); // this has not been mpi'd yet
        ARRAY<T,VECTOR<int,3> > grid_vorticity_particles_magnitude(grid.Domain_Indices(2),false);
    
        // form grid vorticity from vortex particles
        scattered_interpolation.Transfer_To_Grid(vorticity_particles.X,vorticity_particles.vorticity,grid,grid_vorticity_particles);
        if(remove_grid_vorticity_from_particle_vorticity) for(int i=0;i<=grid.counts.x+1;i++) for(int j=0;j<=grid.counts.y+1;j++) for(int ij=0;ij<=grid.counts.z+1;ij++) 
            grid_vorticity_particles(i,j,ij)-=grid_vorticity(i,j,ij);
    
        // find vorticity magnitudes
        for(int i=grid_vorticity.domain.min_corner.x;i<=grid_vorticity.domain.max_corner.x;i++)for(int j=grid_vorticity.domain.min_corner.y;j<=grid_vorticity.domain.max_corner.y;j++)for(int ij=grid_vorticity.domain.min_corner.z;ij<=grid_vorticity.domain.max_corner.z;ij++)
            grid_vorticity_magnitude(i,j,ij)=grid_vorticity(i,j,ij).Magnitude();
        for(int i=0;i<=grid.counts.x+1;i++) for(int j=0;j<=grid.counts.y+1;j++) for(int ij=0;ij<=grid.counts.z+1;ij++)
            grid_vorticity_particles_magnitude(i,j,ij)=grid_vorticity_particles(i,j,ij).Magnitude();
    
        // compute confinement force
        T one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y),one_over_two_dz=1/(2*grid.dX.z);
        for(int i=force.domain.min_corner.x;i<=force.domain.max_corner.x;i++) for(int j=force.domain.min_corner.y;j<=force.domain.max_corner.y;j++) for(int ij=force.domain.min_corner.z;ij<=force.domain.max_corner.z;ij++){
            VECTOR<T,3> vortex_normal_vector((grid_vorticity_magnitude(i+1,j,ij)-grid_vorticity_magnitude(i-1,j,ij))*one_over_two_dx,
                                              (grid_vorticity_magnitude(i,j+1,ij)-grid_vorticity_magnitude(i,j-1,ij))*one_over_two_dy,
                                              (grid_vorticity_magnitude(i,j,ij+1)-grid_vorticity_magnitude(i,j,ij-1))*one_over_two_dz);
            VECTOR<T,3> particle_vortex_normal_vector((grid_vorticity_particles_magnitude(i+1,j,ij)-grid_vorticity_particles_magnitude(i-1,j,ij))*one_over_two_dx,
                                                       (grid_vorticity_particles_magnitude(i,j+1,ij)-grid_vorticity_particles_magnitude(i,j-1,ij))*one_over_two_dy,
                                                       (grid_vorticity_particles_magnitude(i,j,ij+1)-grid_vorticity_particles_magnitude(i,j,ij-1))*one_over_two_dz);
            vortex_normal_vector.Normalize();particle_vortex_normal_vector.Normalize();
            force(i,j,ij)=grid_confinement_parameter*VECTOR<T,3>::Cross_Product(vortex_normal_vector,grid_vorticity(i,j,ij))
                        +particle_confinement_parameter*VECTOR<T,3>::Cross_Product(particle_vortex_normal_vector,grid_vorticity_particles(i,j,ij));}}
}
//#####################################################################
// Function Compute_Body_Force
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Compute_Body_Force(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& force,const T dt,const T time)
{
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > cell_force(grid.Domain_Indices(1),false);
    Compute_Body_Force(face_velocities_ghost,cell_force,dt,time);
    for(T_FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        force.Component(axis)(iterator.Face_Index())+=(T).5*(cell_force(iterator.First_Cell_Index())[axis]+cell_force(iterator.Second_Cell_Index())[axis]);}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Euler_Step(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T dt,const T time)
{
    LOG::Time("Advancing vorticity particles");
    
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > two_times_V_ghost(grid.Domain_Indices(2));
    for(int i=-1;i<=grid.counts.x+2;i++) for(int j=-1;j<=grid.counts.y+2;j++) for(int ij=-1;ij<=grid.counts.z+2;ij++){
        two_times_V_ghost(i,j,ij).x=face_velocities_ghost.Component(1)(i,j,ij)+face_velocities_ghost.Component(1)(i+1,j,ij);
        two_times_V_ghost(i,j,ij).y=face_velocities_ghost.Component(2)(i,j,ij)+face_velocities_ghost.Component(2)(i,j+1,ij);
        two_times_V_ghost(i,j,ij).z=face_velocities_ghost.Component(3)(i,j,ij)+face_velocities_ghost.Component(3)(i,j,ij+1);}
    
    // vortex stretching/tilting term  - omega dot grad V
    ARRAY<MATRIX<T,3> ,VECTOR<int,3> > VX(grid.Domain_Indices(1),false);LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,MATRIX<T,3> > VX_interpolation;
    T one_over_four_dx=1/(4*grid.dX.x),one_over_four_dy=1/(4*grid.dX.y),one_over_four_dz=1/(4*grid.dX.z);
    for(int i=0;i<=grid.counts.x+1;i++) for(int j=0;j<=grid.counts.y+1;j++) for(int ij=0;ij<=grid.counts.z+1;ij++)
        VX(i,j,ij)=MATRIX<T,3>(one_over_four_dx*(two_times_V_ghost(i+1,j,ij)-two_times_V_ghost(i-1,j,ij)),
                                 one_over_four_dy*(two_times_V_ghost(i,j+1,ij)-two_times_V_ghost(i,j-1,ij)),
                                 one_over_four_dz*(two_times_V_ghost(i,j,ij+1)-two_times_V_ghost(i,j,ij-1)));

    if(renormalize_vorticity_after_stretching_tilting) for(int p=1;p<=vorticity_particles.array_collection->Size();p++){
        T old_magnitude=vorticity_particles.vorticity(p).Magnitude();
        vorticity_particles.vorticity(p)+=dt*VX_interpolation.Clamped_To_Array(grid,VX,vorticity_particles.X(p))*vorticity_particles.vorticity(p);
        vorticity_particles.vorticity(p).Normalize();vorticity_particles.vorticity(p)*=old_magnitude;}
    else for(int p=1;p<=vorticity_particles.array_collection->Size();p++){
        vorticity_particles.vorticity(p)+=dt*VX_interpolation.Clamped_To_Array(grid,VX,vorticity_particles.X(p))*vorticity_particles.vorticity(p);}
    // advance vortex particle position
    EULER_STEP_PARTICLES<GRID<TV> >::Euler_Step_Face(vorticity_particles.X,grid,face_velocities_ghost,dt);
    if(mpi_grid) Exchange_Boundary_Particles_Flat(*mpi_grid,vorticity_particles);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
   LOG::Time("Writing Vortex Specific Data");
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/vorticity_particles",output_directory.c_str(),frame),vorticity_particles);
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/grid_vorticity",output_directory.c_str(),frame),grid_vorticity);
   FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/grid_vorticity_particles",output_directory.c_str(),frame),grid_vorticity_particles);
   LOG::Stop_Time();
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_3D<T>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& input_directory,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/vorticity_particles",input_directory.c_str(),frame),vorticity_particles);
}
//#####################################################################
template class VORTEX_PARTICLE_EVOLUTION_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTEX_PARTICLE_EVOLUTION_3D<double>;
#endif
