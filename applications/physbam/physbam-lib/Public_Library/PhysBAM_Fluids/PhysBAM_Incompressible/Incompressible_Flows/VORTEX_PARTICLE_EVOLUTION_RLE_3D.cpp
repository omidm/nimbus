#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE_Computations/VORTICITY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE_HELPER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/DELETE_POINTS.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_RLE_3D.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_RLE_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
using namespace PhysBAM;
//#####################################################################
// Function Set_Radius
//#####################################################################
template<class T> inline T VORTEX_PARTICLE_EVOLUTION_RLE_3D<T>::
Gaussian_Kernel(const T distance_squared)
{
    static const T normalization=1/(T)pow(2*pi,1.5);
    return one_over_radius_cubed*normalization*exp(-(T).5*one_over_radius_squared*distance_squared);
}
//#####################################################################
// Function Set_Radius
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_RLE_3D<T>::
Set_Radius(const T radius_input)
{
    radius=radius_input;radius_squared=sqr(radius);one_over_radius_squared=(T)1/sqr(radius);one_over_radius_cubed=(T)1/cube(radius);
}
//#####################################################################
// Function Compute_Body_Force
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_RLE_3D<T>::
Compute_Body_Force(const ARRAY<T>& V_ghost,ARRAY<T>& force,const T dt,const T time)
{
    VORTICITY_RLE<TV>::Vorticity(grid,V_ghost,grid_vorticity);

/*
    ARRAY<T,FACE_INDEX<3> > V_uniform;Transfer_Faces_To_Uniform(grid,V_ghost,V_uniform);
    ARRAY<T_SPIN,VECTOR<int,3> > vorticity_uniform(grid.uniform_grid,grid.number_of_ghost_cells);
    INCOMPRESSIBLE_3D<T>::Vorticity(grid.uniform_grid,V_uniform,vorticity_uniform);
    ARRAY<T_SPIN,VECTOR<int,3> > vorticity_bad;Transfer_Cells_To_Uniform(grid,grid_vorticity,vorticity_bad);
    typedef typename GRID<TV>::CELL_ITERATOR UNIFORM_CELL_ITERATOR;
    for(UNIFORM_CELL_ITERATOR iterator(grid.uniform_grid,3);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(vorticity_bad(cell)==vorticity_uniform(cell) || vorticity_bad(cell)==TV()) continue;
        std::stringstream ss;ss<<"  cell "<<cell<<", rle "<<vorticity_bad(cell)<<", uniform "<<vorticity_uniform(cell)<<std::endl;LOG::filecout(ss.str());}
*/

    // compute missing vorticity per particle
    VORTICITY_PARTICLES<TV> missing_vorticity_particles;missing_vorticity_particles.array_collection->Add_Elements(vorticity_particles.array_collection->Size());
    if(force_scaling_from_age) missing_vorticity_particles.array_collection->template Add_Array<T>(ATTRIBUTE_ID_AGE);
    LINEAR_INTERPOLATION_RLE<T_GRID,T_SPIN> vorticity_interpolation;
    for(int p=1;p<=vorticity_particles.array_collection->Size();p++){
        // find missing vorticity
        T_BLOCK block(grid,vorticity_particles.X(p));if(!block) continue;
        VECTOR<T,3> local_grid_vorticity=vorticity_interpolation.From_Block_Cell(block,grid_vorticity,vorticity_particles.X(p));
        VECTOR<T,3> missing_vorticity=vorticity_particles.vorticity(p)-local_grid_vorticity;
        VECTOR<T,3> sign_check=missing_vorticity*vorticity_particles.vorticity(p);
        for(int a=1;a<=3;a++)if(sign_check[a]<=0) missing_vorticity[a]=0;
        missing_vorticity_particles.array_collection->Copy_Element(*vorticity_particles.array_collection,p,p);
        missing_vorticity_particles.vorticity(p)=missing_vorticity;}
        
    if(grid.mpi_grid) Exchange_Boundary_Particles(*grid.mpi_grid,missing_vorticity_particles,radius);

    T small_number=(T)1e-4*grid.Minimum_Edge_Length();
    ARRAY_VIEW<T>* missing_vorticity_particle_age=missing_vorticity_particles.array_collection->template Get_Array<T>(ATTRIBUTE_ID_AGE);
    if(force_scaling_from_age && !missing_vorticity_particle_age) PHYSBAM_FATAL_ERROR();
    for(int p=1;p<=missing_vorticity_particles.array_collection->Size();p++){
        // compute force
        RANGE<TV> box=RANGE<TV>(missing_vorticity_particles.X(p)).Thickened(radius);
        TV_INT min_index=grid.uniform_grid.Clamp_To_Cell(box.Minimum_Corner(),grid.number_of_ghost_cells);
        TV_INT max_index=grid.uniform_grid.Clamp_To_Cell(box.Maximum_Corner(),grid.number_of_ghost_cells);
        RANGE<VECTOR<int,1> > vertical_box(min_index.y,max_index.y);
        T_BOX_HORIZONTAL_INT horizontal_box=RANGE<TV_INT>(min_index,max_index).Get_Horizontal_Box();
        // TODO: this shouldn't loop over entire columns
        for(CELL_ITERATOR cell(grid,horizontal_box);cell;cell++) if(cell.Short() && vertical_box.Lazy_Inside(VECTOR<int,1>(cell.j))){
            TV X_minus_Xp=cell.X()-missing_vorticity_particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
            if(distance_squared>small_number&&distance_squared<=radius_squared){
                T particle_force_scaling=force_scaling_from_age?(*force_scaling_from_age)((*missing_vorticity_particle_age)(p)):force_scaling;
                TV force_cell=-particle_force_scaling*Gaussian_Kernel(distance_squared)/sqrt(distance_squared)*TV::Cross_Product(X_minus_Xp,missing_vorticity_particles.vorticity(p));
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    T half_force=(T).5*force_cell[axis];
                    force(cell.First_Face_Index(axis))+=half_force;force(cell.Second_Face_Index(axis))+=half_force;}}}}
    grid.Put_Ghost_Faces((T)0,force);
}
template<class T_GRID,class TV> static void Transfer_Cells_To_Uniform(const T_GRID& grid,const ARRAY<TV>& u,ARRAY<TV,VECTOR<int,3> >& u_uniform)
{
    u_uniform.Resize(grid.uniform_grid.Cell_Indices(grid.number_of_ghost_cells));
    for(typename T_GRID::CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++)if(cell.Short() && u_uniform.domain.min_corner.y<=cell.j && cell.j<=u_uniform.domain.max_corner.y)
        u_uniform(cell.I())=u(cell.Cell());
}
template<class T,class T_GRID> struct Transfer_Faces_To_Uniform_Helper{template<class T_FACE> static void Apply(const T_GRID& grid,const ARRAY<T>& u,ARRAY<typename T_GRID::SCALAR,FACE_INDEX<T_GRID::dimension> >& u_uniform)
{
    ARRAY<T,VECTOR<int,3> >& u_component=u_uniform.Component(T_FACE::Axis());
    for(T_FACE face(grid,grid.number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short() && u_component.domain.min_corner.y<=face.j() && face.j()<=u_component.domain.max_corner.y){
        u_component(face.cell2.I())=u(face.Face());}
}};
template<class T,class T_GRID> static void Transfer_Faces_To_Uniform(const T_GRID& grid,const ARRAY<T>& u,ARRAY<typename T_GRID::SCALAR,FACE_INDEX<T_GRID::dimension> >& u_uniform)
{
    u_uniform.Resize(grid.uniform_grid,grid.number_of_ghost_cells);
    T_GRID::template Face_Loop<Transfer_Faces_To_Uniform_Helper<T,T_GRID> >(grid,u,u_uniform);
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_RLE_3D<T>::
Euler_Step(const ARRAY<T>& V_ghost,const T dt,const T time)
{
    LOG::Time("Advancing vorticity particles");
    
    ARRAY<TV> V_ghost_cell(grid.number_of_cells);LINEAR_INTERPOLATION_RLE_HELPER<T_GRID>::Interpolate_From_Faces_To_Short_Cells(grid,V_ghost,V_ghost_cell);
    
    // vortex stretching/tilting term  - omega dot grad V
    ARRAY<MATRIX<T,3> > VX(grid.number_of_cells);
    TV one_over_two_DX=(T).5*grid.uniform_grid.one_over_dX;
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& short_cell_neighbors=grid.Short_Cell_Neighbors();
    for(int c=1;c<=grid.number_of_cells;c++){
        int neighbors[T_GRID::number_of_neighbors_per_cell];bool has_all_neighbors=true;
        for(int n=1;n<=T_GRID::number_of_neighbors_per_cell;n++){
            neighbors[n-1]=short_cell_neighbors(c)(n);
            if(!neighbors[n-1]){has_all_neighbors=false;break;}}
        if(!has_all_neighbors) continue;
        for(int axis=1;axis<=T_GRID::dimension;axis++)
            VX(c).Column(axis)=one_over_two_DX[axis]*(V_ghost_cell(neighbors[2*axis-1])-V_ghost_cell(neighbors[2*axis-2]));}

    LINEAR_INTERPOLATION_RLE<T_GRID,T> V_interpolation;
    LINEAR_INTERPOLATION_RLE<T_GRID,MATRIX<T,3> > VX_interpolation;
    for(int p=vorticity_particles.array_collection->Size();p>=1;p--){
        T_BLOCK block(grid,vorticity_particles.X(p));if(!block){vorticity_particles.array_collection->Delete_Element(p);continue;}
        // update vorticity
        T_SPIN new_vorticity=vorticity_particles.vorticity(p)+dt*VX_interpolation.From_Block_Cell(block,VX,vorticity_particles.X(p))*vorticity_particles.vorticity(p);
        if(renormalize_vorticity_after_stretching_tilting){new_vorticity.Normalize();new_vorticity*=vorticity_particles.vorticity(p).Magnitude();}
        vorticity_particles.vorticity(p)=new_vorticity;
        // advance position
        vorticity_particles.X(p)+=dt*V_interpolation.From_Block_Face(block,V_ghost,vorticity_particles.X(p));}
    if(grid.mpi_grid) Exchange_Boundary_Particles(*grid.mpi_grid,vorticity_particles);

    // delete particles outside grid
    int deleted_count=POINT_CLOUD_COMPUTATIONS::Delete_Points_Outside_Range(vorticity_particles,grid.Domain());
    std::stringstream ss;ss<<"deleted "<<deleted_count<<" particles"<<std::endl;LOG::filecout(ss.str());
    LOG::Stop_Time();
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_RLE_3D<T>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& input_directory,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/vorticity_particles",input_directory.c_str(),frame),vorticity_particles);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void VORTEX_PARTICLE_EVOLUTION_RLE_3D<T>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
   std::string f=FILE_UTILITIES::Number_To_String(frame);
   FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/vorticity_particles",vorticity_particles);
   FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid_vorticity",grid_vorticity);
}
//#####################################################################
template class VORTEX_PARTICLE_EVOLUTION_RLE_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTEX_PARTICLE_EVOLUTION_RLE_3D<double>;
#endif
#endif
