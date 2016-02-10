//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT_LAPLACE.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID>::
IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM(PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input,const T_ARRAYS_SCALAR& variable_viscosity_input,const ARRAY<T>& densities_input,const ARRAY<T>& viscosities_input,T_MPI_GRID* mpi_grid_input,const int axis_input,bool use_variable_viscosity_input)
    :IMPLICIT_VISCOSITY_UNIFORM<T_GRID>(*projection_input.elliptic_solver,variable_viscosity_input,(T)0,(T)0,mpi_grid_input,axis_input,use_variable_viscosity_input,false),projection(projection_input),densities(densities_input),viscosities(viscosities_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID>::
~IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM()
{}
//#####################################################################
// Function Allocate_Heat_Solver
//#####################################################################
template<class T_GRID> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID>::
Allocate_Heat_Solver()
{
    heat_solver=new HEAT_LAPLACE<POISSON_COLLIDABLE_UNIFORM<T_GRID> >(face_grid,u);
}
//#####################################################################
// Function Setup_Viscosity
//#####################################################################
template<class T_GRID> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID>::
Setup_Viscosity(const T dt)
{
    if(use_variable_viscosity) PHYSBAM_NOT_IMPLEMENTED();

    POISSON_COLLIDABLE_UNIFORM<T_GRID>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<T_GRID>&>(*heat_solver);
    PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_dynamics=dynamic_cast<PROJECTION_DYNAMICS_UNIFORM<T_GRID>&>(projection);
    heat_poisson.multiphase=true;
    const int number_of_regions=densities.m;

    // set viscosity coefficients
    ARRAY<T> dt_times_kinematic_viscosity(number_of_regions);
    for(int i=1;i<=number_of_regions;i++) dt_times_kinematic_viscosity(i)=dt*viscosities(i)/densities(i);
    heat_poisson.Set_Constant_beta(dt_times_kinematic_viscosity);

    // set up internal levelset
    heat_poisson.Use_Internal_Level_Set(number_of_regions);
    T_AVERAGING averaging;const T_LEVELSET_MULTIPLE& cell_centered_levelset_multiple=*projection_dynamics.poisson_collidable->levelset_multiple;
    for(CELL_ITERATOR iterator(face_grid,2);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index(),p_face_index=cell_index;
        for(int i=1;i<=number_of_regions;i++) heat_poisson.levelset_multiple->phis(i)(cell_index)=averaging.Cell_To_Face(projection.p_grid,axis,p_face_index,cell_centered_levelset_multiple.phis(i));}
    heat_poisson.levelset_multiple->Project_Levelset(2);

    if(projection_dynamics.flame) Calculate_Velocity_Jump();
    heat_poisson.Find_Constant_beta_Multiphase(heat_poisson.levelset_multiple->phis);
}
//#####################################################################
// Function Setup_Boundary_Conditions
//#####################################################################
template<class T_GRID> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID>::
Setup_Boundary_Conditions(const T_FACE_ARRAYS_SCALAR& face_velocities)
{
    IMPLICIT_VISCOSITY_UNIFORM<T_GRID>::Setup_Boundary_Conditions(face_velocities);
    // set neumann b.c. at zero viscosity faces
    POISSON_COLLIDABLE_UNIFORM<T_GRID>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<T_GRID>&>(*heat_solver);
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){int face_axis=iterator.Axis();TV_INT face=iterator.Face_Index();
        if(!heat_poisson.beta_face(face_axis,face)) heat_poisson.psi_N(face_axis,face)=true;}
}
//#####################################################################
// Function Calculate_Velocity_Jump
//#####################################################################
template<class T_GRID> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID>::
Calculate_Velocity_Jump()
{
    POISSON_COLLIDABLE_UNIFORM<T_GRID>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<T_GRID>&>(*heat_solver);
    PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_dynamics=dynamic_cast<PROJECTION_DYNAMICS_UNIFORM<T_GRID>&>(projection);
    heat_poisson.Set_Jump_Multiphase();
    const ARRAY<TRIPLE<T,T,T> ,VECTOR<int,2> >& flame_speed_constants=projection_dynamics.flame_speed_constants;
    TV_INT axis_offset=TV_INT::Axis_Vector(axis);
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){TV_INT face_index=iterator.Face_Index();TV_INT face_axis_offset=TV_INT::Axis_Vector(iterator.Axis());
        int region_1=heat_poisson.levelset_multiple->Inside_Region(iterator.First_Cell_Index()),region_2=heat_poisson.levelset_multiple->Inside_Region(iterator.Second_Cell_Index());
        // [Vn]=M*[1/density] with M=-density_fuel*flame_speed, [1/density]=(1/density_fuel-1/density_products), flame_speed_constant.z=(-density_fuel*[1/density])
        int fuel_region=region_1,product_region=region_2;
        if(densities(region_1)<densities(region_2)){fuel_region=region_2;product_region=region_1;}
        const TRIPLE<T,T,T>& constants=flame_speed_constants(fuel_region,product_region);if(constants.z==0)continue;
        const T_LEVELSET& levelset=*projection_dynamics.poisson_collidable->levelset_multiple->levelsets(fuel_region);
        T flame_speed=constants.x;
        if(constants.y){T face_curvature;TV_INT p_face_index=face_index;
            if(iterator.Axis()==axis)face_curvature=(*levelset.curvature)(p_face_index-axis_offset);
            else{face_curvature=(T).25*((*levelset.curvature)(p_face_index-axis_offset-face_axis_offset)+(*levelset.curvature)(p_face_index-face_axis_offset)+
                (*levelset.curvature)(p_face_index)+(*levelset.curvature)(p_face_index-axis_offset));}
            flame_speed+=constants.y*face_curvature;}
        T face_normal;TV_INT p_face_index=face_index;
        if(iterator.Axis()==axis)face_normal=(*levelset.normals)(p_face_index-axis_offset)[axis];
        else{face_normal=((*levelset.normals)(p_face_index-axis_offset-face_axis_offset)+(*levelset.normals)(p_face_index-face_axis_offset)+
            (*levelset.normals)(p_face_index)+(*levelset.normals)(p_face_index-axis_offset)).Normalized()[axis];}
        heat_poisson.u_jump_face.Component(iterator.Axis())(face_index)=LEVELSET_MULTIPLE_UNIFORM<T_GRID>::Sign(product_region,fuel_region)*constants.z*flame_speed*face_normal;}
}
//#####################################################################
// Debug_Write
//#####################################################################
template<class T_GRID> void IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID>::
Debug_Write(const std::string& output_directory_input)
{
    POISSON_COLLIDABLE_UNIFORM<T_GRID>& heat_poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<T_GRID>&>(*heat_solver);
    static int frame[3]={0,0,0};
    std::string output_directory=output_directory_input;if(mpi_grid) output_directory+=STRING_UTILITIES::string_sprintf("/processor%d",mpi_grid->rank);
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string output_directory_axis=STRING_UTILITIES::string_sprintf("%s/%d",output_directory.c_str(),axis);FILE_UTILITIES::Create_Directory(output_directory_axis);
    std::string f=STRING_UTILITIES::Value_To_String(frame[axis-1]);
    FILE_UTILITIES::Write_To_File<T>(output_directory_axis+"/grid",face_grid);
    FILE_UTILITIES::Write_To_File<T>(output_directory_axis+"/psi_N."+f,heat_poisson.psi_N);
    FILE_UTILITIES::Write_To_File<T>(output_directory_axis+"/psi_D."+f,heat_poisson.psi_D);
    FILE_UTILITIES::Write_To_File<T>(output_directory_axis+"/colors."+f,heat_poisson.filled_region_colors);
    FILE_UTILITIES::Write_To_File<T>(output_directory_axis+"/beta_face."+f,heat_poisson.beta_face);
    for(int i=1;i<=densities.m;i++){
        std::string filename=STRING_UTILITIES::string_sprintf("/levelset_%d.%s",i,f.c_str());
        FILE_UTILITIES::Write_To_File<T>(output_directory_axis+filename,*heat_poisson.levelset_multiple->levelsets(i));}
    FILE_UTILITIES::Write_To_Text_File(output_directory_axis+"/common/last_frame",frame[axis-1]);frame[axis-1]+=1;
}
//#####################################################################
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<GRID<VECTOR<float,1> > >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<GRID<VECTOR<double,1> > >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >;
template class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
