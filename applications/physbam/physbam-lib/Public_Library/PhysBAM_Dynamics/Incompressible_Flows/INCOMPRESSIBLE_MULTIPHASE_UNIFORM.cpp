//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Frank Losasso, Duc Nguyen, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/minmag.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
#ifdef WIN32
#pragma warning(disable:4723)
#endif
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
INCOMPRESSIBLE_MULTIPHASE_UNIFORM(const T_GRID& grid_input,PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input)
    :INCOMPRESSIBLE_UNIFORM<T_GRID>(grid_input,projection_input),levelset_for_dirichlet_regions(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
~INCOMPRESSIBLE_MULTIPHASE_UNIFORM()
{
}
//#####################################################################
// Function Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Convection(const T dt,const T time,T_FACE_ARRAYS_SCALAR& advecting_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect,const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells)
{
    // TODO: make efficient if advection velocities are same as advected velocities
    T_FACE_ARRAYS_SCALAR advection_face_velocities_ghost;advection_face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Cells_Face(grid,advecting_face_velocities,advection_face_velocities_ghost,time,number_of_ghost_cells);
    T_FACE_ARRAYS_SCALAR face_velocities_to_advect_ghost;face_velocities_to_advect_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,time,number_of_ghost_cells);

    // update convection
    if(pseudo_dirichlet_regions->Number_True()>0){
        T_FACE_ARRAYS_SCALAR face_velocities_liquid=face_velocities_to_advect;T_FACE_ARRAYS_SCALAR advection_face_velocities_ghost_extrapolated=advection_face_velocities_ghost;
        T_FACE_ARRAYS_SCALAR face_velocities_to_advect_ghost_extrapolated=face_velocities_to_advect_ghost;
        T_ARRAYS_SCALAR phi_for_pseudo_dirichlet_regions;T_GRID grid_temp(grid);T_FAST_LEVELSET levelset_for_pseudo_dirichlet_regions(grid_temp,phi_for_pseudo_dirichlet_regions);
        projection.poisson_collidable->levelset_multiple->Get_Single_Levelset(*pseudo_dirichlet_regions,levelset_for_pseudo_dirichlet_regions,false);
        Extrapolate_Velocity_Across_Interface(advection_face_velocities_ghost_extrapolated,phi_for_pseudo_dirichlet_regions);
        Extrapolate_Velocity_Across_Interface(face_velocities_to_advect_ghost_extrapolated,phi_for_pseudo_dirichlet_regions);
        advection->Update_Advection_Equation_Face(grid,face_velocities_liquid,face_velocities_to_advect_ghost_extrapolated,advection_face_velocities_ghost_extrapolated,*boundary,dt,time);
        advection->Update_Advection_Equation_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,advection_face_velocities_ghost,*boundary,dt,time);
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            int region1=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.First_Cell_Index());
            int region2=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.Second_Cell_Index());
            if(!(*pseudo_dirichlet_regions)(region1)||!(*pseudo_dirichlet_regions)(region2))
                face_velocities_to_advect.Component(axis)(face_index)=face_velocities_liquid.Component(axis)(face_index);}}
    else advection->Update_Advection_Equation_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,advection_face_velocities_ghost,*boundary,dt,time);
}
//#####################################################################
// Function Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity,const ARRAY<T_ARRAYS_SCALAR>* phi_ghost,const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost(grid.Domain_Indices(number_of_ghost_cells),false);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);

    // update strain and apply elastic forces
    for(int i=1;i<=strains.m;i++)if(strains(i)){
        // extrapolate the velocity across the interface to get better strain boundaries
        T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(number_of_ghost_cells));projection.poisson_collidable->levelset_multiple->levelsets(i)->boundary->Fill_Ghost_Cells(grid,projection.poisson_collidable->levelset_multiple->phis(i),phi_ghost,dt,time,number_of_ghost_cells);
        T_FACE_ARRAYS_SCALAR face_velocities_temp=face_velocities_ghost;
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            T_GRID face_grid=grid.Get_Face_Grid(axis);T_ARRAYS_SCALAR phi_face(face_grid.Domain_Indices(),false);T_ARRAYS_BASE& face_velocity=face_velocities_temp.Component(axis);
            for(FACE_ITERATOR iterator(grid,0,T_GRID::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next())
                phi_face(iterator.Face_Index())=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
            int extrapolation_bandwidth=3;
            T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,number_of_ghost_cells);extrapolate.Set_Band_Width((T)extrapolation_bandwidth);
            extrapolate.Extrapolate();}
        strains(i)->Update_Strain_Equation_Multiphase(dt,time,projection.densities(i),face_velocities,face_velocities_temp,*projection.poisson_collidable->levelset_multiple,i,number_of_ghost_cells);}

    // update gravity
    if(gravity) for(int axis=1;axis<=T_GRID::dimension;axis++) if(downward_direction[axis]) face_velocities.Component(axis)+=dt*gravity*downward_direction[axis];
    
    // update body force
    if(use_force) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) 
        face_velocities.Component(iterator.Axis())(iterator.Face_Index())+=dt*force.Component(iterator.Axis())(iterator.Face_Index());

    if((viscosity || use_variable_viscosity) && (!implicit_viscosity || use_explicit_part_of_implicit_viscosity))
        Discretize_Explicit_Viscous_Terms(dt);

    if(!GFM && nonzero_surface_tension){
        T half_width=(T).5*number_of_interface_cells*grid.Minimum_Edge_Length();T dt_over_four=(T).25*dt;
        LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;
        levelset_multiple.Compute_Normals();levelset_multiple.Compute_Curvature();
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
            T sign=(region_1==region_2)?(T)1:(T)-1;
            T phi_face=(T).5*(phi_1+sign*phi_2);
            T twice_curvature=((*levelset_multiple.levelsets(region_1)->curvature)(cell_1)+sign*(*levelset_multiple.levelsets(region_2)->curvature)(cell_2));
            T twice_face_normal_component=((*levelset_multiple.levelsets(region_1)->normals)(cell_1)[axis]+sign*(*levelset_multiple.levelsets(region_2)->normals)(cell_2)[axis]);
            face_velocities.Component(axis)(index)+=dt_over_four*LEVELSET_UTILITIES<T>::Delta(phi_face,half_width)*surface_tensions(region_1,region_2)*
                twice_curvature*twice_face_normal_component/LEVELSET_UTILITIES<T>::Heaviside(phi_face,projection.densities(region_1),projection.densities(region_2),half_width);}}

    if(use_variable_vorticity_confinement){
        T_ARRAYS_VECTOR F(grid.Cell_Indices(1),false);
        Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
        F*=dt*(T).5;F*=variable_vorticity_confinement;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
            face_velocities.Component(axis)(iterator.Face_Index())+=F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis];}}

    if(vorticity_confinements.Count_Matches(0)!=vorticity_confinements.m){
        T_ARRAYS_VECTOR F(grid.Cell_Indices(1),false);
        Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
        T half_dt=(T).5*dt;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();
            if(pseudo_dirichlet_regions){
                int region1=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.First_Cell_Index());
                int region2=projection.poisson_collidable->levelset_multiple->Inside_Region(iterator.Second_Cell_Index());
                if(!(*pseudo_dirichlet_regions)(region1) || !(*pseudo_dirichlet_regions)(region2))
                    if((!(*pseudo_dirichlet_regions)(region1) && vorticity_confinements(region1)==0) || (!(*pseudo_dirichlet_regions)(region2) && vorticity_confinements(region2)==0))
                        continue;}
            int region=projection.poisson_collidable->levelset_multiple->Inside_Region_Face(iterator.Axis(),iterator.Face_Index());
            T face_F=F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis];
            face_velocities.Component(axis)(iterator.Face_Index())+=face_F*half_dt*vorticity_confinements(region);}}

    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Implicit_Part(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity,BOUNDARY_UNIFORM<T_GRID,T>* projection_boundary,
    bool use_levelset_viscosity,BOUNDARY_CONDITIONS_CALLBACKS<TV>* bc_callbacks,bool print_viscosity_matrix)
{
    // boundary conditions
    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);

    // set up poisson equation
    ARRAY<T> one_over_densities(projection.densities.m);for(int i=1;i<=projection.densities.m;i++)one_over_densities(i)=1/projection.densities(i);
    projection.poisson->Set_Constant_beta(one_over_densities);

    if(!GFM){projection.poisson->Use_Delta_Function_Method(number_of_interface_cells);projection.poisson->Smear_One_Over_beta();}

// SAVE THIS HERE...
//    if(viscosity_minus || viscosity_plus){ // update velocity to include viscosity on the interior
//        int ghost_cells=3;
//        ARRAY<T,VECTOR<int,3> > phi_ghost(grid.Domain_Indices(3));levelset.boundary->Fill_Ghost_Cells(grid,levelset.phi,phi_ghost,dt,time,ghost_cells);
//        projection.poisson->Find_Constant_beta(phi_ghost);
//        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
//            T dt_over_rho=.125*dt*(projection.poisson->beta_right(i-1,j,ij)+projection.poisson->beta_right(i,j,ij)+projection.poisson->beta_top(i,j-1,ij)+projection.poisson->beta_top(i,j,ij)
//                +projection.poisson->beta_back(i,j,ij-1)+projection.poisson->beta_back(i,j,ij));
//            V(i,j,ij).x+=dt_over_rho*u_viscosity(i,j,ij);V(i,j,ij).y+=dt_over_rho*v_viscosity(i,j,ij);V(i,j,ij).z+=dt_over_rho*w_viscosity(i,j,ij);}

    if((GFM&&nonzero_surface_tension)||projection.flame) projection.poisson_collidable->Set_Jump_Multiphase();
    
    if(GFM && nonzero_surface_tension){
        LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;levelset_multiple.Compute_Curvature();
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
            if(region_1!=region_2){
                T sign=LEVELSET_MULTIPLE_UNIFORM<T_GRID>::Sign(region_1,region_2);
                T curvature=LEVELSET_UTILITIES<T>::Average(phi_1,-sign*(*levelset_multiple.levelsets(region_1)->curvature)(cell_1),
                                                       -phi_2,sign*(*levelset_multiple.levelsets(region_2)->curvature)(cell_2));
                projection.poisson_collidable->u_jump_face.Component(iterator.Axis())(iterator.Face_Index())+=dt*surface_tensions(region_1,region_2)*curvature;}}}

    if(projection.flame) Calculate_Pressure_Jump(dt,time);

    // viscosity
    if(nonzero_viscosity && implicit_viscosity){
        projection.Make_Divergence_Free(face_velocities,dt,time);Implicit_Viscous_Update(face_velocities,dt,time);}

    // make divergence free
    projection.Make_Divergence_Free(face_velocities,dt,time);
}
//#####################################################################
// Function Calculate_Pressure_Jump
//#####################################################################
// flame_speed must be up to date
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Calculate_Pressure_Jump(const T dt,const T time)
{
    assert(projection.poisson_collidable->levelset_multiple);
    LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        int region1=levelset_multiple.Inside_Region(cell1),region2=levelset_multiple.Inside_Region(cell2);
        const TRIPLE<T,T,T>& constants=projection.flame_speed_constants(region1,region2);
        if(constants.z==0)continue;
        if(projection.densities(region1)<projection.densities(region2))exchange(region1,region2); //region1 is now the fuel region
        // [p]=dt*[1/density]*sqr(M) with M=-density_fuel*flame_speed, [1/density]=(1/density_fuel-1/density_products)
        // flame_speed_constant.z is (-density_fuel*[1/density])
        projection.poisson_collidable->u_jump_face.Component(iterator.Axis())(iterator.Face_Index())+=LEVELSET_MULTIPLE_UNIFORM<T_GRID>::Sign(region1,region2)*
            dt*constants.z*projection.densities(region1)*sqr(projection.Flame_Speed_Face_Multiphase(iterator.Axis(),iterator.Face_Index(),region1,region2));}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
CFL(T_FACE_ARRAYS_SCALAR& face_velocities,const bool inviscid,const bool viscous_only) const
{
    TV DX=grid.dX,sqr_DX=DX*DX,max_abs_V;
    // convection
    T dt_convection=0;
    if(!viscous_only){
        if(projection.flame){
            FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> face_velocities_fire(face_velocities,projection,projection.poisson_collidable->levelset_multiple);
            typename FIRE_INTERPOLATION_POLICY<T_GRID>::AVERAGING_FIRE_MULTIPHASE averaging;
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV V=averaging.Face_To_Face_Vector(grid,iterator.Axis(),iterator.Face_Index(),face_velocities_fire);
                dt_convection=max(dt_convection,TV::Dot_Product(grid.one_over_dX,abs(V)));}}
        else{
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                T local_V_norm=0;for(int axis=1;axis<=T_GRID::dimension;axis++){
                    TV_INT first_face_index=grid.First_Face_Index_In_Cell(axis,cell),second_face_index=grid.Second_Face_Index_In_Cell(axis,cell);
                    local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,first_face_index),face_velocities(axis,second_face_index));}
                dt_convection=max(dt_convection,local_V_norm);}}}
    // viscosity
    T dt_viscosity=0;
    if(!inviscid){
        T norm_2_over_sqr_DX=2*Inverse(sqr_DX).L1_Norm();
        for(int i=1;i<=viscosities.m;i++)dt_viscosity=max(dt_viscosity,viscosities(i)/projection.densities(i));
        dt_viscosity*=norm_2_over_sqr_DX;
        if(use_variable_viscosity) PHYSBAM_NOT_IMPLEMENTED();}
    // surface tension
    T dt_surface_tension=0;
    if(nonzero_surface_tension){
        LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;levelset_multiple.Compute_Curvature();
        int ghost_cells=1;
        levelset_multiple.Fill_Ghost_Cells(levelset_multiple.phis,0,ghost_cells);T kappa_cfl=0;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
            if(region_1!=region_2 && surface_tensions(region_1,region_2)){
                T curvature=LEVELSET_UTILITIES<T>::Average(phi_1,LEVELSET_MULTIPLE_UNIFORM<T_GRID>::Sign(region_2,region_1)*(*levelset_multiple.levelsets(region_1)->curvature)(cell_1),
                    -phi_2,LEVELSET_MULTIPLE_UNIFORM<T_GRID>::Sign(region_1,region_2)*(*levelset_multiple.levelsets(region_2)->curvature)(cell_2));
                kappa_cfl=max(kappa_cfl,abs(curvature*surface_tensions(region_1,region_2)/
                          LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_1-phi_2),projection.densities(region_1),projection.densities(region_2))));}}
        dt_surface_tension=sqrt(kappa_cfl)/grid.Minimum_Edge_Length();}
    TV max_force;
    if(use_force) max_force=force.Maxabs();
    T dt_force=0;
    if(use_force) dt_force=(max_force/DX).L1_Norm();
    if(gravity) dt_force+=abs(gravity)*(downward_direction/DX).L1_Norm();
    T strain_cfl=FLT_MAX;
    if(strain){
        for(int i=1;i<=strains.m;i++)if(strains(i))
            strain_cfl=min(strain_cfl,strains(i)->CFL(projection.densities(i)));}
    dt_force+=1/strain_cfl;
    T dt_overall=(dt_convection+dt_viscosity+sqrt(sqr(dt_convection+dt_viscosity)+4*dt_force+4*sqr(dt_surface_tension)))/2; 
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Set_Dirichlet_Boundary_Conditions(ARRAY<T_ARRAYS_SCALAR>& phis,const ARRAY<bool>& dirichlet_regions,const ARRAY<T>* pressures)
{
    LEVELSET_MULTIPLE_UNIFORM<T_GRID> levelset_multiple(grid,phis);
    if(dirichlet_regions.Number_True()>0){
        if(!pressures){for(CELL_ITERATOR iterator(projection.p_grid);iterator.Valid();iterator.Next()) if(dirichlet_regions(levelset_multiple.Inside_Region(iterator.Cell_Index()))){
            projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=0;}}
        else{for(CELL_ITERATOR iterator(projection.p_grid);iterator.Valid();iterator.Next()){
            int region=levelset_multiple.Inside_Region(iterator.Cell_Index());
            if(dirichlet_regions(levelset_multiple.Inside_Region(iterator.Cell_Index()))){
                projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=(*pressures)(region);}}}}
    if(mpi_grid){
        mpi_grid->Exchange_Boundary_Cell_Data(projection.elliptic_solver->psi_D,1,nonzero_viscosity); // need to exchange corners only in the case of implicit viscosity
        mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
}
//#####################################################################
// Function Add_Surface_Tension
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Add_Surface_Tension(T_LEVELSET& levelset,const T time)
{
    LAPLACE_UNIFORM<T_GRID>& elliptic_solver=*projection.elliptic_solver;T_GRID& p_grid=elliptic_solver.grid;
    LAPLACE_COLLIDABLE<T_GRID>& collidable_solver=*projection.collidable_solver;
    LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple=*projection.poisson_collidable->levelset_multiple;
    T_ARRAYS_SCALAR& phi=levelset.phi; 

   if(collidable_solver.second_order_cut_cell_method) for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index(),face_index=iterator.Face_Index();int axis=iterator.Axis();
        if(!projection.elliptic_solver->psi_N.Component(axis)(face_index) && LEVELSET_UTILITIES<T>::Interface(phi(first_cell_index),phi(second_cell_index))){
            int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(first_cell_index,second_cell_index,region_1,region_2,phi_1,phi_2);
            if(!surface_tensions(region_1,region_2))continue;
            T theta=LEVELSET_UTILITIES<T>::Theta(phi(first_cell_index),phi(second_cell_index));
            TV location=theta*(grid.Center(second_cell_index)-grid.Center(first_cell_index))+grid.Center(first_cell_index);
            T curvature_at_interface=levelset.Compute_Curvature(location);
            collidable_solver.u_interface.Component(axis)(face_index)=-surface_tensions(region_1,region_2)*curvature_at_interface;}}
   else{
       for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()) if(elliptic_solver.psi_D(iterator.Cell_Index()) && phi(iterator.Cell_Index()) < 5*grid.dX.Max()){
           int minimum_region,second_minimum_region;T minimum_phi,second_minimum_phi;
           levelset_multiple.Two_Minimum_Regions(iterator.Cell_Index(),minimum_region,second_minimum_region,minimum_phi,second_minimum_phi);
           projection.p(iterator.Cell_Index())=-surface_tensions(minimum_region,second_minimum_region)*levelset.Compute_Curvature(phi,iterator.Cell_Index());}}
}
//#####################################################################
// Function Implicit_Viscous_Update
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Implicit_Viscous_Update(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM<T_GRID> implicit_viscosity(projection,variable_viscosity,projection.densities,viscosities,mpi_grid,axis,use_variable_viscosity);
        implicit_viscosity.Viscous_Update(grid,face_velocities,face_velocities,dt,time,maximum_implicit_viscosity_iterations);}
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Compute_Vorticity_Confinement_Force(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_ARRAYS_VECTOR& F)
{
    T_ARRAYS_SPIN vorticity(grid.Cell_Indices(2),false);T_ARRAYS_SCALAR vorticity_magnitude(grid.Cell_Indices(2));
    if(projection.flame){FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> face_velocities_lookup(face_velocities_ghost,projection,projection.poisson_collidable->levelset_multiple);
        VORTICITY_UNIFORM<TV>::Vorticity(grid,face_velocities_lookup,vorticity,vorticity_magnitude);}
    else VORTICITY_UNIFORM<TV>::Vorticity(grid,FACE_LOOKUP_UNIFORM<T_GRID>(face_velocities_ghost),vorticity,vorticity_magnitude);
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        TV vortex_normal_vector=T_LEVELSET::Normal_At_Node(grid,vorticity_magnitude,iterator.Cell_Index());
        F(iterator.Cell_Index())=TV::Cross_Product(vortex_normal_vector,vorticity(iterator.Cell_Index()));}
}
//#####################################################################
// Function Discretize_Viscous_Terms ------ OLD 3D CODE!!!!!!
//#####################################################################
/*
template<class T_GRID> void INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>::
Discretize_Viscous_Terms(const ARRAY<T,VECTOR<int,3> >& phi,const T dt)
{

    int m=grid.m,n=grid.n,mn=grid.mn;T dx=grid.dx,dy=grid.dy,dz=grid.dz;
    T half_width=number_of_interface_cells*grid.max_dx_dy_dz/2;
    int i,j,ij;

    ARRAY<T,VECTOR<int,3> > viscosity_x_half(0,m,1,n,1,mn),ux_x_half(0,m,1,n,1,mn),vx_x_half(0,m,1,n,1,mn),wx_x_half(0,m,1,n,1,mn);
    ARRAY<T,VECTOR<int,3> > viscosity_y_half(1,m,0,n,1,mn),uy_y_half(1,m,0,n,1,mn),vy_y_half(1,m,0,n,1,mn),wy_y_half(1,m,0,n,1,mn);
    ARRAY<T,VECTOR<int,3> > viscosity_z_half(1,m,1,n,0,mn),uz_z_half(1,m,1,n,0,mn),vz_z_half(1,m,1,n,0,mn),wz_z_half(1,m,1,n,0,mn);
    for(i=0;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++){
        viscosity_x_half(i,j,ij)=LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside((T).5*(phi(i,j,ij)+phi(i+1,j,ij)),viscosity_minus,viscosity_plus,half_width);
        ux_x_half(i,j,ij)=(V_ghost(i+1,j,ij).x-V_ghost(i,j,ij).x)/dx;
        vx_x_half(i,j,ij)=(V_ghost(i+1,j,ij).y-V_ghost(i,j,ij).y)/dx;
        wx_x_half(i,j,ij)=(V_ghost(i+1,j,ij).z-V_ghost(i,j,ij).z)/dx;}
    for(i=1;i<=m;i++) for(j=0;j<=n;j++) for(ij=1;ij<=mn;ij++){
        viscosity_y_half(i,j,ij)=LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside((T).5*(phi(i,j,ij)+phi(i,j+1,ij)),viscosity_minus,viscosity_plus,half_width);
        uy_y_half(i,j,ij)=(V_ghost(i,j+1,ij).x-V_ghost(i,j,ij).x)/dy;
        vy_y_half(i,j,ij)=(V_ghost(i,j+1,ij).y-V_ghost(i,j,ij).y)/dy;
        wy_y_half(i,j,ij)=(V_ghost(i,j+1,ij).z-V_ghost(i,j,ij).z)/dy;}
    for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=0;ij<=mn;ij++){
        viscosity_z_half(i,j,ij)=LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside((T).5*(phi(i,j,ij)+phi(i,j,ij+1)),viscosity_minus,viscosity_plus,half_width);
        uz_z_half(i,j,ij)=(V_ghost(i,j,ij+1).x-V_ghost(i,j,ij).x)/dz;
        vz_z_half(i,j,ij)=(V_ghost(i,j,ij+1).y-V_ghost(i,j,ij).y)/dz;
        wz_z_half(i,j,ij)=(V_ghost(i,j,ij+1).z-V_ghost(i,j,ij).z)/dz;}

    if(!GFM){
        ARRAY<T,VECTOR<int,3> > uy_x_half(0,m,1,n,1,mn),uz_x_half(0,m,1,n,1,mn),vx_y_half(1,m,0,n,1,mn),vz_y_half(1,m,0,n,1,mn),wx_z_half(1,m,1,n,0,mn),wy_z_half(1,m,1,n,0,mn);
        {ARRAY<T,VECTOR<int,3> > uy_nodes(0,m+1,1,n,1,mn),uz_nodes(0,m+1,1,n,1,mn);
        for(i=0;i<=m+1;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++){
            uy_nodes(i,j,ij)=(V_ghost(i,j+1,ij).x-V_ghost(i,j-1,ij).x)/(2*dy);uz_nodes(i,j,ij)=(V_ghost(i,j,ij+1).x-V_ghost(i,j,ij-1).x)/(2*dz);}
        for(i=0;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++){
            uy_x_half(i,j,ij)=(T).5*(uy_nodes(i,j,ij)+uy_nodes(i+1,j,ij));uz_x_half(i,j,ij)=(T).5*(uz_nodes(i,j,ij)+uz_nodes(i+1,j,ij));}}
        {ARRAY<T,VECTOR<int,3> > vx_nodes(1,m,0,n+1,1,mn),vz_nodes(1,m,0,n+1,1,mn);
        for(i=1;i<=m;i++) for(j=0;j<=n+1;j++) for(ij=1;ij<=mn;ij++){
            vx_nodes(i,j,ij)=(V_ghost(i+1,j,ij).y-V_ghost(i-1,j,ij).y)/(2*dx);vz_nodes(i,j,ij)=(V_ghost(i,j,ij+1).y-V_ghost(i,j,ij-1).y)/(2*dz);}
        for(i=1;i<=m;i++) for(j=0;j<=n;j++) for(ij=1;ij<=mn;ij++){
            vx_y_half(i,j,ij)=(T).5*(vx_nodes(i,j,ij)+vx_nodes(i,j+1,ij));vz_y_half(i,j,ij)=(T).5*(vz_nodes(i,j,ij)+vz_nodes(i,j+1,ij));}}
        {ARRAY<T,VECTOR<int,3> > wx_nodes(1,m,1,n,0,mn+1),wy_nodes(1,m,1,n,0,mn+1);
        for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=0;ij<=mn+1;ij++){
            wx_nodes(i,j,ij)=(V_ghost(i+1,j,ij).z-V_ghost(i-1,j,ij).z)/(2*dx);wy_nodes(i,j,ij)=(V_ghost(i,j+1,ij).z-V_ghost(i,j-1,ij).z)/(2*dy);}     
        for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=0;ij<=mn;ij++){
            wx_z_half(i,j,ij)=(T).5*(wx_nodes(i,j,ij)+wx_nodes(i,j,ij+1));wy_z_half(i,j,ij)=(T).5*(wy_nodes(i,j,ij)+wy_nodes(i,j,ij+1));}

        // (2*viscosity*ux)x + (viscosity*(uy+vx))y + (viscosity*(uz+wx))z
        for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++)
            u_viscosity(i,j,ij)=2*(viscosity_x_half(i,j,ij)*ux_x_half(i,j,ij)-viscosity_x_half(i-1,j,ij)*ux_x_half(i-1,j,ij))/dx+
                                       (viscosity_y_half(i,j,ij)*(uy_y_half(i,j,ij)+vx_y_half(i,j,ij))-viscosity_y_half(i,j-1,ij)*(uy_y_half(i,j-1,ij)+vx_y_half(i,j-1,ij)))/dy+
                                       (viscosity_z_half(i,j,ij)*(uz_z_half(i,j,ij)+wx_z_half(i,j,ij))-viscosity_z_half(i,j,ij-1)*(uz_z_half(i,j,ij-1)+wx_z_half(i,j,ij-1)))/dz;
        // (viscosity*(uy+vx))x + (2*viscosity*vy)y + (viscosity*(vz+wy))z
        for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++)
            v_viscosity(i,j,ij)=(viscosity_x_half(i,j,ij)*(uy_x_half(i,j,ij)+vx_x_half(i,j,ij))-viscosity_x_half(i-1,j,ij)*(uy_x_half(i-1,j,ij)+vx_x_half(i-1,j,ij)))/dx+
                                       2*(viscosity_y_half(i,j,ij)*vy_y_half(i,j,ij)-viscosity_y_half(i,j-1,ij)*vy_y_half(i,j-1,ij))/dy+
                                       (viscosity_z_half(i,j,ij)*(vz_z_half(i,j,ij)+wy_z_half(i,j,ij))-viscosity_z_half(i,j,ij-1)*(vz_z_half(i,j,ij-1)+wy_z_half(i,j,ij-1)))/dz;
        // (viscosity*(uz+wx))x + (viscosity(vz+wy))y + (2*viscosity*wz)z
        for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++)
            w_viscosity(i,j,ij)=(viscosity_x_half(i,j,ij)*(uz_x_half(i,j,ij)+wx_x_half(i,j,ij))-viscosity_x_half(i-1,j,ij)*(uz_x_half(i-1,j,ij)+wx_x_half(i-1,j,ij)))/dx+
                                        (viscosity_y_half(i,j,ij)*(vz_y_half(i,j,ij)+wy_y_half(i,j,ij))-viscosity_y_half(i,j-1,ij)*(vz_y_half(i,j-1,ij)+wy_y_half(i,j-1,ij)))/dy+
                                        2*(viscosity_z_half(i,j,ij)*wz_z_half(i,j,ij)-viscosity_z_half(i,j,ij-1)*wz_z_half(i,j,ij-1))/dz;}
    else{ // GFM
        levelset.Compute_Normals();
        ARRAY<T,VECTOR<int,3> > ux_jump(0,m+1,0,n+1,0,mn+1),uy_jump(0,m+1,0,n+1,0,mn+1),uz_jump(0,m+1,0,n+1,0,mn+1),vx_jump(0,m+1,0,n+1,0,mn+1),vy_jump(0,m+1,0,n+1,0,mn+1),vz_jump(0,m+1,0,n+1,0,mn+1),wx_jump(0,m+1,0,n+1,0,mn+1),wy_jump(0,m+1,0,n+1,0,mn+1),wz_jump(0,m+1,0,n+1,0,mn+1);
        T viscosity_jump=viscosity_plus-viscosity_minus;
        for(i=0;i<=m+1;i++) for(j=0;j<=n+1;j++) for(ij=0;ij<=mn+1;ij++){
            MATRIX<T,3> UX((V_ghost(i+1,j,ij)-V_ghost(i-1,j,ij))/(2*dx),(V_ghost(i,j+1,ij)-V_ghost(i,j-1,ij))/(2*dy),(V_ghost(i,j,ij+1)-V_ghost(i,j,ij-1))/(2*dz));
            VECTOR<T,3> N((*levelset.normals)(i,j,ij)),T1(N.Unit_Orthogonal_Vector()),T2(Cross_Product(N,T1));
            SYMMETRIC_MATRIX<T,3> T_transpose_T(MATRIX<T,3>(VECTOR<T,3>(),T1,T2).Transposed().Normal_Equations_Matrix()),N_transpose_N(SYMMETRIC_MATRIX<T,3>::Outer_Product(N));
            MATRIX<T,3> JUMP=viscosity_jump*(UX*T_transpose_T+(N_transpose_N*UX-T_transpose_T*UX.Transposed())*N_transpose_N);
            ux_jump(i,j,ij)=JUMP.x[0];uy_jump(i,j,ij)=JUMP.x[3];uz_jump(i,j,ij)=JUMP.x[6];
            vx_jump(i,j,ij)=JUMP.x[1];vy_jump(i,j,ij)=JUMP.x[4];vz_jump(i,j,ij)=JUMP.x[7];
            wx_jump(i,j,ij)=JUMP.x[2];wy_jump(i,j,ij)=JUMP.x[5];wz_jump(i,j,ij)=JUMP.x[8];}

        // find pressure jump condition
        projection.poisson->Set_Jump(); // allocates the memory for projection.poisson->u_jump
        projection.poisson->levelset->Compute_Normals(); // cell center normals
        for(int i=1;i<=m-1;i++) for(int j=1;j<=n-1;j++) for(int ij=1;ij<=mn-1;ij++){
            T ux=(T).25*(ux_x_half(i,j,ij)+ux_x_half(i,j+1,ij)+ux_x_half(i,j,ij+1)+ux_x_half(i,j+1,ij+1)),uy=(T).25*(uy_y_half(i,j,ij)+uy_y_half(i+1,j,ij)+uy_y_half(i,j,ij+1)+uy_y_half(i+1,j,ij+1)),uz=(T).25*(uz_z_half(i,j,ij)+uz_z_half(i+1,j,ij)+uz_z_half(i,j+1,ij)+uz_z_half(i+1,j+1,ij)),
              vx=(T).25*(vx_x_half(i,j,ij)+vx_x_half(i,j+1,ij)+vx_x_half(i,j,ij+1)+vx_x_half(i,j+1,ij+1)),vy=(T).25*(vy_y_half(i,j,ij)+vy_y_half(i+1,j,ij)+vy_y_half(i,j,ij+1)+vy_y_half(i+1,j,ij+1)),vz=(T).25*(vz_z_half(i,j,ij)+vz_z_half(i+1,j,ij)+vz_z_half(i,j+1,ij)+vz_z_half(i+1,j+1,ij)),
              wx=(T).25*(wx_x_half(i,j,ij)+wx_x_half(i,j+1,ij)+wx_x_half(i,j,ij+1)+wx_x_half(i,j+1,ij+1)),wy=(T).25*(wy_y_half(i,j,ij)+wy_y_half(i+1,j,ij)+wy_y_half(i,j,ij+1)+wy_y_half(i+1,j,ij+1)),wz=(T).25*(wz_z_half(i,j,ij)+wz_z_half(i+1,j,ij)+wz_z_half(i,j+1,ij)+wz_z_half(i+1,j+1,ij));
            MATRIX<T,2> UX(ux,vx,wx,uy,vy,wy,uz,vz,wz);
            VECTOR<T,3> N((*projection.poisson->levelset->normals)(i,j,ij));
            projection.poisson->u_jump(i,j,ij)=2*dt*viscosity_jump*VECTOR<T,3>::Dot_Product(UX*N,N);}

        // update u_viscosity = viscosity*uxx + viscosity*uyy + viscosity*uzz
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){ 
            T viscosity_ux_left,viscosity_ux_right,viscosity_uy_bottom,viscosity_uy_top,viscosity_uz_front,viscosity_uz_back;
            T phi_mid=phi(i,j,ij),phi_left=phi(i-1,j,ij),phi_right=phi(i+1,j,ij),phi_bottom=phi(i,j-1,ij),phi_top=phi(i,j+1,ij),phi_front=phi(i,j,ij-1),phi_back=phi(i,j,ij+1);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_left,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_left,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_left,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_left,ux_jump(i-1,j,ij),phi_mid,ux_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_left,phi_mid);
                viscosity_ux_left=jump_viscosity*ux_x_half(i-1,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_ux_left=viscosity_x_half(i-1,j,ij)*ux_x_half(i-1,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_right)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_right,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_right,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_right,ux_jump(i+1,j,ij),phi_mid,ux_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_right,phi_mid);
                viscosity_ux_right=jump_viscosity*ux_x_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_ux_right=viscosity_x_half(i,j,ij)*ux_x_half(i,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_bottom,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_bottom,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_bottom,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_bottom,uy_jump(i,j-1,ij),phi_mid,uy_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_bottom,phi_mid);
                viscosity_uy_bottom=jump_viscosity*uy_y_half(i,j-1,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_uy_bottom=viscosity_y_half(i,j-1,ij)*uy_y_half(i,j-1,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_top)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_top,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_top,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_top,uy_jump(i,j+1,ij),phi_mid,uy_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_top,phi_mid);
                viscosity_uy_top=jump_viscosity*uy_y_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_uy_top=viscosity_y_half(i,j,ij)*uy_y_half(i,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_front,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_front,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_front,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_front,uz_jump(i,j,ij-1),phi_mid,uz_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_front,phi_mid);
                viscosity_uz_front=jump_viscosity*uz_z_half(i,j,ij-1)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_uz_front=viscosity_z_half(i,j,ij-1)*uz_z_half(i,j,ij-1);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_back)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_back,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_back,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_back,uz_jump(i,j,ij+1),phi_mid,uz_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_back,phi_mid);
                viscosity_uz_back=jump_viscosity*uz_z_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_uz_back=viscosity_z_half(i,j,ij)*uz_z_half(i,j,ij);
            u_viscosity(i,j,ij)=(viscosity_ux_right-viscosity_ux_left)/dx+(viscosity_uy_top-viscosity_uy_bottom)/dy
                               +(viscosity_uz_back-viscosity_uz_front)/dz;}
            
        // update v_viscosity = viscosity*vxx + viscosity*vyy + viscosity*vzz
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
            T viscosity_vx_left,viscosity_vx_right,viscosity_vy_bottom,viscosity_vy_top,viscosity_vz_front,viscosity_vz_back;
            T phi_mid=phi(i,j,ij),phi_left=phi(i-1,j,ij),phi_right=phi(i+1,j,ij),phi_bottom=phi(i,j-1,ij),phi_top=phi(i,j+1,ij),phi_front=phi(i,j,ij-1),phi_back=phi(i,j,ij+1);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_left,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_left,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_left,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_left,vx_jump(i-1,j,ij),phi_mid,vx_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_left,phi_mid);
                viscosity_vx_left=jump_viscosity*vx_x_half(i-1,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_vx_left=viscosity_x_half(i-1,j,ij)*vx_x_half(i-1,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_right)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_right,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_right,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_right,vx_jump(i+1,j,ij),phi_mid,vx_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_right,phi_mid);
                viscosity_vx_right=jump_viscosity*vx_x_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_vx_right=viscosity_x_half(i,j,ij)*vx_x_half(i,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_bottom,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_bottom,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_bottom,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_bottom,vy_jump(i,j-1,ij),phi_mid,vy_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_bottom,phi_mid);
                viscosity_vy_bottom=jump_viscosity*vy_y_half(i,j-1,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_vy_bottom=viscosity_y_half(i,j-1,ij)*vy_y_half(i,j-1,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_top)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_top,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_top,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_top,vy_jump(i,j+1,ij),phi_mid,vy_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_top,phi_mid);
                viscosity_vy_top=jump_viscosity*vy_y_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_vy_top=viscosity_y_half(i,j,ij)*vy_y_half(i,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_front,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_front,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_front,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_front,vz_jump(i,j,ij-1),phi_mid,vz_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_front,phi_mid);
                viscosity_vz_front=jump_viscosity*vz_z_half(i,j,ij-1)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_vz_front=viscosity_z_half(i,j,ij-1)*vz_z_half(i,j,ij-1);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_back)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_back,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_back,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_back,vz_jump(i,j,ij+1),phi_mid,vz_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_back,phi_mid);
                viscosity_vz_back=jump_viscosity*vz_z_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_vz_back=viscosity_z_half(i,j,ij)*vz_z_half(i,j,ij);
            v_viscosity(i,j,ij)=(viscosity_vx_right-viscosity_vx_left)/dx+(viscosity_vy_top-viscosity_vy_bottom)/dy
                               +(viscosity_vz_back-viscosity_vz_front)/dz;}

        // update w_viscosity = viscosity*wxx + viscosity*wyy + viscosity*wzz
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
            T viscosity_wx_left,viscosity_wx_right,viscosity_wy_bottom,viscosity_wy_top,viscosity_wz_front,viscosity_wz_back;
            T phi_mid=phi(i,j,ij),phi_left=phi(i-1,j,ij),phi_right=phi(i+1,j,ij),phi_bottom=phi(i,j-1,ij),phi_top=phi(i,j+1,ij),phi_front=phi(i,j,ij-1),phi_back=phi(i,j,ij+1);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_left,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_left,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_left,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_left,wx_jump(i-1,j,ij),phi_mid,wx_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_left,phi_mid);
                viscosity_wx_left=jump_viscosity*wx_x_half(i-1,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_wx_left=viscosity_x_half(i-1,j,ij)*wx_x_half(i-1,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_right)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_right,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_right,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_right,wx_jump(i+1,j,ij),phi_mid,wx_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_right,phi_mid);
                viscosity_wx_right=jump_viscosity*wx_x_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_wx_right=viscosity_x_half(i,j,ij)*wx_x_half(i,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_bottom,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_bottom,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_bottom,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_bottom,wy_jump(i,j-1,ij),phi_mid,wy_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_bottom,phi_mid);
                viscosity_wy_bottom=jump_viscosity*wy_y_half(i,j-1,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_wy_bottom=viscosity_y_half(i,j-1,ij)*wy_y_half(i,j-1,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_top)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_top,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_top,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_top,wy_jump(i,j+1,ij),phi_mid,wy_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_top,phi_mid);
                viscosity_wy_top=jump_viscosity*wy_y_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_wy_top=viscosity_y_half(i,j,ij)*wy_y_half(i,j,ij);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_front,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_front,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_front,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_front,wz_jump(i,j,ij-1),phi_mid,wz_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_front,phi_mid);
                viscosity_wz_front=jump_viscosity*wz_z_half(i,j,ij-1)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_wz_front=viscosity_z_half(i,j,ij-1)*wz_z_half(i,j,ij-1);
            if(LEVELSET_UTILITIES<VECTOR<T,3> >::Interface(phi_mid,phi_back)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,3> >::Convex_Average(phi_mid,phi_back,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,3> >::Heaviside(phi_back,viscosity_minus,viscosity_plus)*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Average(phi_back,wz_jump(i,j,ij+1),phi_mid,wz_jump(i,j,ij))*
                            LEVELSET_UTILITIES<VECTOR<T,3> >::Theta(phi_back,phi_mid);
                viscosity_wz_back=jump_viscosity*wz_z_half(i,j,ij)+LEVELSET_UTILITIES<VECTOR<T,3> >::Sign(phi_mid)*jump;}
            else viscosity_wz_back=viscosity_z_half(i,j,ij)*wz_z_half(i,j,ij);
            w_viscosity(i,j,ij)=(viscosity_wx_right-viscosity_wx_left)/dx+(viscosity_wy_top-viscosity_wy_bottom)/dy
                               +(viscosity_wz_back-viscosity_wz_front)/dz;}}
}*/
//#####################################################################
// Function Discretize_Viscous_Terms ------ OLD 2D CODE!!!!!!
//#####################################################################
/*
template<class T> void INCOMPRESSIBLE_MULTIPHASE_2D<T>::
Discretize_Viscous_Terms(const ARRAY<T,VECTOR<int,2> >& phi,const T dt)
{
    int m=grid.m,n=grid.n;T dx=grid.dx,dy=grid.dy;
    T half_width=number_of_interface_cells*grid.max_dx_dy/2;

    ARRAY<T,VECTOR<int,2> > viscosity_x_half(0,m,1,n),ux_x_half(0,m,1,n),vx_x_half(0,m,1,n);
    ARRAY<T,VECTOR<int,2> > viscosity_y_half(1,m,0,n),uy_y_half(1,m,0,n),vy_y_half(1,m,0,n);
    for(int i=0;i<=m;i++) for(int j=1;j<=n;j++){
        viscosity_x_half(i,j)=LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside((T).5*(phi(i,j)+phi(i+1,j)),viscosity_minus,viscosity_plus,half_width);
        ux_x_half(i,j)=(V_ghost(i+1,j).x-V_ghost(i,j).x)/dx;
        vx_x_half(i,j)=(V_ghost(i+1,j).y-V_ghost(i,j).y)/dx;}
    for(int i=1;i<=m;i++) for(int j=0;j<=n;j++){
        viscosity_y_half(i,j)=LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside((T).5*(phi(i,j)+phi(i,j+1)),viscosity_minus,viscosity_plus,half_width);
        uy_y_half(i,j)=(V_ghost(i,j+1).x-V_ghost(i,j).x)/dy;
        vy_y_half(i,j)=(V_ghost(i,j+1).y-V_ghost(i,j).y)/dy;}

    if(!GFM){ //delta function method
        ARRAY<T,VECTOR<int,2> > uy_x_half(0,m,1,n),vx_y_half(1,m,0,n);
        {ARRAY<T,VECTOR<int,2> > uy_nodes(0,m+1,1,n);for(int i=0;i<=m+1;i++) for(int j=1;j<=n;j++) uy_nodes(i,j)=(V_ghost(i,j+1).x-V_ghost(i,j-1).x)/(2*dy);
        for(int i=0;i<=m;i++) for(int j=1;j<=n;j++) uy_x_half(i,j)=(T).5*(uy_nodes(i,j)+uy_nodes(i+1,j));}
        {ARRAY<T,VECTOR<int,2> > vx_nodes(1,m,0,n+1);for(int i=1;i<=m;i++) for(int j=0;j<=n+1;j++) vx_nodes(i,j)=(V_ghost(i+1,j).y-V_ghost(i-1,j).y)/(2*dx);
        for(int i=1;i<=m;i++) for(int j=0;j<=n;j++) vx_y_half(i,j)=(T).5*(vx_nodes(i,j)+vx_nodes(i,j+1));}
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
            u_viscosity(i,j)=2*(viscosity_x_half(i,j)*ux_x_half(i,j)-viscosity_x_half(i-1,j)*ux_x_half(i-1,j))/dx+(viscosity_y_half(i,j)*(uy_y_half(i,j)+vx_y_half(i,j))-viscosity_y_half(i,j-1)*(uy_y_half(i,j-1)+vx_y_half(i,j-1)))/dy;
            v_viscosity(i,j)=(viscosity_x_half(i,j)*(uy_x_half(i,j)+vx_x_half(i,j))-viscosity_x_half(i-1,j)*(uy_x_half(i-1,j)+vx_x_half(i-1,j)))/dx+2*(viscosity_y_half(i,j)*vy_y_half(i,j)-viscosity_y_half(i,j-1)*vy_y_half(i,j-1))/dy;}}
    else{ // GFM      
        levelset.Compute_Normals();
        ARRAY<T,VECTOR<int,2> > ux_jump(0,m+1,0,n+1),uy_jump(0,m+1,0,n+1),vx_jump(0,m+1,0,n+1),vy_jump(0,m+1,0,n+1);
        T viscosity_jump=viscosity_plus-viscosity_minus;
        for(int i=0;i<=m+1;i++) for(int j=0;j<=n+1;j++){
            MATRIX<T,2> UX((V_ghost(i+1,j)-V_ghost(i-1,j))/(2*dx),(V_ghost(i,j+1)-V_ghost(i,j-1))/(2*dy));
            VECTOR<T,2> N((*levelset.normals)(i,j)),T1(-N.y,N.x);
            SYMMETRIC_MATRIX<T,2> T_transpose_T(MATRIX<T,2>(VECTOR<T,2>(),T1).Transposed().Normal_Equations_Matrix()),N_transpose_N(SYMMETRIC_MATRIX<T,2>::Outer_Product(N));
            MATRIX<T,2> JUMP=viscosity_jump*(UX*T_transpose_T+(N_transpose_N*UX-T_transpose_T*UX.Transposed())*N_transpose_N);
            ux_jump(i,j)=JUMP.x[0];uy_jump(i,j)=JUMP.x[2];
            vx_jump(i,j)=JUMP.x[1];vy_jump(i,j)=JUMP.x[3];}

        // find pressure jump condition
        projection.poisson->Set_Jump(); // allocates the memory for projection.poisson->u_jump
        projection.poisson->levelset->Compute_Normals(); // cell center normals
        for(int i=1;i<=m-1;i++) for(int j=1;j<=n-1;j++){
            T ux=(T).5*(ux_x_half(i,j)+ux_x_half(i,j+1)),uy=(T).5*(uy_y_half(i,j)+uy_y_half(i+1,j));
            T vx=(T).5*(vx_x_half(i,j)+vx_x_half(i,j+1)),vy=(T).5*(vy_y_half(i,j)+vy_y_half(i+1,j));
            MATRIX<T,2> UX(ux,vx,uy,vy);
            VECTOR<T,2> N((*projection.poisson->levelset->normals)(i,j));
            projection.poisson->u_jump(i,j)=2*dt*viscosity_jump*VECTOR<T,2>::Dot_Product(UX*N,N);}

        // update u_viscosity = viscosity*uxx + viscosity*uyy
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
            T viscosity_ux_left,viscosity_ux_right,viscosity_uy_bottom,viscosity_uy_top;
            T phi_mid=phi(i,j),phi_left=phi(i-1,j),phi_right=phi(i+1,j),phi_bottom=phi(i,j-1),phi_top=phi(i,j+1);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_left,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_left,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_left,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_left,ux_jump(i-1,j),phi_mid,ux_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_left,phi_mid);
                viscosity_ux_left=jump_viscosity*ux_x_half(i-1,j)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_ux_left=viscosity_x_half(i-1,j)*ux_x_half(i-1,j);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_mid,phi_right)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_mid,phi_right,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_right,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_right,ux_jump(i+1,j),phi_mid,ux_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_right,phi_mid);
                viscosity_ux_right=jump_viscosity*ux_x_half(i,j)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_ux_right=viscosity_x_half(i,j)*ux_x_half(i,j);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_bottom,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_bottom,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_bottom,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_bottom,uy_jump(i,j-1),phi_mid,uy_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_bottom,phi_mid);
                viscosity_uy_bottom=jump_viscosity*uy_y_half(i,j-1)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_uy_bottom=viscosity_y_half(i,j-1)*uy_y_half(i,j-1);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_mid,phi_top)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_mid,phi_top,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_top,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_top,uy_jump(i,j+1),phi_mid,uy_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_top,phi_mid);
                viscosity_uy_top=jump_viscosity*uy_y_half(i,j)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_uy_top=viscosity_y_half(i,j)*uy_y_half(i,j);
            u_viscosity(i,j)=(viscosity_ux_right-viscosity_ux_left)/dx+(viscosity_uy_top-viscosity_uy_bottom)/dy;}

        // update v_viscosity = viscosity*vxx + viscosity*vyy
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
            T viscosity_vx_left,viscosity_vx_right,viscosity_vy_bottom,viscosity_vy_top;
            T phi_mid=phi(i,j),phi_left=phi(i-1,j),phi_right=phi(i+1,j),phi_bottom=phi(i,j-1),phi_top=phi(i,j+1);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_left,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_left,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_left,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_left,vx_jump(i-1,j),phi_mid,vx_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_left,phi_mid);
                viscosity_vx_left=jump_viscosity*vx_x_half(i-1,j)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_vx_left=viscosity_x_half(i-1,j)*vx_x_half(i-1,j);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_mid,phi_right)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_mid,phi_right,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_right,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_right,vx_jump(i+1,j),phi_mid,vx_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_right,phi_mid);
                viscosity_vx_right=jump_viscosity*vx_x_half(i,j)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_vx_right=viscosity_x_half(i,j)*vx_x_half(i,j);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_bottom,phi_mid)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_bottom,phi_mid,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_bottom,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_bottom,vy_jump(i,j-1),phi_mid,vy_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_bottom,phi_mid);
                viscosity_vy_bottom=jump_viscosity*vy_y_half(i,j-1)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_vy_bottom=viscosity_y_half(i,j-1)*vy_y_half(i,j-1);
            if(LEVELSET_UTILITIES<VECTOR<T,2> >::Interface(phi_mid,phi_top)){
                T jump_viscosity=viscosity_minus*viscosity_plus/LEVELSET_UTILITIES<VECTOR<T,2> >::Convex_Average(phi_mid,phi_top,viscosity_minus,viscosity_plus);
                T jump=jump_viscosity/LEVELSET_UTILITIES<VECTOR<T,2> >::Heaviside(phi_top,viscosity_minus,viscosity_plus)*LEVELSET_UTILITIES<VECTOR<T,2> >::Average(phi_top,vy_jump(i,j+1),phi_mid,vy_jump(i,j))*LEVELSET_UTILITIES<VECTOR<T,2> >::Theta(phi_top,phi_mid);
                viscosity_vy_top=jump_viscosity*vy_y_half(i,j)+LEVELSET_UTILITIES<VECTOR<T,2> >::Sign(phi_mid)*jump;}
            else viscosity_vy_top=viscosity_y_half(i,j)*vy_y_half(i,j);
            v_viscosity(i,j)=(viscosity_vx_right-viscosity_vx_left)/dx+(viscosity_vy_top-viscosity_vy_bottom)/dy;}}
}*/
//#####################################################################
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,1> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,1> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >;
template class INCOMPRESSIBLE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
