//#####################################################################
// Copyright 2002-2010, Mridul Aanjaneya, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_FACE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VISCOSITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_UNIFORM<T_GRID>::
INCOMPRESSIBLE_UNIFORM(const T_GRID& grid_input,PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input,THREAD_QUEUE* thread_queue_input)
    :grid(grid_input.Get_MAC_Grid()),mpi_grid(0),projection(projection_input),strain(0),collision_body_list(0),momentum_conserving_vorticity(false),use_analytic_energy(false),use_vorticity_weights(false),energy_clamp(0),vc_projection_direction(0),buoyancy_constant(0),thread_queue(thread_queue_input),
    boundary_default(*new BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>),advection_maccormack(0)
{ 
    boundary=&boundary_default;
    Initialize_Grids(grid);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_UNIFORM<T_GRID>::
~INCOMPRESSIBLE_UNIFORM()
{
    delete strain;delete &boundary_default;
    delete advection_maccormack;
}
//#####################################################################
// Function Advance_One_Time_Step_Convection
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Convection(const T dt,const T time,const T_FACE_ARRAYS_SCALAR& advecting_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect,const int number_of_ghost_cells)
{
    // TODO: make efficient if advection velocities are same as advected velocities
    // find ghost cells
    T_FACE_ARRAYS_SCALAR advection_face_velocities_ghost(grid,number_of_ghost_cells,false);
    T_FACE_ARRAYS_SCALAR face_velocities_to_advect_ghost(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,time,number_of_ghost_cells);
    boundary->Fill_Ghost_Cells_Face(grid,advecting_face_velocities,advection_face_velocities_ghost,time,number_of_ghost_cells);

    // update convection
    advection->Update_Advection_Equation_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,advection_face_velocities_ghost,*boundary,dt,time);
}

//#####################################################################
// Function Advance_One_Time_Step_Forces
// This version of the function gets as input the face_velocity_ghost.
// This is meant to work with Nimbus applications. -omidm
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T dt,const T time,const bool implicit_viscosity,const T_ARRAYS_SCALAR* phi_ghost,const int number_of_ghost_cells)
{
  // The next two lines are the only difference that we made in addition to
  // signature, instead of building the ghost velocity it will be passed -omidm

  //  T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
  //  boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);

    // update strain and apply elastic forces
    if(strain){
        assert(!projection.flame);assert(phi_ghost);strain->Update_Strain_Equation(dt,time,projection.density,face_velocities,face_velocities_ghost,*phi_ghost,number_of_ghost_cells);}

    // update gravity
    if(gravity) for(int axis=1;axis<=TV::dimension;axis++)
        DOMAIN_ITERATOR_THREADED_ALPHA<INCOMPRESSIBLE_UNIFORM<T_GRID>,TV>(grid.Face_Indices()[axis],thread_queue).template Run<T_FACE_ARRAYS_SCALAR&,const T,int>(*this,&INCOMPRESSIBLE_UNIFORM<T_GRID>::Add_Gravity_Threaded,face_velocities,dt,axis);

    // update body force
    if(use_force){
        for(int axis=1;axis<=TV::dimension;axis++)
            DOMAIN_ITERATOR_THREADED_ALPHA<INCOMPRESSIBLE_UNIFORM<T_GRID>,TV>(grid.Face_Indices()[axis],thread_queue).template Run<T_FACE_ARRAYS_SCALAR&,const T,int>(*this,&INCOMPRESSIBLE_UNIFORM<T_GRID>::Add_Body_Force_Threaded,face_velocities,dt,axis);
        boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
        boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);}

    // update viscosity explicitly
    //if(dt && (viscosity || use_variable_viscosity) && (!implicit_viscosity || use_explicit_part_of_implicit_viscosity)){
    //    if(!implicit_viscosity) PHYSBAM_NOT_IMPLEMENTED();
    //    IMPLICIT_VISCOSITY_UNIFORM<T_GRID>::Variable_Viscosity_Explicit_Part(projection.density,variable_viscosity,grid,face_velocities,face_velocities_ghost,dt,time);}
    
    T_FACE_ARRAYS_SCALAR face_velocities_old=face_velocities;
    if(vorticity_confinement || use_variable_vorticity_confinement){
        T tolerance=0; //TODO (mlentine): Look into what are good values here
        T_ARRAYS_VECTOR F(grid.Cell_Indices(1),false);
        Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
        if(collision_body_list){
            if(use_variable_vorticity_confinement){F*=dt;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement;}
        else{
            if(use_variable_vorticity_confinement){F*=dt*(T).5;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement*(T).5;}
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            for(int i=1;i<=TV::dimension;i++) if(abs(F(cell)(i))<tolerance) F(cell)(i)=0;}
        Apply_Vorticity_Confinement_Force(face_velocities,F);}
    Update_Potential_Energy(face_velocities,face_velocities_old,dt,time);

    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
}
//#####################################################################

//#####################################################################
// Function Advance_One_Time_Step_Forces
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity,const T_ARRAYS_SCALAR* phi_ghost,const int number_of_ghost_cells)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);

    // update strain and apply elastic forces
    if(strain){
        assert(!projection.flame);assert(phi_ghost);strain->Update_Strain_Equation(dt,time,projection.density,face_velocities,face_velocities_ghost,*phi_ghost,number_of_ghost_cells);}

    // update gravity
    if(gravity) for(int axis=1;axis<=TV::dimension;axis++)
        DOMAIN_ITERATOR_THREADED_ALPHA<INCOMPRESSIBLE_UNIFORM<T_GRID>,TV>(grid.Face_Indices()[axis],thread_queue).template Run<T_FACE_ARRAYS_SCALAR&,const T,int>(*this,&INCOMPRESSIBLE_UNIFORM<T_GRID>::Add_Gravity_Threaded,face_velocities,dt,axis);

    // update body force
    if(use_force){
        for(int axis=1;axis<=TV::dimension;axis++)
            DOMAIN_ITERATOR_THREADED_ALPHA<INCOMPRESSIBLE_UNIFORM<T_GRID>,TV>(grid.Face_Indices()[axis],thread_queue).template Run<T_FACE_ARRAYS_SCALAR&,const T,int>(*this,&INCOMPRESSIBLE_UNIFORM<T_GRID>::Add_Body_Force_Threaded,face_velocities,dt,axis);
        boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
        boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);}

    // update viscosity explicitly
    //if(dt && (viscosity || use_variable_viscosity) && (!implicit_viscosity || use_explicit_part_of_implicit_viscosity)){
    //    if(!implicit_viscosity) PHYSBAM_NOT_IMPLEMENTED();
    //    IMPLICIT_VISCOSITY_UNIFORM<T_GRID>::Variable_Viscosity_Explicit_Part(projection.density,variable_viscosity,grid,face_velocities,face_velocities_ghost,dt,time);}
    
    T_FACE_ARRAYS_SCALAR face_velocities_old=face_velocities;
    if(vorticity_confinement || use_variable_vorticity_confinement){
        T tolerance=0; //TODO (mlentine): Look into what are good values here
        T_ARRAYS_VECTOR F(grid.Cell_Indices(1),false);
        Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
        if(collision_body_list){
            if(use_variable_vorticity_confinement){F*=dt;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement;}
        else{
            if(use_variable_vorticity_confinement){F*=dt*(T).5;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement*(T).5;}
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            for(int i=1;i<=TV::dimension;i++) if(abs(F(cell)(i))<tolerance) F(cell)(i)=0;}
        Apply_Vorticity_Confinement_Force(face_velocities,F);}
    Update_Potential_Energy(face_velocities,face_velocities_old,dt,time);

    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Add_Gravity_Threaded(RANGE<TV_INT>& domain,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,int axis)
{
    for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
        if(downward_direction[axis]) face_velocities.Component(axis)(iterator.Face_Index())+=dt*gravity*downward_direction[axis];}
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Add_Body_Force_Threaded(RANGE<TV_INT>& domain,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,int axis)
{
    for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
        face_velocities.Component(iterator.Axis())(iterator.Face_Index())+=dt*force.Component(iterator.Axis())(iterator.Face_Index());}
}
//#####################################################################
// Function Calculate_Kinetic_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Calculate_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(conserve_kinetic_energy) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        kinetic_energy(face)=(T).5*projection.density*face_velocities(face)*face_velocities(face);}
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Apply_Pressure_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,const T dt,const T time)
{
    T_FACE_ARRAYS_SCALAR face_velocities(grid);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();face_velocities(face)=(face_velocities_new(face)+face_velocities_old(face))/(T)2.;}
    if(conserve_kinetic_energy) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        if(!projection.elliptic_solver->psi_N(face) && !(projection.elliptic_solver->psi_D(first_cell) && projection.elliptic_solver->psi_D(second_cell)))
            kinetic_energy(face)-=face_velocities(face)*(projection.p(second_cell)-projection.p(first_cell))/grid.dX[iterator.Axis()];}
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Update_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,const T dt,const T time)
{
    if(conserve_kinetic_energy) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        kinetic_energy(face)+=(T).5*projection.density*(face_velocities_new(face)*face_velocities_new(face)-face_velocities_old(face)*face_velocities_old(face));}
}
//#####################################################################
// Function Update_Potential_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(conserve_kinetic_energy) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();if(projection.elliptic_solver->psi_N(face)) continue;
        potential_energy(face)+=kinetic_energy(face)-(T).5*projection.density*(face_velocities(face)*face_velocities(face));}
}
//#####################################################################
// Function Correct_Kinetic_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,const T dt,const T time)
{
    int number_of_faces=0;
    if(conserve_kinetic_energy) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) number_of_faces++;
    if(conserve_kinetic_energy) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        potential_energy(face)-=(T).5*projection.density*(face_velocities_new(face)*face_velocities_new(face)-face_velocities_old(face)*face_velocities_old(face));
        potential_energy(face)+=allowed_energy_gained/(T)number_of_faces;}
    allowed_energy_gained=0;
}
//#####################################################################
// Function Correct_Kinetic_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Correct_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(conserve_kinetic_energy) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(kinetic_energy(face)<0) kinetic_energy(face)=0;
        face_velocities(face)=sign_nonzero(face_velocities(face))*sqrt(kinetic_energy(face)*(T)2./projection.density);}
}
//#####################################################################
// Function Correct_Kinetic_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Advect_With_Vorticity(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const int iterations,const int number_of_ghost_cells)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Advect with Vorticity START",0,0);
    //calculate vorticity directions
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);
    T_ARRAYS_SPIN vorticity(grid.Cell_Indices(2),false);
    T_ARRAYS_SCALAR vorticity_magnitude(grid.Cell_Indices(2));
    T_ARRAYS_VECTOR vortex_normal_vector(grid.Cell_Indices(1));
    if(collision_body_list){
        FACE_LOOKUP_UNIFORM<T_GRID> face_velocities_lookup_uniform(face_velocities_ghost);
        FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> face_velocities_lookup(face_velocities_lookup_uniform,*collision_body_list,&valid_mask);
        VORTICITY_UNIFORM<TV>::Vorticity(grid,face_velocities_lookup,vorticity,vorticity_magnitude);}
    else VORTICITY_UNIFORM<TV>::Vorticity(grid,T_FACE_LOOKUP(face_velocities_ghost),vorticity,vorticity_magnitude);
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){ // do collision awareness when these are averaged to faces
        vortex_normal_vector(iterator.Cell_Index())=T_LEVELSET::Normal_At_Node(grid,vorticity_magnitude,iterator.Cell_Index());}
    T_FACE_ARRAYS_SCALAR face_normal_vectors(grid);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();int axis=full_index.axis;
        face_normal_vectors(full_index)=(vortex_normal_vector(iterator.First_Cell_Index())(axis)+vortex_normal_vector(iterator.Second_Cell_Index())(axis))/(T)2.;}

    //Advect energy using vorticity vectors
    for(int i=1;i<=iterations;i++){
        boundary->Set_Fixed_Boundary(true,0);
        ARRAY<T,FACE_INDEX<TV::dimension> > kinetic_energy_ghost(grid,7),face_normals_ghost(grid,4);
        boundary->Fill_Ghost_Cells_Face(grid,kinetic_energy,kinetic_energy_ghost,time,7);
        boundary->Fill_Ghost_Cells_Face(grid,face_normal_vectors,face_normals_ghost,time,4);
        advection->Update_Advection_Equation_Face(grid,kinetic_energy,kinetic_energy_ghost,face_normals_ghost,*boundary,dt,time);
        ARRAY<T,FACE_INDEX<TV::dimension> > potential_energy_ghost(grid,7);
        boundary->Fill_Ghost_Cells_Face(grid,potential_energy,potential_energy_ghost,time,7);
        advection->Update_Advection_Equation_Face(grid,potential_energy,potential_energy_ghost,face_normals_ghost,*boundary,dt,time);
        boundary->Set_Fixed_Boundary(false);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Advect with Vorticity END",0,0);}
}
//#####################################################################
// Function Correct_Kinetic_Energy
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Add_Energy_With_Vorticity(T_FACE_ARRAYS_SCALAR& face_velocities,const VECTOR<VECTOR<bool,2>,TV::dimension>& domain_boundary,const T dt,const T time,const int number_of_ghost_cells,T_LEVELSET* lsv,T_ARRAYS_SCALAR* density)
{
    use_vorticity_weights=true;
    vorticity_weights.Resize(grid);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> index=iterator.Full_Index();
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(iterator.Axis())++;
        for(int i=1;i<=TV::dimension;i++){if(domain_boundary(i)(1)) domain.min_corner(i)++;if(domain_boundary(i)(2)) domain.max_corner(i)--;}
        vorticity_weights(index)=(conserve_kinetic_energy && lsv)?(-1*lsv->Phi(iterator.Location())):1;
        if(vorticity_weights(index)<energy_clamp) vorticity_weights(index)=0;if(projection.elliptic_solver->psi_N(index) || !domain.Lazy_Inside(index.index)) vorticity_weights(index)=0;}
    if(!conserve_kinetic_energy) return;
    //int iterations=1;
    //Advect_With_Vorticity(face_velocities,dt,time,iterations);

    T_ARRAYS_VECTOR F(grid.Cell_Indices(1));
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    T_FACE_ARRAYS_SCALAR potential_energy_ghost;potential_energy_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);
    boundary->Fill_Ghost_Cells_Face(grid,potential_energy,potential_energy_ghost,time,number_of_ghost_cells);
    Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
    if(!use_variable_vorticity_confinement) Use_Variable_Vorticity_Confinement(true);
    if(true){
        std::stringstream ss;
        T uF=0,PE=0,Fsqr=0;
        if(use_analytic_energy){
            T KE=0,potential=0;
            for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                if(lsv && lsv->phi(iterator.Cell_Index())>0) continue;
                if(density) potential+=(*density)(iterator.Cell_Index())*(*density)(iterator.Cell_Index())*buoyancy_constant*(1-iterator.Location()(2));
                else potential+=projection.density*grid.dX.Product()*iterator.Location()(2)*(T)9.8;}
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> index=iterator.Full_Index();if(lsv && lsv->Phi(iterator.Location())>0) continue;
                T mass=density?(((*density)(iterator.First_Cell_Index())+((*density)(iterator.Second_Cell_Index())))*(T).5):projection.density*grid.dX.Product();
                KE+=(T)0.5*mass*face_velocities(index)*face_velocities(index);
                potential_energy(index)=0;}
            ss<<"Confinement analytical "<<analytic_energy<<" potential "<<potential<<" KE "<<KE<<std::endl;
            PE=(analytic_energy-KE)+(analytic_potential-potential);
            if(!density) PE/=grid.dX.Product();}
        else for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> index=iterator.Full_Index();
            PE+=potential_energy(index);}
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> index=iterator.Full_Index();
            T mass=density?(((*density)(iterator.First_Cell_Index())+((*density)(iterator.Second_Cell_Index())))*(T).5):projection.density;
            T force=vorticity_weights(index)*(F(iterator.First_Cell_Index())(iterator.Axis())+F(iterator.Second_Cell_Index())(iterator.Axis()))*dt/(T)2.;
            Fsqr+=mass*force*force;uF+=mass*face_velocities(index)*force;}
        ss<<"Confinement debug uF "<<uF<<" Fsqr "<<Fsqr<<" PE "<<PE<<" allowed "<<allowed_energy_gained<<std::endl;
        if(mpi_grid){
            T allowed_mpi=mpi_grid->Reduce_Add(allowed_energy_gained),PE_mpi=mpi_grid->Reduce_Add(PE),uF_mpi=mpi_grid->Reduce_Add(uF),Fsqr_mpi=mpi_grid->Reduce_Add(Fsqr);
            T sqr_term=uF_mpi*uF_mpi+2*(PE_mpi+allowed_mpi)*Fsqr_mpi;if(sqr_term<0) sqr_term=0;
            vorticity_confinement=Fsqr_mpi>1e-2?(sign(uF_mpi)*(sqrt(sqr_term)-abs(uF_mpi))/(Fsqr_mpi)):0;}
        else{
            T sqr_term=uF*uF+2*(PE+allowed_energy_gained)*Fsqr;if(sqr_term<0) sqr_term=0;
            vorticity_confinement=Fsqr>1e-2?(sign(uF)*(sqrt(sqr_term)-abs(uF))/(Fsqr)):0;}
        ss<<"Confinement force "<<vorticity_confinement<<std::endl;LOG::filecout(ss.str());
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
            variable_vorticity_confinement(index)=vorticity_confinement;}}
    else{
        potential_energy_ghost*=(T)0.5;
        T_ARRAYS_VECTOR V(grid.Cell_Indices(1),false);
        T_ARRAYS_SCALAR PE(grid.Domain_Indices(1),false);
        AVERAGING_UNIFORM<T_GRID> averaging;
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            V(cell)=averaging.Face_To_Cell_Vector(grid,cell,face_velocities_ghost);
            PE(cell)=0;for(int axis=1;axis<=TV::dimension;axis++) PE(cell)+=potential_energy_ghost(iterator.Full_First_Face_Index(axis))+potential_energy_ghost(iterator.Full_Second_Face_Index(axis));}
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            T energy_sqr=2*(PE(cell)+(T)0.5*projection.density*V(cell).Magnitude_Squared())/projection.density;
            if(energy_sqr<0) energy_sqr=0;T energy=sqrt(energy_sqr)-V(cell).Magnitude();
            T mag=F(cell).Magnitude();
            if(mag<1e-6) variable_vorticity_confinement(cell)=0;
            else variable_vorticity_confinement(cell)=energy/mag;
            if(abs(variable_vorticity_confinement(cell))>vorticity_confinement*mag) variable_vorticity_confinement(iterator.Cell_Index())=vorticity_confinement*mag*sign(variable_vorticity_confinement(iterator.Cell_Index()));}}
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Advance_One_Time_Step_Implicit_Part(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity,BOUNDARY_UNIFORM<T_GRID,T>* projection_boundary,
    bool use_levelset_viscosity,BOUNDARY_CONDITIONS_CALLBACKS<TV>* bc_callbacks,bool print_viscosity_matrix)
{
    int ghost_cells=3;
    // boundary conditions
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before apply boundary",0,0);
    if(!projection_boundary) projection_boundary=boundary;
    projection_boundary->Apply_Boundary_Condition_Face(projection.p_grid,face_velocities,time+dt);
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(projection.p_grid,ghost_cells);
    projection_boundary->Fill_Ghost_Cells_Face(projection.p_grid,face_velocities,face_velocities_ghost,time,ghost_cells);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After apply boundary",0,0);

    assert(Consistent_Boundary_Conditions(face_velocities));

    // viscosity
    if(use_levelset_viscosity && implicit_viscosity && dt && viscosity!=0){
        LEVELSET_VISCOSITY_UNIFORM<TV> levelset_viscosity(bc_callbacks,grid,dt,projection.density,viscosity);
        levelset_viscosity.print_matrix=print_viscosity_matrix;
        levelset_viscosity.Apply_Viscosity(face_velocities,false,true,false);}
    else if(implicit_viscosity && dt && (use_variable_viscosity || viscosity!=0)){
        projection.Make_Divergence_Free(face_velocities,dt,time);
        Implicit_Viscous_Update(face_velocities,dt,time);}
        //VISCOSITY<T_GRID> viscosity_helper(*projection.elliptic_solver,variable_viscosity,projection.density,viscosity,implicit_viscosity,use_explicit_part_of_implicit_viscosity,use_variable_viscosity,maximum_implicit_viscosity_iterations);
        //viscosity_helper.Add_Implicit_Forces_Before_Projection(grid,face_velocities,face_velocities,dt,time);}

    projection.Make_Divergence_Free(face_velocities,dt,time);
    if(conserve_kinetic_energy) projection.Update_Potential_Energy(face_velocities_ghost,face_velocities,potential_energy,allowed_energy_gained,dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After final projection",0,0);
}
//#####################################################################
// Function Implicit_Viscous_Update
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Implicit_Viscous_Update(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{ 
    int number_of_ghost_cells=3;
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);

    for(int axis=1;axis<=T_GRID::dimension;axis++){
        IMPLICIT_VISCOSITY_UNIFORM<T_GRID> implicit_viscosity(*projection.elliptic_solver,variable_viscosity,projection.density,viscosity,0,axis,false,false);
        implicit_viscosity.Viscous_Update(grid,face_velocities,face_velocities_ghost,dt,time,maximum_implicit_viscosity_iterations);}
    if(mpi_grid) mpi_grid->Copy_Common_Face_Data(face_velocities);
    Update_Kinetic_Energy(face_velocities,face_velocities_ghost,dt,time);
}
//#####################################################################
// Function Real_CFL
//#####################################################################
template<class T_GRID> int INCOMPRESSIBLE_UNIFORM<T_GRID>::
Real_CFL(T_FACE_ARRAYS_SCALAR& face_velocities,const bool inviscid,const bool viscous_only,T input_dt) const
{
    T dt=CFL(face_velocities,inviscid,viscous_only);
    return (int) (input_dt/dt + 1);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR INCOMPRESSIBLE_UNIFORM<T_GRID>::
CFL(T_FACE_ARRAYS_SCALAR& face_velocities,const bool inviscid,const bool viscous_only) const
{
    TV DX=grid.dX,sqr_DX=DX*DX,max_abs_V,one_over_DX=grid.one_over_dX;
    int ghost_cells=3;
    // convection
    T dt_convection=0;
    if(!viscous_only){
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            T local_V_norm=0;for(int axis=1;axis<=T_GRID::dimension;axis++)
                local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),
                    face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
            dt_convection=max(dt_convection,local_V_norm);}}
    // surface tension
    T dt_surface_tension=0;
    if(nonzero_surface_tension){
        T_LEVELSET& levelset=*projection.collidable_solver->levelset;levelset.Compute_Curvature();
        T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(ghost_cells));levelset.boundary->Fill_Ghost_Cells(grid,levelset.phi,phi_ghost,0,0,ghost_cells); // TODO: use real dt, time
        T kappa_cfl=0;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            T phi_1=phi_ghost(cell_1);T phi_2=phi_ghost(cell_2);
            if(LEVELSET_UTILITIES<T>::Interface(phi_1,phi_2)){
                T surface_tension_coefficient=surface_tension;
                if(use_variable_surface_tension) surface_tension_coefficient=(T).5*(variable_surface_tension(cell_1)+variable_surface_tension(cell_2));
                if(surface_tension_coefficient){
                    T curvature=LEVELSET_UTILITIES<T>::Average(phi_1,(*levelset.curvature)(cell_1),phi_2,(*levelset.curvature)(cell_2));
                    kappa_cfl=max(kappa_cfl,abs(curvature*surface_tension_coefficient/projection.density));}}}
        dt_surface_tension=sqrt(kappa_cfl)/grid.Minimum_Edge_Length();}
    T dt_viscosity=0;
    /*if(!inviscid && nonzero_viscosity){
        T norm_2_over_sqr_DX=2*Inverse(sqr_DX).L1_Norm();
        if(viscosity) dt_viscosity=viscosity/projection.density*norm_2_over_sqr_DX;
        if(use_variable_viscosity){assert(!projection.flame);dt_viscosity=variable_viscosity.Maxabs()/projection.density*norm_2_over_sqr_DX;}}*/
    if(viscous_only) return 1/max(dt_viscosity,1/max_time_step);
    TV max_force;
    if(use_force) max_force=force.Maxabs();
    T dt_force=0;
    if(use_force) dt_force=(max_force*one_over_DX).L1_Norm();
    if(gravity) dt_force+=abs(gravity)*(downward_direction*one_over_DX).L1_Norm();
    if(strain) dt_force+=1/strain->CFL(projection.density);
    T dt_overall=(dt_convection+dt_viscosity+sqrt(sqr(dt_convection+dt_viscosity)+4*dt_force+4*sqr(dt_surface_tension)))/2; 
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Extrapolate_Velocity_Across_Interface(T_FACE_ARRAYS_SCALAR& face_velocities,T_ARRAYS_SCALAR& phi_ghost,const bool enforce_divergence_free,const T band_width,
    const T damping,const TV& air_speed,const T_FACE_ARRAYS_BOOL_DIMENSION* face_neighbors_visible,const T_FACE_ARRAYS_BOOL* fixed_faces_input)
{
    T extrapolation_band_width=enforce_divergence_free?band_width+2:band_width;
    T delta=extrapolation_band_width*grid.dX.Max();
    const int ghost_cells=2*(int)ceilf((float)extrapolation_band_width)+1;
    if(mpi_grid){
        T_FACE_ARRAYS_SCALAR phi_faces(grid,ghost_cells);
        T_FACE_ARRAYS_SCALAR face_velocities_ghost(grid,ghost_cells,false);boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,0,ghost_cells); // TODO: use real time
        T_FACE_ARRAYS_BOOL fixed_faces=fixed_faces_input?*fixed_faces_input:T_FACE_ARRAYS_BOOL(grid,ghost_cells);
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            T_ARRAYS_BASE &phi_face=phi_faces.Component(axis),&face_velocity=face_velocities_ghost.Component(axis);T_ARRAYS_BOOL_BASE& fixed_face=fixed_faces.Component(axis);
            for(FACE_ITERATOR iterator(grid,ghost_cells,T_GRID::INTERIOR_REGION,0,axis);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Face_Index();phi_face(index)=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
                if(phi_face(index)<=0) fixed_face(index)=true;if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0;}}
        //mpi_grid->Exchange_Boundary_Face_Data(fixed_faces);
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            T_ARRAYS_BASE &phi_face=phi_faces.Component(axis),&face_velocity=face_velocities_ghost.Component(axis);T_ARRAYS_BOOL_BASE& fixed_face=fixed_faces.Component(axis);
            T_GRID face_grid=grid.Get_Face_Grid(axis);T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,ghost_cells);
            extrapolate.Set_Band_Width(extrapolation_band_width);extrapolate.Set_Custom_Seed_Done(&fixed_face);
            if(face_neighbors_visible) extrapolate.Set_Collision_Aware_Extrapolation(face_neighbors_visible->Component(axis));
            extrapolate.Extrapolate(0,false);
            T_ARRAYS_SCALAR::Get(face_velocities.Component(axis),face_velocity);
            if(damping) for(FACE_ITERATOR iterator(grid,0,T_GRID::INTERIOR_REGION,ghost_cells,axis);iterator.Valid();iterator.Next()){TV_INT index=iterator.Face_Index();
                if(!fixed_face(index) && phi_face(index)<delta) face_velocity(index)=(1-damping)*face_velocity(index)+damping*air_speed[axis];}}}
    else{
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            T_GRID face_grid=grid.Get_Face_Grid(axis);T_ARRAYS_SCALAR phi_face(face_grid.Domain_Indices(),false);T_ARRAYS_BASE& face_velocity=face_velocities.Component(axis);
            T_ARRAYS_BOOL fixed_face;
            if(fixed_faces_input) fixed_face=fixed_faces_input->Component(axis); else fixed_face=T_ARRAYS_BOOL(face_grid.Domain_Indices());
            for(FACE_ITERATOR iterator(grid,0,T_GRID::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Face_Index();phi_face(index)=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
                if(phi_face(index)<=0) fixed_face(index)=true;if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0;}
            std::stringstream ss;ss<<"something..."<<std::endl;LOG::filecout(ss.str());  // TODO(jontg): If this log statement doesn't appear, the code crashes in release mode...
            T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,ghost_cells);extrapolate.Set_Band_Width(extrapolation_band_width);extrapolate.Set_Custom_Seed_Done(&fixed_face);
            if(face_neighbors_visible) extrapolate.Set_Collision_Aware_Extrapolation(face_neighbors_visible->Component(axis));
            extrapolate.Extrapolate();
            if(damping) for(FACE_ITERATOR iterator(grid,0,T_GRID::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next()){TV_INT index=iterator.Face_Index();
                if(!fixed_face(index) && phi_face(index)<delta) face_velocity(index)=(1-damping)*face_velocity(index)+damping*air_speed[axis];}}}

    // make extrapolated velocity divergence free
    if(enforce_divergence_free){
        extrapolation_band_width-=2;delta=extrapolation_band_width*grid.dX.Max();
        T_ARRAYS_SCALAR p_new(grid.Domain_Indices(1));T_ARRAYS_BOOL psi_D_new(grid.Domain_Indices(1));T_FACE_ARRAYS_BOOL psi_N_new(grid,1);
        T_ARRAYS_SCALAR::Exchange_Arrays(p_new,projection.p);T_ARRAYS_BOOL::Exchange_Arrays(psi_D_new,projection.elliptic_solver->psi_D);
        T_FACE_ARRAYS_BOOL::Exchange_Arrays(psi_N_new,projection.elliptic_solver->psi_N);
        VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary_local;
        for(int i=1;i<=TV::dimension;i++){domain_boundary_local(i)(1)=true;domain_boundary_local(i)(2)=true;}        
        if(mpi_grid) mpi_grid->Initialize(domain_boundary_local);
        for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
            TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
            TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
            if(domain_boundary_local(axis)(axis_side)){
                for(typename GRID<TV>::FACE_ITERATOR iterator(grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT face=iterator.Face_Index()+boundary_face_offset;
                    if(phi_ghost(face+interior_cell_offset)<=0){
                        if(face_velocities.Component(axis).Valid_Index(face)){projection.elliptic_solver->psi_N.Component(axis)(face)=true;face_velocities.Component(axis)(face)=0;}}
                    else{TV_INT cell=face+exterior_cell_offset;projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}
            else for(typename GRID<TV>::FACE_ITERATOR iterator(grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            phi_ghost(iterator.Cell_Index())-=delta;}
        Set_Dirichlet_Boundary_Conditions(&phi_ghost,0);
        //if(projection.elliptic_solver->mpi_grid){
        //    projection.elliptic_solver->mpi_grid->Exchange_Boundary_Face_Data(projection.elliptic_solver->psi_N,1);}
        //projection.elliptic_solver->Set_Dirichlet_Outer_Boundaries();
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            phi_ghost(iterator.Cell_Index())+=delta;}
        /*for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
            //if(phi_ghost(index) >= delta) projection.elliptic_solver->psi_D(index)=true;
            //else if(phi_ghost(index) <= 0) projection.elliptic_solver->psi_N.Set_All_Faces(true,index);
            if(phi_ghost(index) <= 0) projection.elliptic_solver->psi_N.Set_All_Faces(true,index);
            //else;
            {
                bool local_maximum=true;
                for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++) if(phi_ghost(index)<phi_ghost(iterator.Cell_Neighbor(i))){local_maximum=false;break;}
                if(local_maximum)projection.elliptic_solver->psi_D(index)=true;}}*/
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            if(phi_ghost(iterator.First_Cell_Index())<=0 || phi_ghost(iterator.Second_Cell_Index())<=0) projection.elliptic_solver->psi_N(iterator.Full_Index())=true;}
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            phi_ghost(iterator.Cell_Index())-=delta;}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("After setting boundary conditions in divergence free solve.",0,0);
        projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(true);

        /*projection.elliptic_solver->Find_Solution_Regions();
        T_ARRAYS_SCALAR p_back=projection.p;
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            projection.p(iterator.Cell_Index())=projection.elliptic_solver->filled_region_touches_dirichlet(projection.elliptic_solver->filled_region_colors(iterator.Cell_Index()))?1:0;
        }
        projection.p=p_back;*/
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Viewing touches dirichlet",0,0);
        
        Advance_One_Time_Step_Implicit_Part(face_velocities,0,0);
        //projection.Make_Divergence_Free(face_velocities,0,0); // TODO: use real dt/time
        projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            phi_ghost(iterator.Cell_Index())+=delta;}
        T_ARRAYS_SCALAR::Exchange_Arrays(p_new,projection.p);T_ARRAYS_BOOL::Exchange_Arrays(psi_D_new,projection.elliptic_solver->psi_D);
        T_FACE_ARRAYS_BOOL::Exchange_Arrays(psi_N_new,projection.elliptic_solver->psi_N);} // restore pressure for use as initial guess for incompressible projection
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T_ARRAYS_SCALAR* phi,const T pressure)
{
    T_ARRAYS_BOOL& psi_D=projection.elliptic_solver->psi_D;
    if(phi) for(CELL_ITERATOR iterator(projection.p_grid);iterator.Valid();iterator.Next()) if((*phi)(iterator.Cell_Index())>0){
        psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=pressure;}
    if(projection.elliptic_solver->mpi_grid){
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(psi_D,1,false);
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After exchange dirichlet",0,0);
}
//#####################################################################
// Function Set_Boundary_Conditions_With_Extrapolation
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Set_Boundary_Conditions_With_Extrapolation(const T_ARRAYS_SCALAR* phi,const T pressure,const T_ARRAYS_BOOL &extrapolated_psi_D)
{
    T_ARRAYS_BOOL& psi_D=projection.elliptic_solver->psi_D;
    if(phi){for(CELL_ITERATOR iterator(projection.p_grid);iterator.Valid();iterator.Next()){
            if((*phi)(iterator.Cell_Index())>0&&(!extrapolated_psi_D(iterator.Cell_Index()))){
                psi_D(iterator.Cell_Index())=true;
                projection.p(iterator.Cell_Index())=pressure;}}}
    if(projection.elliptic_solver->mpi_grid){
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(psi_D,1,false);
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T_ARRAYS_SCALAR* phi,const T_ARRAYS_SCALAR& pressure)
{
    T_ARRAYS_BOOL& psi_D=projection.elliptic_solver->psi_D;
    if(phi) for(CELL_ITERATOR iterator(projection.p_grid);iterator.Valid();iterator.Next()) if((*phi)(iterator.Cell_Index())>0){
        psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=pressure(iterator.Cell_Index());}
    if(mpi_grid){
        mpi_grid->Exchange_Boundary_Cell_Data(psi_D,1,false);
        mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After exchange dirichlet",0,0);
}
//#####################################################################
// Function Add_Surface_Tension
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Add_Surface_Tension(T_LEVELSET& levelset,const T time)
{
    typedef typename T_GRID::VECTOR_T TV;
    LAPLACE_UNIFORM<T_GRID>& elliptic_solver=*projection.elliptic_solver;T_GRID& p_grid=elliptic_solver.grid;
    T_ARRAYS_SCALAR& phi=levelset.phi;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;

    if(projection.collidable_solver->second_order_cut_cell_method) for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
        if(!projection.elliptic_solver->psi_N.Component(iterator.Axis())(iterator.Face_Index()) && LEVELSET_UTILITIES<T>::Interface(phi(iterator.First_Cell_Index()),phi(iterator.Second_Cell_Index()))){
            T theta=LEVELSET_UTILITIES<T>::Theta(phi(iterator.First_Cell_Index()),phi(iterator.Second_Cell_Index()));
            TV location=theta*(grid.Center(iterator.Second_Cell_Index())-grid.Center(iterator.First_Cell_Index()))+grid.Center(iterator.First_Cell_Index());
            T curvature_at_interface=levelset.Compute_Curvature(location);
            T surface_tension_coefficient=surface_tension;
            if(use_variable_surface_tension)surface_tension_coefficient=interpolation.Clamped_To_Array(grid,variable_surface_tension,location);
            projection.collidable_solver->u_interface.Component(iterator.Axis())(iterator.Face_Index())=-surface_tension_coefficient*curvature_at_interface;}}
    else{
        for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()) if(elliptic_solver.psi_D(iterator.Cell_Index()) && phi(iterator.Cell_Index()) < 5*grid.dX.Max()){
            T surface_tension_coefficient=surface_tension;
            if(use_variable_surface_tension) surface_tension_coefficient=variable_surface_tension(iterator.Cell_Index());
            projection.p(iterator.Cell_Index())=-surface_tension_coefficient*levelset.Compute_Curvature(phi,iterator.Cell_Index());}}
}
//#####################################################################
// Function Apply_Vorticity_Confinement_Force
//#####################################################################
template<class T_GRID,class T_ARRAYS_TV> static void
Apply_Vorticity_Confinement_Force_Helper(const T_GRID& grid,typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& face_velocities,T_ARRAYS_TV& F,const INCOMPRESSIBLE_UNIFORM<T_GRID>& incompressible)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    // want cells to face averaging here
    if(incompressible.collision_body_list){
        AVERAGING_COLLIDABLE_UNIFORM<T_GRID,FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> > vorticity_averaging_collidable(*incompressible.collision_body_list,T());
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities.Component(axis)(iterator.Face_Index())+=vorticity_averaging_collidable.Cell_To_Face(grid,axis,iterator.Face_Index(),F);}}
    else
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
            face_velocities.Component(axis)(iterator.Face_Index())+=(incompressible.use_vorticity_weights?incompressible.vorticity_weights(iterator.Full_Index()):1)*(F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis]);}
}
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<float,1> >&,ARRAY<float,FACE_INDEX<1> >&,ARRAY<VECTOR<float,1> ,VECTOR<int,1> >&,const INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<float,1> > >&)
{PHYSBAM_NOT_IMPLEMENTED();}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<double,1> >&,ARRAY<double,FACE_INDEX<1> >&,ARRAY<VECTOR<double,1> ,VECTOR<int,1> >&,const INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<double,1> > >&)
{PHYSBAM_NOT_IMPLEMENTED();}
#endif
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Apply_Vorticity_Confinement_Force(T_FACE_ARRAYS_SCALAR& face_velocities,T_ARRAYS_VECTOR& F)
{
    Apply_Vorticity_Confinement_Force_Helper(grid,face_velocities,F,*this);
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force_Helper
//#####################################################################
template<class T_GRID,class T_FACE_ARRAYS,class T_ARRAYS_TV> static void 
Compute_Vorticity_Confinement_Force_Helper(const T_GRID& grid,const T_FACE_ARRAYS& face_velocities_ghost,T_ARRAYS_TV& F,const INCOMPRESSIBLE_UNIFORM<T_GRID>& incompressible)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_ARRAYS_TV::template REBIND<T>::TYPE T_ARRAYS_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_ARRAYS_TV::template REBIND<typename TV::SPIN>::TYPE T_ARRAYS_SPIN;typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    T_ARRAYS_SPIN vorticity(grid.Cell_Indices(2),false);
    T_ARRAYS_SCALAR vorticity_magnitude(grid.Cell_Indices(2));
    if(incompressible.collision_body_list){
        FACE_LOOKUP_UNIFORM<T_GRID> face_velocities_lookup_uniform(face_velocities_ghost);
        FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> face_velocities_lookup(face_velocities_lookup_uniform,*incompressible.collision_body_list,&incompressible.valid_mask);
        VORTICITY_UNIFORM<TV>::Vorticity(grid,face_velocities_lookup,vorticity,vorticity_magnitude);}
    else VORTICITY_UNIFORM<TV>::Vorticity(grid,T_FACE_LOOKUP(face_velocities_ghost),vorticity,vorticity_magnitude);
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){ // do collision awareness when these are averaged to faces
        TV vortex_normal_vector=T_LEVELSET::Normal_At_Node(grid,vorticity_magnitude,iterator.Cell_Index());
        F(iterator.Cell_Index())=TV::Cross_Product(vortex_normal_vector,vorticity(iterator.Cell_Index()));}

    if(incompressible.vc_projection_direction){
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) F(iterator.Cell_Index())(incompressible.vc_projection_direction)=0;}
    if(incompressible.momentum_conserving_vorticity){
        for(int axis=1;axis<=TV::dimension;axis++){T sum=0; int count=0;
            for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                if(incompressible.variable_vorticity_confinement(cell)!=(T)0){sum+=F(cell)[axis];count++;}}
            PHYSBAM_ASSERT(count>0);
            sum/=count;T final_sum=0;
            for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index(); 
                if(incompressible.variable_vorticity_confinement(cell)!=(T)0){F(cell)[axis]-=sum;final_sum+=F(cell)[axis];}}}}
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force_Helper
//#####################################################################
template<class T,class T_FACE_ARRAYS,class T_ARRAYS_TV> static void
Compute_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<T,1> >&,const T_FACE_ARRAYS&,T_ARRAYS_TV&,const INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<T,1> > >&)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Compute_Vorticity_Confinement_Force(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_ARRAYS_VECTOR& F)
{
    Compute_Vorticity_Confinement_Force_Helper(grid,face_velocities_ghost,F,*this);
}
//#####################################################################
// Function Use_Maccormack_Advection
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Use_Maccormack_Advection(const T_ARRAYS_BOOL* node_mask, const T_ARRAYS_BOOL* cell_mask, const T_FACE_ARRAYS_BOOL* face_mask)
{
    advection_maccormack=new ADVECTION_MACCORMACK_UNIFORM<T_GRID,T,ADVECTION<T_GRID,T> >(*advection,node_mask,cell_mask,face_mask,thread_queue);
    Set_Custom_Advection(*advection_maccormack);
}
//#####################################################################
// Function Consistent_Boundary_Conditions
//#####################################################################
template<class T_GRID> bool INCOMPRESSIBLE_UNIFORM<T_GRID>::
Consistent_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& face_velocities) const
{
    if(projection.elliptic_solver->mpi_grid){
        std::stringstream ss;ss<<"checking for consistent mpi boundaries"<<std::endl;LOG::filecout(ss.str());
        T_ARRAYS_BOOL psi_D_ghost(projection.elliptic_solver->psi_D);
        T_FACE_ARRAYS_SCALAR face_velocities_ghost(face_velocities);T_FACE_ARRAYS_SCALAR psi_N_ghost(projection.p_grid);
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(psi_D_ghost,1);
        for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
            if(projection.elliptic_solver->mpi_grid->Neighbor(axis,axis_side)){TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
                for(FACE_ITERATOR iterator(projection.p_grid,0,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT face=iterator.Face_Index(),cell=face+exterior_cell_offset;int axis=iterator.Axis();
                    psi_N_ghost(axis,face)=(T)projection.elliptic_solver->psi_N(axis,face);
                    face_velocities_ghost(axis,face)=face_velocities(axis,face);
                    assert(projection.elliptic_solver->psi_D(cell)==psi_D_ghost(cell));}}}
        projection.elliptic_solver->mpi_grid->Assert_Common_Face_Data(face_velocities_ghost);projection.elliptic_solver->mpi_grid->Assert_Common_Face_Data(psi_N_ghost);}
    return true;
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Initialize_Grids(const T_GRID& grid_input)
{
    INCOMPRESSIBLE<T_GRID>::Initialize_Grids(grid_input);
    grid=grid_input.Get_MAC_Grid();projection.Initialize_Grid(grid);
    if(strain) strain->Initialize_Grid(grid);
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Set_Body_Force(const bool use_force_input)
{
    use_force=use_force_input;
    if(use_force) force.Resize(grid,1);
    else force.Clean_Memory();
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Set_Variable_Surface_Tension(const bool use_variable_surface_tension_input)
{
    use_variable_surface_tension=use_variable_surface_tension_input;
    if(use_variable_surface_tension){
        nonzero_surface_tension=true;
        variable_surface_tension.Resize(grid.Cell_Indices(1));
        surface_tension=0;}
    else variable_surface_tension.Clean_Memory();
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Set_Variable_Viscosity(const bool use_variable_viscosity_input)
{
    use_variable_viscosity=use_variable_viscosity_input;
    if(use_variable_viscosity){
        nonzero_viscosity=true;
        variable_viscosity.Resize(grid.Cell_Indices(1));
        viscosity=0;}
    else variable_viscosity.Clean_Memory();
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Use_Variable_Vorticity_Confinement(T_GRID& grid,const bool use_variable_vorticity_confinement_input)
{
    use_variable_vorticity_confinement=use_variable_vorticity_confinement_input;
    if(use_variable_vorticity_confinement) variable_vorticity_confinement.Resize(grid.Cell_Indices(1));
    else variable_vorticity_confinement.Clean_Memory();
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Use_Variable_Vorticity_Confinement(const bool use_variable_vorticity_confinement_input)
{
    Use_Variable_Vorticity_Confinement(grid,use_variable_vorticity_confinement_input);
}
template<class T_GRID> void INCOMPRESSIBLE_UNIFORM<T_GRID>::
Use_Strain()
{
    delete strain;
    strain=new FLUID_STRAIN_UNIFORM<T_GRID>(grid);
}
//#####################################################################
template class INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<float,1> > >;
template class INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<float,2> > >;
template class INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<double,1> > >;
template class INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<double,2> > >;
template class INCOMPRESSIBLE_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
