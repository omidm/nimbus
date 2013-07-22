//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
PROJECTION_DYNAMICS_UNIFORM(const T_GRID& mac_grid,const bool flame_input,const bool multiphase,const bool use_variable_beta,const bool use_poisson,THREAD_QUEUE* thread_queue)
    :PROJECTION_COLLIDABLE_UNIFORM<T_GRID>(mac_grid,multiphase,flame_input || multiphase || use_variable_beta || use_poisson,use_variable_beta,thread_queue),PROJECTION_DYNAMICS<T>(flame_input),
    use_flame_speed_multiplier(0),dsd(0),use_divergence_multiplier_save_for_sph(false),use_non_zero_divergence_save_for_sph(false),
    p_save_for_sph(0),divergence_save_for_sph(0),divergence_multiplier_save_for_sph(0),face_velocities_save_for_sph(0),elliptic_solver_save_for_sph(0),laplace_save_for_sph(0),poisson_save_for_sph(0),
    collidable_solver_save_for_sph(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
PROJECTION_DYNAMICS_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input)
    :PROJECTION_COLLIDABLE_UNIFORM<T_GRID>(mac_grid,levelset_input),PROJECTION_DYNAMICS<T>(false),use_flame_speed_multiplier(0),dsd(0),collidable_solver_save_for_sph(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
~PROJECTION_DYNAMICS_UNIFORM()
{
    delete dsd;
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Initialize_Grid(const T_GRID& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
    if(dsd)dsd->Initialize_Grid();
}
//#####################################################################
// Function Initialize_Dsd
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Initialize_Dsd(const LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple,const ARRAY<bool>& fuel_region)
{
    int region=0;if(!fuel_region.Find(true,region)) PHYSBAM_FATAL_ERROR();//TODO: multiple fuel regions
    delete dsd;dsd=new DETONATION_SHOCK_DYNAMICS<T_GRID>(p_grid,*levelset_multiple.levelsets(region));
    if(elliptic_solver->mpi_grid)
        dsd->Set_Custom_Boundary(new BOUNDARY_MPI<T_GRID>(elliptic_solver->mpi_grid,*(dsd->boundary)),
            new BOUNDARY_MPI<T_GRID,TV>(elliptic_solver->mpi_grid,*(dsd->boundary_vector)));
}
//#####################################################################
// Function Initialize_Dsd
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Initialize_Dsd(const T_LEVELSET& levelset,const ARRAY<bool>& fuel_region)
{
    int region=0;if(!fuel_region.Find(true,region)) PHYSBAM_FATAL_ERROR();//TODO: multiple fuel regions
    delete dsd;dsd=new DETONATION_SHOCK_DYNAMICS<T_GRID>(p_grid,levelset);
    if(elliptic_solver->mpi_grid)
        dsd->Set_Custom_Boundary(new BOUNDARY_MPI<T_GRID>(elliptic_solver->mpi_grid,*dsd->boundary),
            new BOUNDARY_MPI<T_GRID,TV>(elliptic_solver->mpi_grid,*dsd->boundary_vector));
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    // find f - divergence of the velocity
    if(flame) Compute_Divergence(T_FACE_LOOKUP_FIRE_MULTIPHASE(face_velocities,*this,poisson_collidable->levelset_multiple),elliptic_solver);
    else Compute_Divergence(T_FACE_LOOKUP(face_velocities),elliptic_solver);

    // find the pressure
    elliptic_solver->Find_Solution_Regions(); // flood fill
    elliptic_solver->Compute_beta_And_Add_Jumps_To_b(dt,time); // only does something for poisson solver
    if(elliptic_solver->solve_neumann_regions) Enforce_Velocity_Compatibility(face_velocities); // make all the right hand sides compatible
    elliptic_solver->Solve(time,true); // solve all regions

    Apply_Pressure(face_velocities,dt,time);
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_SCALAR& potential_energy,const T dt)
{
    for(typename GRID<TV>::FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        if(elliptic_solver->psi_N(full_index)) continue;
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        potential_energy(full_index)-=density*(face_velocities(full_index))*(p(first_cell)-p(second_cell))/p_grid.dX(full_index.axis);}
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,T_ARRAYS_SCALAR& potential_energy,const T dt)
{
    ARRAY<TV,TV_INT> energy_gained(p_grid.Domain_Indices());
    for(typename GRID<TV>::CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()) energy_gained(iterator.Cell_Index())=TV();
    for(typename GRID<TV>::CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();for(int i=1;i<=TV::dimension;i++){
        if(elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i))) && elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i)))) continue;
        if(elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i))))
            energy_gained(cell)(i)+=face_velocities(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i)))*p(cell)*p_grid.dX(i);
        else if(elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i)))) 
            energy_gained(cell)(i)-=face_velocities(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i)))*p(cell)*p_grid.dX(i);
        else energy_gained(cell)(i)+=(face_velocities(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i)))-face_velocities(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i))))*p(cell)*p_grid.dX(i);}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()) energy_gained(iterator.Cell_Index())/=(p_grid.dX.Product());
    for(typename GRID<TV>::CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        potential_energy(cell)+=energy_gained(cell).Sum();}
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,T_ARRAYS_SCALAR& potential_energy,const T dt)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,T_FACE_ARRAYS_SCALAR& potential_energy,const T dt)
{
    for(typename GRID<TV>::FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        potential_energy(full_index)-=(T)0.5*density*(face_velocities_new(full_index)*face_velocities_new(full_index)-face_velocities_old(full_index)*face_velocities_old(full_index));}
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> template<class T_POTENTIAL> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_old,T_FACE_ARRAYS_SCALAR& face_velocities_new,T_POTENTIAL& potential_energy,T& allowed_energy_gained,const T dt)
{
    return;
    PHYSBAM_ASSERT(allowed_energy_gained==0);
    TV one_over_dx=p_grid.one_over_dX;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities(p_grid);
    for(typename GRID<TV>::FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        face_velocities(full_index)=(face_velocities_old(full_index)+face_velocities_new(full_index))/(T)2.;}
    Update_Potential_Energy(face_velocities,potential_energy,dt);
    
    for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()) for(int i=1;i<=TV::dimension;i++){
        if(elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i))) && elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i)))) continue;
        if(!elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i))) && !elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i)))) continue;
        if(elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i)))) 
            allowed_energy_gained+=density*p(iterator.Cell_Index())*face_velocities(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i)))*one_over_dx[i];
        else allowed_energy_gained-=density*p(iterator.Cell_Index())*face_velocities(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i)))*one_over_dx[i];}
    T allowed_energy_gained_tmp=allowed_energy_gained;
    for(FACE_ITERATOR iterator(p_grid,0,GRID<TV>::BOUNDARY_REGION);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        if(!elliptic_solver->psi_N(full_index)){
            if(p_grid.Domain_Indices().Lazy_Inside(iterator.First_Cell_Index())) allowed_energy_gained-=density*p(iterator.Second_Cell_Index())*face_velocities(full_index)*one_over_dx[iterator.Axis()];
            else allowed_energy_gained+=density*p(iterator.First_Cell_Index())*face_velocities(full_index)*one_over_dx[iterator.Axis()];}}
    PHYSBAM_ASSERT(abs(allowed_energy_gained-allowed_energy_gained_tmp)<1e-5);
    for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        T divergence_u_hat=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++)divergence_u_hat+=(face_velocities.Component(axis)(iterator.Second_Face_Index(axis))-face_velocities.Component(axis)(iterator.First_Face_Index(axis)))*one_over_dx[axis];
        T divergence_cell=divergence.Valid_Index(cell)?divergence(cell):0;
        allowed_energy_gained+=density*(divergence_u_hat-divergence_cell)*p(cell);}
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,2> PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Compare_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_old,T_FACE_ARRAYS_SCALAR& face_velocities_new,const T dt)
{
    //Assumes density of 1 for now
    T energy_gained_cell_total=0,energy_gained_face_total=0;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities(p_grid);
    for(typename GRID<TV>::FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        face_velocities(full_index)=(face_velocities_old(full_index)+face_velocities_new(full_index))/(T)2.;}
    
    ARRAY<T,TV_INT> energy_gained(p_grid.Domain_Indices());
    for(typename GRID<TV>::CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()) energy_gained(iterator.Cell_Index())=0;
    Update_Potential_Energy(face_velocities,energy_gained,dt);
    for(typename GRID<TV>::CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        energy_gained_cell_total+=energy_gained(cell);}

    ARRAY<T,FACE_INDEX<TV::dimension> > energy_gained_face(p_grid);
    for(typename GRID<TV>::FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()) energy_gained_face(iterator.Full_Index())=0;
    Update_Potential_Energy(face_velocities,energy_gained_face,dt);
    for(typename GRID<TV>::FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        energy_gained_face_total+=energy_gained_face(full_index);}
 
    return VECTOR<T,2>(energy_gained_cell_total,energy_gained_face_total);
}
//#####################################################################
// Function Compute_Divergence
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Compute_Divergence_For_Energy_Correction(const T_FACE_ARRAYS_SCALAR& face_velocities)
{
    divergence.Resize(p_grid.Domain_Indices());
    TV one_over_dx=p_grid.one_over_dX;
    for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
        T local_divergence=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++)local_divergence+=(face_velocities.Component(axis)(iterator.Second_Face_Index(axis))-face_velocities.Component(axis)(iterator.First_Face_Index(axis)))*one_over_dx[axis];
        divergence(iterator.Cell_Index())=local_divergence;}
}
//#####################################################################
// Function Compute_Divergence
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Compute_Divergence(const T_FACE_LOOKUP_FIRE_MULTIPHASE& face_lookup,LAPLACE_UNIFORM<T_GRID>* solver)
{
    TV one_over_dx=p_grid.one_over_dX;
    for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
        const typename T_FACE_LOOKUP_FIRE_MULTIPHASE::LOOKUP& lookup=face_lookup.Starting_Point_Cell(iterator.Cell_Index());T divergence=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++)divergence+=(lookup(axis,iterator.Second_Face_Index(axis))-lookup(axis,iterator.First_Face_Index(axis)))*one_over_dx[axis];
        solver->f(iterator.Cell_Index())=divergence;}
    
    if(use_non_zero_divergence) for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next())
        solver->f(iterator.Cell_Index())-=divergence(iterator.Cell_Index());
    if(use_divergence_multiplier) for(CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next())
        solver->f(iterator.Cell_Index())*=divergence_multiplier(iterator.Cell_Index());
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt)
{
    BASE::Apply_Pressure(face_velocities,dt,time,scale_by_dt);

    // fix the jump in pressure - interior only
    if(poisson && poisson->u_jumps){
        T_ARRAYS_BOOL& psi_D=elliptic_solver->psi_D;
        T_FACE_ARRAYS_BOOL& psi_N=elliptic_solver->psi_N;
        TV dx=p_grid.dX,one_over_dx=Inverse(dx);
        int ghost_cells=1;
        if(poisson->multiphase){
            ARRAY<T_ARRAYS_SCALAR> phis_ghost;phis_ghost.Resize(poisson_collidable->levelset_multiple->levelsets.m);
            for(int i=1;i<=poisson_collidable->levelset_multiple->levelsets.m;i++){phis_ghost(i).Resize(p_grid.Domain_Indices(ghost_cells),false);
                poisson_collidable->levelset_multiple->levelsets(i)->boundary->Fill_Ghost_Cells(p_grid,poisson_collidable->levelset_multiple->levelsets(i)->phi,phis_ghost(i),dt,time,ghost_cells);}
            LEVELSET_MULTIPLE_UNIFORM<T_GRID> levelset_multiple(p_grid,phis_ghost);
            for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(levelset_multiple.Interface(second_cell,first_cell) && !psi_N.Component(axis)(face_index) && !(psi_D(second_cell)&&psi_D(first_cell))){
                    int region_1,region_2;T phi_1,phi_2;levelset_multiple.Minimum_Regions(second_cell,first_cell,region_1,region_2,phi_1,phi_2);
                    face_velocities.Component(axis)(face_index)+=poisson->beta_face.Component(axis)(face_index)*one_over_dx[axis]*
                        LEVELSET_MULTIPLE_UNIFORM<T_GRID>::Sign(region_1,region_2)*poisson_collidable->u_jump_face.Component(axis)(face_index);}}}
        else{
            T_ARRAYS_SCALAR phi_ghost(p_grid.Domain_Indices(ghost_cells));poisson_collidable->levelset->boundary->Fill_Ghost_Cells(p_grid,poisson_collidable->levelset->phi,phi_ghost,dt,time,ghost_cells);
            for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(second_cell),phi_ghost(first_cell)) && !psi_N.Component(axis)(face_index) && !(psi_D(second_cell)&&psi_D(first_cell))){
                    face_velocities.Component(axis)(face_index)+=poisson->beta_face.Component(axis)(face_index)*one_over_dx[axis]*LEVELSET_UTILITIES<T>::Sign(phi_ghost(second_cell))*
                        LEVELSET_UTILITIES<T>::Average(phi_ghost(second_cell),poisson_collidable->u_jump(second_cell),phi_ghost(first_cell),poisson_collidable->u_jump(first_cell));}}}}
}
//#####################################################################
// Function Set_Up_For_SPH
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Set_Up_For_SPH(T_FACE_ARRAYS_SCALAR& face_velocities,const bool use_variable_density_solve,const bool use_one_way_coupling)
{
    if(use_variable_density_solve){
        POISSON_COLLIDABLE_UNIFORM<T_GRID>* poisson_for_sph=new POISSON_COLLIDABLE_UNIFORM<T_GRID>(p_grid,p,true,false,true);
        laplace_save_for_sph=laplace_collidable;elliptic_solver_save_for_sph=elliptic_solver;poisson_save_for_sph=poisson_collidable;
        collidable_solver_save_for_sph=collidable_solver;
        elliptic_solver=poisson=poisson_collidable=poisson_for_sph;laplace=laplace_collidable=0;
        collidable_solver=poisson_collidable;
        poisson->Solve_Neumann_Regions(elliptic_solver_save_for_sph->solve_neumann_regions);
        poisson->Set_Relative_Tolerance(elliptic_solver_save_for_sph->relative_tolerance);
        poisson->pcg.Set_Maximum_Iterations(elliptic_solver_save_for_sph->pcg.maximum_iterations);
        poisson->pcg.Show_Results();
        poisson->psi_N=elliptic_solver_save_for_sph->psi_N;poisson_for_sph->psi_D=elliptic_solver_save_for_sph->psi_D;
        poisson->mpi_grid=elliptic_solver_save_for_sph->mpi_grid;
        poisson_collidable->levelset=collidable_solver_save_for_sph->levelset;
        if(use_one_way_coupling){
            use_non_zero_divergence_save_for_sph=use_non_zero_divergence;
            use_divergence_multiplier_save_for_sph=use_divergence_multiplier;
            for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
                TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
                if(!elliptic_solver_save_for_sph->psi_D(cell_1) || !elliptic_solver_save_for_sph->psi_D(cell_2)) elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;}}}
    else if(use_one_way_coupling){
        face_velocities_save_for_sph=new T_FACE_ARRAYS_SCALAR(face_velocities);
        p_save_for_sph=new T_ARRAYS_SCALAR(p);
        divergence_save_for_sph=new T_ARRAYS_SCALAR(divergence);
        divergence_multiplier_save_for_sph=new T_ARRAYS_SCALAR(divergence_multiplier);
        use_divergence_multiplier_save_for_sph=use_divergence_multiplier;
        use_non_zero_divergence_save_for_sph=use_non_zero_divergence;
        elliptic_solver->psi_D_save_for_sph=new T_ARRAYS_BOOL(elliptic_solver->psi_D);
        elliptic_solver->psi_N_save_for_sph=new T_FACE_ARRAYS_BOOL(elliptic_solver->psi_N);
        for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            if(!(*elliptic_solver->psi_D_save_for_sph)(cell_1) || !(*elliptic_solver->psi_D_save_for_sph)(cell_2)) elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;}}

    Use_Divergence_Multiplier(true);Use_Non_Zero_Divergence(true);
}
//#####################################################################
// Function Restore_After_SPH
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Restore_After_SPH(T_FACE_ARRAYS_SCALAR& face_velocities,const bool use_variable_density_solve,const bool use_one_way_coupling)
{
    if(use_variable_density_solve){
        delete poisson;poisson=poisson_collidable=poisson_save_for_sph;
        laplace=laplace_collidable=laplace_save_for_sph;elliptic_solver=elliptic_solver_save_for_sph;collidable_solver=collidable_solver_save_for_sph;
        if(use_one_way_coupling){Use_Divergence_Multiplier(use_divergence_multiplier_save_for_sph);Use_Non_Zero_Divergence(use_non_zero_divergence_save_for_sph);}}
    else if(use_one_way_coupling){
        face_velocities=*face_velocities_save_for_sph;
        delete face_velocities_save_for_sph;face_velocities_save_for_sph=0;
        p=*p_save_for_sph;
        delete p_save_for_sph;p_save_for_sph=0;
        divergence=*divergence_save_for_sph;
        delete divergence_save_for_sph;divergence_save_for_sph=0;
        divergence_multiplier=*divergence_multiplier_save_for_sph;
        delete divergence_multiplier_save_for_sph;divergence_multiplier_save_for_sph=0;
        Use_Divergence_Multiplier(use_divergence_multiplier_save_for_sph);Use_Non_Zero_Divergence(use_non_zero_divergence_save_for_sph);
        elliptic_solver->psi_D=*elliptic_solver->psi_D_save_for_sph;
        delete elliptic_solver->psi_D_save_for_sph;elliptic_solver->psi_D_save_for_sph=0;
        elliptic_solver->psi_N=*elliptic_solver->psi_N_save_for_sph;
        delete elliptic_solver->psi_N_save_for_sph;elliptic_solver->psi_N_save_for_sph=0;}
}
//#####################################################################
// Function Update_Phi_And_Move_Velocity_Discontinuity
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Update_Phi_And_Move_Velocity_Discontinuity(T_FACE_ARRAYS_SCALAR& face_velocities,LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple,const T time,const bool update_phi_only)
{
    assert(flame);
    int ghost_cells=3;
    levelset_multiple.Fill_Ghost_Cells(levelset_multiple.phis,time,ghost_cells);
    if(!update_phi_only){
        LEVELSET_MULTIPLE<T_GRID>& levelset_multiple_old=*poisson_collidable->levelset_multiple;
        for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            int region_old=levelset_multiple_old.Inside_Region_Face(axis,face),region_new=levelset_multiple.Inside_Region_Face(axis,face);
            if(region_old!=region_new) face_velocities.Component(axis)(face)-=Face_Jump_Multiphase(axis,face,region_new,region_old);}}
    poisson_collidable->Update_Internal_Level_Set(levelset_multiple);poisson_collidable->levelset_multiple->Compute_Normals();
}
//#####################################################################
// Function Update_Phi_And_Move_Velocity_Discontinuity
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Flame_Speed_Face_Multiphase(const int axis,const TV_INT& face_index,const int fuel_region,const int product_region) const
{
    T multiplier=1;if(use_flame_speed_multiplier) multiplier=flame_speed_multiplier.Component(axis)(face_index);
    if(dsd) return multiplier*dsd->Normal_Flame_Speed(axis,face_index);
    TV_INT offset=TV_INT::Axis_Vector(axis);const TRIPLE<T,T,T>& constants=flame_speed_constants(fuel_region,product_region);
    const T normal_flame_speed=constants.x;const T curvature_flame_speed=constants.y;
    if(!curvature_flame_speed) return multiplier*normal_flame_speed;
    const T_LEVELSET* levelset=poisson_collidable->levelset_multiple->levelsets(fuel_region);
    T face_curvature=(T).5*((*levelset->curvature)(face_index)+(*levelset->curvature)(face_index-offset));
    return multiplier*(normal_flame_speed+curvature_flame_speed*face_curvature);
}
//#####################################################################
// Function Use_Flame_Speed_Multiplier
//#####################################################################
template<class T_GRID> void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Use_Flame_Speed_Multiplier(const bool use_flame_speed_multiplier_input)
{
    use_flame_speed_multiplier=use_flame_speed_multiplier_input;
    if(use_flame_speed_multiplier) flame_speed_multiplier.Resize(p_grid,3);
    else flame_speed_multiplier.Clean_Memory();
}
//#####################################################################
// Function Face_Jump_Multiphase
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Face_Jump_Multiphase(const int axis,const TV_INT& face_index,const int current_region,const int face_region) const
{
    const TRIPLE<T,T,T>& constants=flame_speed_constants(current_region,face_region);
    if(constants.z==0) return 0; // getting the diagonal if the two regions are the same. i.e. .z will be zero and the jump will return as nothing
    TV_INT offset=TV_INT::Axis_Vector(axis);
    const T_LEVELSET* levelset=poisson_collidable->levelset_multiple->levelsets(current_region);
    T face_normal=(levelset->phi(face_index)-levelset->phi(face_index-offset))*p_grid.one_over_dX[axis];
    return constants.z*Flame_Speed_Face_Multiphase(axis,face_index,current_region,face_region)*face_normal;
}
//#####################################################################
#define INSTANTIATION_HELPER(T_GRID) \
    template class PROJECTION_DYNAMICS_UNIFORM<T_GRID>; \
    template void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::Update_Potential_Energy(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& face_velocities_old,GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& face_velocities_new,GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& potential_energy,T&,const T dt); \
    template void PROJECTION_DYNAMICS_UNIFORM<T_GRID>::Update_Potential_Energy(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& face_velocities_old,GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS& face_velocities_new,GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR& potential_energy,T&,const T dt);

#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >));
#endif
