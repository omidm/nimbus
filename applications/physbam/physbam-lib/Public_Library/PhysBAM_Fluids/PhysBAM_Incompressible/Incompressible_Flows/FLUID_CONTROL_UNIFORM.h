//#####################################################################
// Copyright 2006, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_CONTROL_UNIFORM  
//#####################################################################
#ifndef __FLUID_CONTROL_UNIFORM__
#define __FLUID_CONTROL_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Computations/SMOOTH_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Solids_And_Fluids/FLUID_CONTROL_CALLBACKS.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class FLUID_CONTROL_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
public:
    const FLUID_CONTROL_CALLBACKS<T_GRID>* callbacks;
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    T_LEVELSET& levelset;
    T_FACE_ARRAYS_SCALAR& face_velocities;
    T_FACE_ARRAYS_BOOL& psi_N;
    PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection;
    T alpha; // shape control gain
    T beta; // derivative control gain
    T C; // strength of potential field
    int shape_smoothing_steps;
    int velocity_smoothing_steps;
    int potential_smoothing_steps;
    int pcg_iterations;
    int pcg_iterations_mpi;

    bool use_control_falloff;
    bool use_custom_shape_distance_falloff;
    bool use_custom_velocity_distance_falloff;
    bool use_custom_potential_distance_falloff;
    
    bool left_wall,right_wall,bottom_wall,top_wall,front_wall,back_wall;

    FLUID_CONTROL_UNIFORM(T_LEVELSET& levelset_input,T_FACE_ARRAYS_SCALAR& face_velocities_input,T_FACE_ARRAYS_BOOL& psi_N_input,PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input,
                          bool left_wall_input,bool right_wall_input,bool bottom_wall_input,bool top_wall_input,bool front_wall_input,bool back_wall_input)
        :mpi_grid(0),levelset(levelset_input),face_velocities(face_velocities_input),psi_N(psi_N_input),projection(projection_input),
         alpha(625),beta(25),shape_smoothing_steps(0),velocity_smoothing_steps(0),potential_smoothing_steps(0),pcg_iterations(20),pcg_iterations_mpi(80),
         use_control_falloff(false),use_custom_shape_distance_falloff(false),use_custom_velocity_distance_falloff(false),
         left_wall(left_wall_input),right_wall(right_wall_input),bottom_wall(bottom_wall_input),top_wall(top_wall_input),front_wall(front_wall_input),back_wall(back_wall_input)
    {}

    void Set_Fluid_Control_Callbacks(const FLUID_CONTROL_CALLBACKS<T_GRID>& callbacks_input)
    {callbacks=&callbacks_input;}

    void Smooth_Face_Array(BOUNDARY_UNIFORM<T_GRID,TV>* boundary,const T_FACE_ARRAYS_SCALAR& array,T_FACE_ARRAYS_SCALAR& smoothed_array,const int smoothing_steps)
    {int ghost_cells=3;
    T_ARRAYS_VECTOR cell_array(levelset.grid.Domain_Indices(ghost_cells));
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){
        for(int axis=1;axis<=levelset.grid.dimension;axis++) cell_array(iterator.Cell_Index())[axis]=(T).5*(array(axis,iterator.First_Face_Index(axis))+array(axis,iterator.Second_Face_Index(axis)));}
    for(int i=1;i<=smoothing_steps;i+=1){
        boundary->Fill_Ghost_Cells(levelset.grid,cell_array,cell_array,0,0,ghost_cells); // TODO: use real time/dt
        SMOOTH::Smooth<T_GRID>(cell_array,1,0);}
    boundary->Fill_Ghost_Cells(levelset.grid,cell_array,cell_array,0,0,ghost_cells); // TODO: use real dt/time
    for(FACE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        smoothed_array(axis,iterator.Face_Index())=(T).5*(cell_array(iterator.First_Cell_Index())[axis]+cell_array(iterator.Second_Cell_Index())[axis]);}}

    void Get_Fluid_Control_Force(T_FACE_ARRAYS_SCALAR& force,const T time)
    {int ghost_cells=3;
    T_FACE_ARRAYS_SCALAR f_velocity(levelset.grid);
    T_FACE_ARRAYS_SCALAR* processed_face_velocities=&face_velocities;
    BOUNDARY_UNIFORM<T_GRID,TV>* boundary=new BOUNDARY_UNIFORM<T_GRID,TV>();
    if(mpi_grid)boundary=new BOUNDARY_MPI<T_GRID,TV>(mpi_grid,*boundary);
    BOUNDARY_UNIFORM<T_GRID,T>* phi_boundary=new BOUNDARY_UNIFORM<T_GRID,T>();
    if(mpi_grid)phi_boundary=new BOUNDARY_MPI<T_GRID>(mpi_grid,*phi_boundary);
    phi_boundary->Fill_Ghost_Cells(levelset.grid,levelset.phi,levelset.phi,0,time,ghost_cells); // TODO: use real dt
    if(velocity_smoothing_steps){
        processed_face_velocities=new T_FACE_ARRAYS_SCALAR(levelset.grid,ghost_cells);
        Smooth_Face_Array(boundary,face_velocities,*processed_face_velocities,velocity_smoothing_steps);}
    projection.elliptic_solver->mpi_grid=mpi_grid;
    T_FACE_ARRAYS_SCALAR& f_shape=face_velocities; //TODO: Was this ever different?
    f_shape.Fill(0);
    projection.elliptic_solver->psi_N.Fill(false);
    projection.elliptic_solver->psi_D.Fill(false);
    T_ARRAYS_SCALAR potential(levelset.grid.Domain_Indices(ghost_cells));
    T_FACE_ARRAYS_SCALAR f_potential(levelset.grid);
    TV one_over_dx=Inverse(levelset.grid.DX());
    // potential field
    for(CELL_ITERATOR it(levelset.grid,ghost_cells);it.Valid();it.Next()){TV pos = it.Location();
        T target_phi=callbacks->Phi(pos,time);
        potential(it.Cell_Index())=C*target_phi;}
    phi_boundary->Fill_Ghost_Cells(levelset.grid,potential,potential,0,time,ghost_cells); // TODO: use real dt
    for(FACE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        T phi1=levelset.phi(iterator.First_Cell_Index()),phi2=levelset.phi(iterator.Second_Cell_Index());
        // derivative control
        f_velocity(axis,face_index)+=-beta*((*processed_face_velocities)(axis,face_index)-callbacks->Face_Velocity(iterator.Location(),iterator.Axis(),time));
        if(use_custom_velocity_distance_falloff) f_velocity(axis,face_index)*=callbacks->Velocity_Distance_Falloff(callbacks->Phi(iterator.Location(),time),time);
        // shape control
        if(phi1>0 || phi2>0){
            projection.laplace->psi_N(axis,face_index)=true;
            T target_phi=callbacks->Phi(iterator.Location(),time)-0.5*(phi1+phi2);
            T falloff = (use_custom_shape_distance_falloff)?callbacks->Shape_Distance_Falloff(target_phi,time):fabs(target_phi);
            if(target_phi>0) f_shape(axis,face_index)+=-alpha*falloff*callbacks->Normal(iterator.Location(),time)[axis];
            else f_shape(axis,face_index)+=alpha*falloff*(phi2-phi1)*one_over_dx[axis];}
        // potential control
        f_potential(axis,face_index)+=-(potential(iterator.Second_Cell_Index())-potential(iterator.First_Cell_Index()))*one_over_dx[axis];
        if(use_custom_potential_distance_falloff) f_potential(axis,face_index)*=callbacks->Potential_Distance_Falloff(callbacks->Phi(iterator.Location(),time),time);}
    if(velocity_smoothing_steps){
        // deallocate previously allocated face velocities.
        delete processed_face_velocities;
    }
    if(use_control_falloff){
        T_FACE_ARRAYS_SCALAR control_falloff_value(levelset.grid);
        for(FACE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next())
            control_falloff_value(iterator.Axis(),iterator.Face_Index()) = callbacks->Control_Falloff(iterator.Location(),time);
        if(mpi_grid)mpi_grid->Average_Common_Face_Data(control_falloff_value);
        for(FACE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            T control_falloff = control_falloff_value(axis,face_index);
            f_shape(axis,face_index)*=control_falloff;f_velocity(axis,face_index)*=control_falloff;f_potential(axis,face_index)*=control_falloff;
            if(control_falloff<1e-6)projection.laplace->psi_N(axis,face_index)=true;}}
    for(FACE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        if(psi_N(axis,face_index)){
            f_velocity(axis,face_index)=0;f_shape(axis,face_index)=0;f_potential(axis,face_index)=0;
            projection.laplace->psi_N(axis,face_index)=true;}}

    if (alpha!=0){
        for(CELL_ITERATOR iterator(levelset.grid,1);iterator.Valid();iterator.Next()) 
            if(levelset.phi(iterator.Cell_Index())>0){projection.laplace->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=0;}

        VECTOR<VECTOR<bool,2>,T_GRID::dimension> domain_walls;
        for(int i=1;i<=3;i++) for(int j=1;j<=2;j++) domain_walls[i][j]=true;
        if(mpi_grid)mpi_grid->Initialize(domain_walls);
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            for(int side=1;side<=2;side++){
                int side_number=(axis-1)*2+side;
                if(domain_walls[axis][side]) for(CELL_ITERATOR iterator(levelset.grid,1,T_GRID::GHOST_REGION,side_number);iterator.Valid();iterator.Next()){
                    projection.laplace->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=0;}}}

        if(mpi_grid){mpi_grid->Exchange_Boundary_Cell_Data(projection.elliptic_solver->psi_D,1,false);mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
        projection.laplace->Solve_Neumann_Regions();

        if(mpi_grid) projection.elliptic_solver->pcg.Set_Maximum_Iterations(pcg_iterations_mpi);
        else projection.elliptic_solver->pcg.Set_Maximum_Iterations(pcg_iterations);

        projection.Make_Divergence_Free(face_velocities,0,time); // TODO: use real dt

        if(shape_smoothing_steps) Smooth_Face_Array(boundary,f_shape,f_shape,shape_smoothing_steps);}
    if(potential_smoothing_steps)
        Smooth_Face_Array(boundary,f_potential,f_potential,potential_smoothing_steps);
    for(FACE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        force(axis,face_index)+=f_velocity(axis,face_index)+f_shape(axis,face_index)+f_potential(axis,face_index);}
    if (mpi_grid) {
        delete &((BOUNDARY_MPI<T_GRID,TV>*)boundary)->boundary;
        delete &((BOUNDARY_MPI<T_GRID,T>*)phi_boundary)->boundary;
    }
    delete boundary;
    delete phi_boundary;
    }
//#####################################################################
};
}
#endif
