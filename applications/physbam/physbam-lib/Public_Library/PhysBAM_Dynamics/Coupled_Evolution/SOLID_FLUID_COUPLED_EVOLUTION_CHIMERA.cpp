//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Sergey Levine, Nick Rasmussen, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/INTERRUPTS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/VOF_ADVECTION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_CHIMERA.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_CHIMERA.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/CHIMERA_SYSTEM.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Dynamics/Meshing/POLYGONAL_TRIANGULATION.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_CHIMERA.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/POLYGON_HYPERPLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/AVERAGING_CHIMERA.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/CELL_LOOKUP_CHIMERA.h>
#include <PhysBAM_Tools/Math_Tools/Is_NaN.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/LINEAR_INTERPOLATION_CHIMERA.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/VORONOI_INTERPOLATION_CHIMERA.h>
#include <iomanip>
using namespace PhysBAM;

//#####################################################################
// Macros for chimera grid and data access
//#####################################################################
#define FOR_EACH_LOCAL_GRID(I) for(int I=1;I<=example.rigid_grid_collection.particles.array_collection->Size();I++)
#define LOCAL_RIGID_GRID_ACCESS(I) example.rigid_grid_collection.Rigid_Grid(I)
#define LOCAL_GRID_ACCESS(I) example.rigid_grid_collection.Rigid_Grid(I).grid
#define FACE_VELOCITIES(I) example.incompressible_fluid_containers(I)->face_velocities
#define PHI(I) example.incompressible_fluid_containers(I)->particle_levelset_evolution.phi
#define DENSITY(I) example.incompressible_fluid_containers(I)->density_container.density
#define TEMPERATURE(I) example.incompressible_fluid_containers(I)->temperature_container.temperature
#define PRESSURE(I) example.incompressible_fluid_containers(I)->pressure
#define PSI_N(I) example.incompressible_fluid_containers(I)->psi_N
#define PSI_D(I) example.incompressible_fluid_containers(I)->psi_D
#define FOR_EACH_GLOBAL_GRID(I) for(int I=1;I<=n_global_grids;I++)
#define GLOBAL_RIGID_GRID_ACCESS(I) example.chimera_grid->global_grid.Rigid_Grid(I)
#define GLOBAL_GRID_ACCESS(I) example.chimera_grid->global_grid.Rigid_Grid(I).grid

//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Advance_One_Time_Step(const T time,const T dt)
{
    Write_Substep("before update grid",0,1);
    int substep=0;
    example.Update_Grid(time,dt);
    example.chimera_grid->Adjustment_After_Updating_Grids();
    n_global_grids=example.chimera_grid->number_of_global_grids;
    n_local_grids=example.rigid_grid_collection.particles.array_collection->Size();
    Write_Substep("after update grid",0,1);

    Setup_Solids(time);
    Setup_Fluids(time);
    
    example.Update_Fluid_Parameters(dt,time); 
    Write_Substep("after setup",0,1);
    
    //Write_Substep("integrate fluid forces for solid coupling",substep,1);
    
    Solid_Position_Update(time,dt);/*S1*/
    
    Write_Substep("object compatibility",substep,1);
    
    //incompressible->projection.Restore_After_Projection();

    example.Adjust_Density_And_Temperature_With_Sources(time);
    
    Write_Substep("after adjust_density_and_temperature",0,1);
    
    Write_Substep("advect fluid",substep,1);
    if(!example.solve_poisson_equation) Advect_Fluid(time,dt);/*F2*/
    Write_Substep("after advect_fluid",0,1);
    
    Solid_Velocity_Update(time,dt);/*S2*/
    Write_Substep("after solid velocity update",substep,1);

    Write_Substep("after get source velocity before projection",substep,1);
    if(!example.run_advection_test){
        Advance_Fluid_One_Time_Step_Implicit_Part(time,dt);/*FS3*/
        Write_Substep("after advance fluid implicit part",substep,1);

        Set_Zero_Air_Velocities();
        Extrapolate_Velocity_Across_Interface_On_All_Local_Grids();
        Coupling_Vector_Fields();

        //doesn't work for mpi
        /*LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
            if(laplace_grid.Local_Grid(grid_index)){
                //int local_grid_index=example.chimera_grid->global_to_local_grid_index_map(grid_index);
                T_GRID& grid=laplace_grid.Grid(grid_index);
                const FRAME<TV>& frame=laplace_grid.Frame(grid_index);
                for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                    //LOG::cout << "iterator " << grid_index << " " << iterator.Full_Index() << std::endl;
                    if(laplace_grid.Boundary_Cell(grid_index,iterator.First_Cell_Index()) || laplace_grid.Boundary_Cell(grid_index,iterator.Second_Cell_Index())){
                        D_FACE_INDEX face_index=iterator.Full_Index();
                        TV location=frame*iterator.Location();
                        int finest_grid_index=example.chimera_grid->Find_Finest_Grid(location);
                        if(finest_grid_index && finest_grid_index!=grid_index){
                            T_GRID& finest_grid=laplace_grid.Grid(finest_grid_index);
                            const FRAME<TV>& finest_frame=laplace_grid.Frame(finest_grid_index);
                            TV normal=frame.r.Rotate(TV::Axis_Vector(face_index.axis));
                            TV interpolated_vector=finest_frame.r.Rotate(interpolation.Clamped_To_Array_Face(finest_grid,FACE_VELOCITIES(finest_grid_index),finest_frame.Inverse_Times(location)));
                            FACE_VELOCITIES(grid_index)(face_index)=TV::Dot_Product(normal,interpolated_vector);}}}}*/
        
        if(example.solve_heat_equation) Coupling_Scalar_Fields(time,dt);
    }
}
//#####################################################################
// Function Setup_Solids
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Calculate_Maximum_Allowable_dt(const T time)
{
    if(example.solve_poisson_equation) return (T)1;

    T max_dt=1,min_dt=0;
    if(example.fixed_dt) max_dt=example.fixed_dt;
    if(!example.cfl_alpha || !example.cfl_beta) PHYSBAM_FATAL_ERROR("More ghost cells should be allocated!");
    //fluid velocity constraint
    ARRAY<T_FACE_ARRAYS_SCALAR*> u_array(n_local_grids);
    FOR_EACH_LOCAL_GRID(grid_index){
        u_array(grid_index)=&FACE_VELOCITIES(grid_index);}
    example.Fill_Ghost_Cells_Face_Chimera((int)std::ceil(example.cfl_alpha),u_array,u_array);
    FOR_EACH_LOCAL_GRID(grid_index){
        T_GRID& grid=LOCAL_GRID_ACCESS(grid_index);T min_dX=grid.dX.Min();
        for(CELL_ITERATOR iterator(grid,(const int)std::ceil(example.cfl_alpha));iterator.Valid();iterator.Next()){
            TV local_V;
            for(int axis=1;axis<=T_GRID::dimension;++axis)
                local_V(axis)=maxabs(FACE_VELOCITIES(grid_index)(axis,iterator.First_Face_Index(axis)),
                    FACE_VELOCITIES(grid_index)(axis,iterator.Second_Face_Index(axis)));
            T local_V_norm=local_V.Magnitude();
            if(local_V_norm>((T)example.cfl_beta*min_dX/max_dt)) max_dt=example.cfl_beta*min_dX/local_V_norm;}}
    //sync dt_fluid in all processes (MPI)
    if(example.chimera_grid->use_mpi){example.chimera_grid->Synchronize_Dt(max_dt);}
    //grid motion constraint (using bisection iteration)
    T tolerance=GLOBAL_GRID_ACCESS(1).dX.Min()*(T)1e-3;
    ARRAY<TV> X(n_global_grids);ARRAY<TV> V(n_global_grids);ARRAY<ROTATION<TV> > rotation(n_global_grids);ARRAY<typename TV::SPIN> angular(n_global_grids);
    //first check whether max_dt satisfy the grid motion constraint
    FOR_EACH_GLOBAL_GRID(i){
        X(i)=GLOBAL_RIGID_GRID_ACCESS(i).X();
        rotation(i)=GLOBAL_RIGID_GRID_ACCESS(i).Rotation();
        V(i)=GLOBAL_RIGID_GRID_ACCESS(i).V();
        angular(i)=GLOBAL_RIGID_GRID_ACCESS(i).Angular();}
    example.Update_Grid(X,V,rotation,angular,time,max_dt);
    bool bounded=true;//This indicates all corners of all grids lie inside bounding boxes. If any corner found outside, this is set to false.
    FOR_EACH_GLOBAL_GRID(grid_index){
        T_GRID& grid=GLOBAL_GRID_ACCESS(grid_index);RIGID_GRID<T_GRID>& time_n_rigid_grid=GLOBAL_RIGID_GRID_ACCESS(grid_index);
        RANGE<TV> domain=grid.Domain();domain.Change_Size(grid.DX()*(T)(example.cfl_alpha-0.5-tolerance));
        ORIENTED_BOX<TV> bounding_box(domain,time_n_rigid_grid.Frame());
        int number_of_corners=1;for(int i=1;i<=TV::dimension;i++){number_of_corners*=2;}
        for(int i=0;i<number_of_corners;i++){
            TV corner_point;
            for(int axis=1;axis<=TV::dimension;axis++){
                int bit=(i>>(axis-1)) & 1;
                corner_point(axis)=bit?grid.domain.max_corner(axis):grid.domain.min_corner(axis);}
            TV time_n_plus_one_location=FRAME<TV>(X(grid_index),rotation(grid_index))*corner_point;
            if(!bounding_box.Lazy_Inside(time_n_plus_one_location)){bounded=false;break;}}
        if(!bounded) break;}
    if(bounded && (example.use_maccormack_advection_scalar || example.use_maccormack_advection_vector)){
        FOR_EACH_GLOBAL_GRID(grid_index){
            T_GRID& grid=GLOBAL_GRID_ACCESS(grid_index);RIGID_GRID<T_GRID>& time_n_rigid_grid=GLOBAL_RIGID_GRID_ACCESS(grid_index);
            RANGE<TV> domain=grid.Domain();domain.Change_Size(grid.DX()*(T)(example.cfl_alpha-0.5-tolerance));
            ORIENTED_BOX<TV> bounding_box(domain,FRAME<TV>(X(grid_index),rotation(grid_index)));
            int number_of_corners=1;for(int i=1;i<=TV::dimension;i++){number_of_corners*=2;}
            for(int i=0;i<number_of_corners;i++){
                TV corner_point;
                for(int axis=1;axis<=TV::dimension;axis++){
                    int bit=(i>>(axis-1)) & 1;
                    corner_point(axis)=bit?grid.domain.max_corner(axis):grid.domain.min_corner(axis);}
                TV time_n_location=time_n_rigid_grid.Frame()*corner_point;
                if(!bounding_box.Lazy_Inside(time_n_location)){bounded=false;break;}}
            if(!bounded) break;}}
    if(bounded){
      ////////////////////////////////////////////////////////////////////////////////////////////DEBUG OUTPUT
      /*FOR_EACH_GLOBAL_GRID(grid_index){
          T_GRID& grid=GLOBAL_GRID_ACCESS(grid_index);RIGID_GRID<T_GRID>& time_n_rigid_grid=GLOBAL_RIGID_GRID_ACCESS(grid_index);
          RANGE<TV> domain=grid.Domain();domain.Change_Size(grid.DX()*(T)(example.cfl_alpha-0.5));
          ORIENTED_BOX<TV> bounding_box(domain,time_n_rigid_grid.Frame());
          int number_of_corners=1;for(int i=1;i<=TV::dimension;i++){number_of_corners*=2;}
          for(int i=0;i<number_of_corners;i++){
              TV corner_point;
              for(int axis=1;axis<=TV::dimension;axis++){
                  int bit=(i>>(axis-1)) & 1;
                  corner_point(axis)=bit?grid.domain.max_corner(axis):grid.domain.min_corner(axis);}
              TV time_n_plus_one_location=FRAME<TV>(X(grid_index),rotation(grid_index))*corner_point;
              T distance=bounding_box.Signed_Distance(time_n_plus_one_location);
              LOG::cout<<"lqiu debug corner distance:"<<distance<<std::endl;
          }
      }*/
      ///////////////////////////////////////////////////////////////////////////////////////////////////////
        LOG::cout<<"dt="<<max_dt<<std::endl;return max_dt;}
  
    int iteration_count=0;
    while(abs(max_dt-min_dt)>1e-6){
        T mid_dt=(max_dt+min_dt)/(T)2.0;
        iteration_count++;
        FOR_EACH_GLOBAL_GRID(i){
            X(i)=GLOBAL_RIGID_GRID_ACCESS(i).X();
            rotation(i)=GLOBAL_RIGID_GRID_ACCESS(i).Rotation();
            V(i)=GLOBAL_RIGID_GRID_ACCESS(i).V();
            angular(i)=GLOBAL_RIGID_GRID_ACCESS(i).Angular();}
        example.Update_Grid(X,V,rotation,angular,time,mid_dt);
        bounded=true;//Set initial value to be true.
        FOR_EACH_GLOBAL_GRID(grid_index){
            T_GRID& grid=GLOBAL_GRID_ACCESS(grid_index);RIGID_GRID<T_GRID>& time_n_rigid_grid=GLOBAL_RIGID_GRID_ACCESS(grid_index);
            RANGE<TV> domain=grid.Domain();domain.Change_Size(grid.DX()*(T)(example.cfl_alpha-0.5-tolerance));
            ORIENTED_BOX<TV> bounding_box(domain,time_n_rigid_grid.Frame());
            int number_of_corners=1;for(int i=1;i<=TV::dimension;i++){number_of_corners*=2;}
            for(int i=0;i<number_of_corners;i++){
                TV corner_point;
                for(int axis=1;axis<=TV::dimension;axis++){
                    int bit=(i>>(axis-1)) & 1;
                    corner_point(axis)=bit?grid.domain.max_corner(axis):grid.domain.min_corner(axis);}
                TV time_n_plus_one_location=FRAME<TV>(X(grid_index),rotation(grid_index))*corner_point;
                if(!bounding_box.Lazy_Inside(time_n_plus_one_location)){bounded=false;break;}}
            if(!bounded) break;}
        if(bounded && (example.use_maccormack_advection_scalar || example.use_maccormack_advection_vector)){
            FOR_EACH_GLOBAL_GRID(grid_index){
                T_GRID& grid=GLOBAL_GRID_ACCESS(grid_index);RIGID_GRID<T_GRID>& time_n_rigid_grid=GLOBAL_RIGID_GRID_ACCESS(grid_index);
                RANGE<TV> domain=grid.Domain();domain.Change_Size(grid.DX()*(T)(example.cfl_alpha-0.5-tolerance));
                ORIENTED_BOX<TV> bounding_box(domain,FRAME<TV>(X(grid_index),rotation(grid_index)));
                int number_of_corners=1;for(int i=1;i<=TV::dimension;i++){number_of_corners*=2;}
                for(int i=0;i<number_of_corners;i++){
                    TV corner_point;
                    for(int axis=1;axis<=TV::dimension;axis++){
                        int bit=(i>>(axis-1)) & 1;
                        corner_point(axis)=bit?grid.domain.max_corner(axis):grid.domain.min_corner(axis);}
                    TV time_n_location=time_n_rigid_grid.Frame()*corner_point;
                    if(!bounding_box.Lazy_Inside(time_n_location)){bounded=false;break;}}
                if(!bounded) break;}}
        if(bounded) min_dt=mid_dt;
        if(!bounded) max_dt=mid_dt;} 
////////////////////////////////////////////////////////////////////////////////////////////DEBUG OUTPUT
  /*LOG::cout<<"lqiu debug number of iterations:"<<iteration_count<<std::endl;
  FOR_EACH_GLOBAL_GRID(i){
      X(i)=GLOBAL_RIGID_GRID_ACCESS(i).X();
      rotation(i)=GLOBAL_RIGID_GRID_ACCESS(i).Rotation();
      V(i)=GLOBAL_RIGID_GRID_ACCESS(i).V();
      angular(i)=GLOBAL_RIGID_GRID_ACCESS(i).Angular();}
  example.Update_Grid(X,V,rotation,angular,time,min_dt);
  FOR_EACH_GLOBAL_GRID(grid_index){
      T_GRID& grid=GLOBAL_GRID_ACCESS(grid_index);RIGID_GRID<T_GRID>& time_n_rigid_grid=GLOBAL_RIGID_GRID_ACCESS(grid_index);
          RANGE<TV> domain=grid.Domain();domain.Change_Size(grid.DX()*(T)(example.cfl_alpha-0.5));
          ORIENTED_BOX<TV> bounding_box(domain,time_n_rigid_grid.Frame());
          int number_of_corners=1;for(int i=1;i<=TV::dimension;i++){number_of_corners*=2;}
          for(int i=0;i<number_of_corners;i++){
              TV corner_point;
              for(int axis=1;axis<=TV::dimension;axis++){
                  int bit=(i>>(axis-1)) & 1;
                  corner_point(axis)=bit?grid.domain.max_corner(axis):grid.domain.min_corner(axis);}
              TV time_n_plus_one_location=FRAME<TV>(X(grid_index),rotation(grid_index))*corner_point;
              T distance=bounding_box.Signed_Distance(time_n_plus_one_location);
              if(distance>0){LOG::cout<<"lqiu debug lazy inside:"<<bounding_box.Lazy_Inside(time_n_plus_one_location)<<std::endl;}
              LOG::cout<<"lqiu debug corner distance:"<<distance<<std::endl;}}
  if(use_maccormack_advection) FOR_EACH_GLOBAL_GRID(grid_index){
      T_GRID& grid=GLOBAL_GRID_ACCESS(grid_index);RIGID_GRID<T_GRID>& time_n_rigid_grid=GLOBAL_RIGID_GRID_ACCESS(grid_index);
      RANGE<TV> domain=grid.Domain();domain.Change_Size(grid.DX()*(T)(example.cfl_alpha-0.5));
      ORIENTED_BOX<TV> bounding_box(domain,FRAME<TV>(X(grid_index),rotation(grid_index)));
      int number_of_corners=1;for(int i=1;i<=TV::dimension;i++){number_of_corners*=2;}
      for(int i=0;i<number_of_corners;i++){
          TV corner_point;
          for(int axis=1;axis<=TV::dimension;axis++){
              int bit=(i>>(axis-1)) & 1;
              corner_point(axis)=bit?grid.domain.max_corner(axis):grid.domain.min_corner(axis);}
          TV time_n_location=time_n_rigid_grid.Frame()*corner_point;
          T distance=bounding_box.Signed_Distance(time_n_location);
          if(distance>0){LOG::cout<<"lqiu debug lazy inside:"<<bounding_box.Lazy_Inside(time_n_location)<<std::endl;}
          LOG::cout<<"lqiu debug corner distance (backward):"<<distance<<std::endl;}}*/
///////////////////////////////////////////////////////////////////////////////////////////////////////
    LOG::cout<<"dt="<<min_dt<<std::endl;
    return min_dt;
}
//#####################################################################
// Function Setup_Solids
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Setup_Solids(const T time)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;
    int substep=0;
    
    if(solids_parameters.triangle_collision_parameters.perform_self_collision && solids_parameters.triangle_collision_parameters.temporary_enable_collisions){
        solids_evolution_callbacks->Self_Collisions_Begin_Callback(time,substep);
        solids_parameters.triangle_collision_parameters.repulsion_pair_update_count=0;
        example.solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Save_Self_Collision_Free_State();
        if((solids_parameters.triangle_collision_parameters.topological_hierarchy_build_count++)%solids_parameters.triangle_collision_parameters.topological_hierarchy_build_frequency==0){
            LOG::SCOPE scope("hierarchybuild","Building Hierarchy Topology");
            example.solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Build_Topological_Structure_Of_Hierarchies();}
        solids_parameters.triangle_collision_parameters.self_collision_free_time=time;}

    if(solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies) example.solid_body_collection.rigid_body_collection.Reset_Impulse_Accumulators();
    solids_evolution_callbacks->Preprocess_Solids_Substep(time,substep);
    if(solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects) // TODO - ANDY - why is this needed??? TODO: move this to the right places inside solids evolution 
        example.solid_body_collection.collision_body_list.Update_Spatial_Partition(solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
            solids_parameters.deformable_object_collision_parameters.spatial_partition_number_of_cells,solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);
    example.solid_body_collection.Update_Time_Varying_Material_Properties(time);
}
//#####################################################################
// Function Setup_Fluids
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Setup_Fluids(const T time)
{
    //const T fictitious_dt=(T)1./example.frame_rate;
    FOR_EACH_LOCAL_GRID(grid_index){
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);
        example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Set_Number_Particles_Per_Cell(example.fluids_parameters.number_particles_per_cell);}
}
//#####################################################################
// Function Advect_Fluid
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Advect_Fluid(const T time,const T dt)
{
    LOG::SCOPE scope("Advect_Fluid");
    ARRAY<T_FACE_ARRAYS_SCALAR> face_velocities_ghost;
    ARRAY<T_ARRAYS_SCALAR> scalar_ghost;

    Fill_Ghost_Regions_Before_Advection(face_velocities_ghost,scalar_ghost,time,dt);
    Advect_Scalar_Fields(face_velocities_ghost,scalar_ghost,time,dt);
    if(!example.run_advection_test){
        Advect_Vector_Fields(face_velocities_ghost,time,dt);
        Revalidation();
        Modify_Levelset_And_Particles_After_Advection(face_velocities_ghost,time,dt);
        Coupling_Vector_Fields();}
    Coupling_Scalar_Fields(time,dt);
    if(!example.run_advection_test){
        Postprocess_After_Advection(time,dt);}
}
//#####################################################################
// Function Fill_Ghost_Regions_Before_Advection
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Fill_Ghost_Regions_Before_Advection(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,ARRAY<T_ARRAYS_SCALAR>& scalar_ghost,const T time,const T dt)
{
    LOG::SCOPE scope("Fill ghost cells before advection");

    int number_of_ghost_cells=example.fluids_parameters.number_of_ghost_cells;
    // fill velocities
    face_velocities_ghost.Resize(n_local_grids);
    FOR_EACH_LOCAL_GRID(grid_index)
        face_velocities_ghost(grid_index).Resize(LOCAL_GRID_ACCESS(grid_index).Domain_Indices(number_of_ghost_cells));
    ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_array(n_local_grids);
    ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_ghost_array(n_local_grids);
    for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
        face_velocities_array(grid_index)=&FACE_VELOCITIES(grid_index);
        face_velocities_ghost_array(grid_index)=&face_velocities_ghost(grid_index);}
    example.Fill_Ghost_Cells_Face_Chimera(number_of_ghost_cells,face_velocities_array,face_velocities_ghost_array);
    // fill vector fields
    if(example.solve_heat_equation_on_faces){
        ARRAY<T_FACE_ARRAYS_SCALAR*> face_temperatures_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            face_temperatures_array(grid_index)=&(example.face_temperatures(grid_index));}
        example.Fill_Ghost_Cells_Face_Chimera(number_of_ghost_cells,face_temperatures_array,face_temperatures_array);}
    // fill scalar fields
    scalar_ghost.Resize(n_local_grids);
    FOR_EACH_LOCAL_GRID(grid_index)
        scalar_ghost(grid_index).Resize(LOCAL_GRID_ACCESS(grid_index).Domain_Indices(number_of_ghost_cells));
    ARRAY<T_ARRAYS_SCALAR*> scalar_array(n_local_grids);
    ARRAY<T_ARRAYS_SCALAR*> scalar_ghost_array(n_local_grids);
    if(example.fluids_parameters.smoke){
        FOR_EACH_LOCAL_GRID(grid_index){
            scalar_array(grid_index)=&DENSITY(grid_index);
            scalar_ghost_array(grid_index)=&scalar_ghost(grid_index);}}
    if(example.solve_heat_equation){
        FOR_EACH_LOCAL_GRID(grid_index){
            scalar_array(grid_index)=&TEMPERATURE(grid_index);
            scalar_ghost_array(grid_index)=&scalar_ghost(grid_index);}}
    example.Fill_Ghost_Cells_Chimera(number_of_ghost_cells,scalar_array,scalar_ghost_array);
    // fill particles in ghost regions
    if(example.fluids_parameters.water){
        FOR_EACH_LOCAL_GRID(grid_index){example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Reseed_Particles_In_Ghost_Region(time,number_of_ghost_cells);}}
    if(example.run_advection_test) example.Fill_Ghost_Regions_Analytic(scalar_ghost,face_velocities_ghost,number_of_ghost_cells,time,dt);
    ////////////////////////////////////////////////////DEBUG
    /*FOR_EACH_LOCAL_GRID(grid_index){
        T_ARRAYS_SCALAR::Exchange_Arrays(DENSITY(grid_index),scalar_ghost(grid_index));
        T_FACE_ARRAYS_SCALAR::Exchange_Arrays(FACE_VELOCITIES(grid_index),face_velocities_ghost(grid_index));
    }
    Write_Substep("filled ghost cell velocity and density",0,2);
    FOR_EACH_LOCAL_GRID(grid_index){
        T_ARRAYS_SCALAR::Exchange_Arrays(DENSITY(grid_index),scalar_ghost(grid_index));
        T_FACE_ARRAYS_SCALAR::Exchange_Arrays(FACE_VELOCITIES(grid_index),face_velocities_ghost(grid_index));
    }*/
    /////////////////////////////////////////////////////
    Write_Substep("after filling ghost cells before advection",0,1);
}
//#####################################################################
// Function Advect_Scalar_Fields
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Advect_Scalar_Fields(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,ARRAY<T_ARRAYS_SCALAR>& scalar_ghost,const T time,const T dt)
{
    LOG::SCOPE scope("Advect scalar fields");

    int number_of_ghost_cells=example.fluids_parameters.number_of_ghost_cells;
    if(example.fluids_parameters.water){
        FOR_EACH_LOCAL_GRID(grid_index){example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Advance_Levelset(dt);}
        Write_Substep("after levelset advection",0,1);
        Revalidate_Fluid_Scalars();
        example.Extrapolate_Phi_Into_Objects(time+dt);
        example.chimera_grid->Adjustment_After_Advecting_Scalar_Fields();
        
        LOG::Time("advecting particles");
        FOR_EACH_LOCAL_GRID(grid_index){
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Euler_Step_Particles(face_velocities_ghost(grid_index),dt,time,true,true,false,0);}
        Write_Substep("after particle advection in the interior region",0,1);
        FOR_EACH_LOCAL_GRID(grid_index){
            T_GRID& grid=LOCAL_GRID_ACCESS(grid_index);
            PARTICLE_LEVELSET_UNIFORM<T_GRID>& pls=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset;
            pls.Advect_Particles_In_Ghost_Region(grid,face_velocities_ghost(grid_index),pls.positive_particles,dt,number_of_ghost_cells);
            pls.Advect_Particles_In_Ghost_Region(grid,face_velocities_ghost(grid_index),pls.negative_particles,dt,number_of_ghost_cells);}
        FOR_EACH_LOCAL_GRID(grid_index){
            Adjust_Particles_After_Updating_Grid(example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.positive_particles,LOCAL_RIGID_GRID_ACCESS(grid_index),dt,time);
            Adjust_Particles_After_Updating_Grid(example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.negative_particles,LOCAL_RIGID_GRID_ACCESS(grid_index),dt,time);}
        Write_Substep("after particle advection",0,1);
        //example.Scalar_Advection_Callback(dt,time);
            
        LOG::Time("updating removed particle velocities");
        //example.Modify_Removed_Particles_Before_Advection(dt,time);
        FOR_EACH_LOCAL_GRID(grid_index){
            RIGID_GRID<T_GRID>& rigid_grid=LOCAL_RIGID_GRID_ACCESS(grid_index);
            PARTICLE_LEVELSET_UNIFORM<T_GRID>& pls=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset;
            LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV> interpolation;
            if(pls.use_removed_positive_particles) for(NODE_ITERATOR iterator(rigid_grid.grid);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
                for(int p=1;p<=particles.array_collection->Size();p++){
                    TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(rigid_grid.grid,face_velocities_ghost(grid_index),X);
                    if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=rigid_grid.previous_state.Object_Space_Vector(example.fluids_parameters.removed_positive_particle_buoyancy_constant*example.fluids_parameters.gravity_direction); // buoyancy
                    particles.V(p)=V;}}
            if(pls.use_removed_negative_particles) for(NODE_ITERATOR iterator(rigid_grid.grid);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
                for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=rigid_grid.previous_state.Object_Space_Vector(dt*example.fluids_parameters.gravity*example.fluids_parameters.gravity_direction); // ballistic
            }}
        FOR_EACH_LOCAL_GRID(grid_index){
            PARTICLE_LEVELSET_UNIFORM<T_GRID>& pls=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset;
            if(pls.use_removed_positive_particles)
                Adjust_Removed_Particles_After_Updating_Grid(example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.removed_positive_particles,example.rigid_grid_collection.Rigid_Grid(grid_index),dt,time);
            if(pls.use_removed_negative_particles)
                Adjust_Removed_Particles_After_Updating_Grid(example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.removed_negative_particles,example.rigid_grid_collection.Rigid_Grid(grid_index),dt,time);}}
    else if(example.fluids_parameters.smoke){
        Write_Substep("before advecting density",0,1);
        ARRAY<T_ARRAYS_SCALAR*> density_array(n_local_grids);
        ARRAY<T_ARRAYS_SCALAR*> density_ghost_array(n_local_grids);
        FOR_EACH_LOCAL_GRID(grid_index){
            density_array(grid_index)=&DENSITY(grid_index);
            density_ghost_array(grid_index)=&scalar_ghost(grid_index);}
        ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_ghost_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            face_velocities_ghost_array(grid_index)=&face_velocities_ghost(grid_index);}
        example.Extrapolate_Scalar_Field_Into_Objects(density_ghost_array,dt,time);
        example.Advect_Scalar_Field(density_array,density_ghost_array,face_velocities_ghost_array,dt,time);
        example.chimera_grid->Adjustment_After_Advecting_Scalar_Fields();
        Write_Substep("after advecting density",0,1);}
    if(example.solve_heat_equation){
        ARRAY<T_ARRAYS_SCALAR*> scalar_array(n_local_grids);
        ARRAY<T_ARRAYS_SCALAR*> scalar_ghost_array(n_local_grids);
        FOR_EACH_LOCAL_GRID(grid_index){
            scalar_array(grid_index)=&TEMPERATURE(grid_index);
            scalar_ghost_array(grid_index)=&scalar_ghost(grid_index);}
        ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_ghost_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            face_velocities_ghost_array(grid_index)=&face_velocities_ghost(grid_index);}
        example.Advect_Scalar_Field(scalar_array,scalar_ghost_array,face_velocities_ghost_array,dt,time);
        example.chimera_grid->Adjustment_After_Advecting_Scalar_Fields();}
    Write_Substep("after scalar advection",0,1);
}
//#####################################################################
// Function Advect_Velocities
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Advect_Vector_Fields(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,const T time,const T dt)
{
    LOG::SCOPE scope("Advect velocities");

    LOG::Time("updating velocity (explicit part)"); // TODO: fix vorticity confinement for thin shells
    Write_Substep("before forces",1,1);
    if(example.solve_heat_equation_on_faces){
        ARRAY<T_FACE_ARRAYS_SCALAR*> vector_field_array(n_local_grids);
        ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_ghost_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            vector_field_array(grid_index)=&(example.face_temperatures(grid_index));
            face_velocities_ghost_array(grid_index)=&face_velocities_ghost(grid_index);}
            example.Advect_Vector_Field(vector_field_array,vector_field_array,face_velocities_ghost_array,dt,time);
    }else{
        ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_array(n_local_grids);
        ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_ghost_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            face_velocities_array(grid_index)=&FACE_VELOCITIES(grid_index);
            face_velocities_ghost_array(grid_index)=&face_velocities_ghost(grid_index);}
        example.Advect_Vector_Field(face_velocities_array,face_velocities_ghost_array,face_velocities_ghost_array,dt,time);}
    example.chimera_grid->Adjustment_After_Advecting_Velocities();
    Write_Substep("after advect velocity",1,1);
    //if(example.fluids_parameters.water)
    FOR_EACH_LOCAL_GRID(grid_index)
        Integrate_Fluid_Non_Advection_Forces(FACE_VELOCITIES(grid_index),LOCAL_RIGID_GRID_ACCESS(grid_index),dt,0);
    Write_Substep("after integrate non advection forces",1,1);
    //fluids_parameters.Blend_In_External_Velocity(face_velocities,dt,time);
    FOR_EACH_LOCAL_GRID(grid_index){
        Adjust_Face_Scalar_Field_After_Updating_Grid(face_velocities_ghost(grid_index),face_velocities_ghost(grid_index),LOCAL_RIGID_GRID_ACCESS(grid_index),dt,time);}
    Write_Substep("after updating velocity (explicit part)",0,1);
}
//#####################################################################
// Function Revalidation
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Revalidation()
{
    if(!example.use_collidable_advection) return;
    LOG::Time("effective velocity acceleration structures");
    // revalidate scalars and velocity in body's new position
    Write_Substep("before scalar revalidation",0,1);
    FOR_EACH_LOCAL_GRID(grid_index){
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
        example.Update_Object_Space_Rigid_Body_Collection_To_New_Grid_Frame(grid_index);
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*LOCAL_GRID_ACCESS(grid_index).Minimum_Edge_Length(),5);  // static occupied blocks
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Compute_Grid_Visibility();} // used in fast marching and extrapolation too... NOTE: this requires that objects_in_cell be current!
    Revalidate_Fluid_Scalars(); // uses visibility
    Revalidate_Fluid_Velocity(); // uses visibility
    Write_Substep("after scalar revalidation",0,1);
}
//#####################################################################
// Function Modify_Levelset_And_Particles_After_Advection
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Modify_Levelset_And_Particles_After_Advection(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,const T time,const T dt)
{
    /*if(example.fluids_parameters.water){
        LOG::Time("modifying levelset");
        example.Fill_Ghost_Cells_Chimera(2*example.fluids_parameters.number_of_ghost_cells+1);
        FOR_EACH_LOCAL_GRID(grid_index){
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Exchange_Overlap_Particles();}
        Write_Substep("exchanging overlap particles",0,1);
        FOR_EACH_LOCAL_GRID(grid_index){
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost(grid_index));}
        Write_Substep("after modify levelset and particles",0,1);
        Revalidate_Phi_After_Modify_Levelset(); // second revalidation -- uses visibility too
        Write_Substep("after revalidate phi",0,1);

        LOG::Time("adding sources");
        if(example.Adjust_Phi_With_Sources(time+dt)){
            example.Fill_Ghost_Cells_Chimera(2*example.fluids_parameters.number_of_ghost_cells+1);
            FOR_EACH_LOCAL_GRID(grid_index) example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Make_Signed_Distance();}
        example.Fill_Ghost_Cells_Chimera(example.fluids_parameters.number_of_ghost_cells);
        LOG::Time("getting sources");
        FOR_EACH_LOCAL_GRID(grid_index){
            T_ARRAYS_BOOL* source_mask=0;example.Get_Source_Reseed_Mask(source_mask,time+dt);
            if(source_mask){LOG::Time("reseeding sources");example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Reseed_Particles(time+dt,0,source_mask);delete source_mask;}}
        Write_Substep("after adding sources",0,1);
        
        FOR_EACH_LOCAL_GRID(grid_index){
            LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time+dt);
            LOG::Time("deleting particles in local maxima");
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);
            LOG::Time("deleting particles far from interface");
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Delete_Particles_Far_From_Interface(); // uses visibility
            Write_Substep("after delete particles far from interface",0,1);
            LOG::Time("re-incorporating removed particles");
            example.Modify_Removed_Particles_Before_Reincorporation(dt,time+dt);
            // TODO: if your particles fall entirely within the grid this shouldn't need ghost cells, but it may need them for MPI
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Identify_And_Remove_Escaped_Particles(face_velocities_ghost(grid_index),1.5,time+dt);
            if(example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.use_removed_positive_particles || example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.use_removed_negative_particles)
                example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Reincorporate_Removed_Particles(1,example.fluids_parameters.removed_particle_mass_scaling,example.fluids_parameters.reincorporate_removed_particle_velocity?&FACE_VELOCITIES(grid_index):0,!example.fluids_parameters.use_sph_for_removed_negative_particles || !example.fluids_parameters.sph_evolution->use_two_way_coupling);
            example.Modify_Removed_Particles_After_Reincorporation(dt,time+dt);}}*/
}
//#####################################################################
// Function Coupling_Scalar_Fields_After_Advection
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Coupling_Scalar_Fields(const T time,const T dt)
{
    LOG::SCOPE scope("Coupling scalar fields");
 
    ARRAY<T_ARRAYS_SCALAR*> scalar_array(n_local_grids);
    for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
        if(example.solve_heat_equation)
            scalar_array(grid_index)=&TEMPERATURE(grid_index);
        else
            scalar_array(grid_index)=&DENSITY(grid_index);}

    example.Coupling_Overlap_Regions_Cell(example.fluids_parameters.number_of_ghost_cells,scalar_array);
    Write_Substep("after scalar field injection",0,1);
    /*if(example.fluids_parameters.water){
        // make sign distance and reseed the injected region
        example.Fill_Ghost_Cells_Chimera(2*example.fluids_parameters.number_of_ghost_cells+1);
        FOR_EACH_LOCAL_GRID(grid_index){example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Make_Signed_Distance();}
        example.chimera_grid->Find_Overlap_Regions_Cell(0);
        FOR_EACH_LOCAL_GRID(grid_index){
            T_ARRAYS_BOOL overlap_cell_mask(LOCAL_GRID_ACCESS(grid_index).Domain_Indices(number_of_ghost_cells));
            overlap_cell_mask.Fill(false);
            for(int pack=1;pack<=example.chimera_grid->packages_cell_recv.Size();pack++){
                int send_grid=0,recv_grid=0,tag=example.chimera_grid->packages_cell_recv(pack).y;
                example.chimera_grid->Get_Send_Recv_Grid(tag,send_grid,recv_grid,n_global_grids*n_global_grids);
                if(grid_index==recv_grid){
                    for(int cell=1;cell<=(example.chimera_grid->packages_cell_recv(pack).x).x.Size();cell++){
                        overlap_cell_mask((example.chimera_grid->packages_cell_recv(pack).x).x(cell))=true;}}}
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Reseed_Delete_Particles_In_Whole_Region(number_of_ghost_cells,example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.positive_particles,1,&overlap_cell_mask);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Reseed_Delete_Particles_In_Whole_Region(number_of_ghost_cells,example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.negative_particles,-1,&overlap_cell_mask);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Reseed_Particles(time,&overlap_cell_mask);}
        Write_Substep("after reseeding in injection region",0,1);}*/
}
//#####################################################################
// Function Fill_Ghost_Regions_After_Advection
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Postprocess_After_Advection(const T time,const T dt)
{
    LOG::SCOPE scope("Fill ghost cells after advection"); 

    //decide reseeding
    if(example.fluids_parameters.water){
        example.particles_counter(1)=0;bool reseed=false;
        T_ARRAYS_PARTICLE_LEVELSET_PARTICLES* particles;
        FOR_EACH_LOCAL_GRID(grid_index) for(NODE_ITERATOR iterator(LOCAL_GRID_ACCESS(grid_index));iterator.Valid();iterator.Next()){
            PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>& particle_levelset_evolution=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution;
            particles=&(particle_levelset_evolution.particle_levelset.positive_particles);
            TV_INT block=iterator.Node_Index();
            if((*particles)(block)){
                typename T_ARRAYS_PARTICLE_LEVELSET_PARTICLES::ELEMENT cell_particles=(*particles)(block);
                while(cell_particles){
                    example.particles_counter(1)+=cell_particles->array_collection->Size();
                    cell_particles=cell_particles->next;}}
            particles=&(particle_levelset_evolution.particle_levelset.negative_particles);
            if((*particles)(block)){
                typename T_ARRAYS_PARTICLE_LEVELSET_PARTICLES::ELEMENT cell_particles=(*particles)(block);
                while(cell_particles){
                    example.particles_counter(1)+=cell_particles->array_collection->Size();
                    cell_particles=cell_particles->next;}}}
        if(example.particles_counter(2)<example.particles_counter(1)) example.particles_counter(2)=example.particles_counter(1);
        else if(example.particles_counter(1)<example.particles_counter(2)/3*2) reseed=true;
        if(reseed){FOR_EACH_LOCAL_GRID(grid_index){
            PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>& particle_levelset_evolution=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution;
            particle_levelset_evolution.Reseed_Particles(time);
            particle_levelset_evolution.Delete_Particles_Outside_Grid();}
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);
            example.particles_counter(2)=0;}}
}
//#####################################################################
// Function Coupling_Vector_Fields_After_Projection
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Coupling_Vector_Fields()
{
    LOG::SCOPE scope("Coupling vector fields");

    ARRAY<T_FACE_ARRAYS_SCALAR*> vector_field_array(n_local_grids);
    for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
        if(example.solve_heat_equation_on_faces || example.solve_heat_equation_on_faces_coupled) vector_field_array(grid_index)=&(example.face_temperatures(grid_index));
        else vector_field_array(grid_index)=&FACE_VELOCITIES(grid_index);}
    
    Write_Substep("before coupling vector fields",0,1);

    example.Coupling_Overlap_Regions_Face(example.fluids_parameters.number_of_ghost_cells,vector_field_array);
  
    Write_Substep("after coupling vector fields",0,1);
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Revalidate_Fluid_Scalars()
{
    if(!example.simulate_solids || !example.use_collidable_advection) return;
    FOR_EACH_LOCAL_GRID(grid_index){
        if(example.fluids_parameters.water){
            T_FAST_LEVELSET& levelset=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset;
            int sign=1;//if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
            example.advection_collidable(grid_index)->Average_To_Invalidated_Cells(LOCAL_GRID_ACCESS(grid_index),sign*example.fluids_parameters.collidable_phi_replacement_value,levelset.phi);}
        else if(example.fluids_parameters.smoke){
            example.advection_collidable(grid_index)->Average_To_Invalidated_Cells(LOCAL_GRID_ACCESS(grid_index),(T)0,example.incompressible_fluid_containers(grid_index)->density_container.density);}}
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Revalidate_Phi_After_Modify_Levelset()
{
    if(!example.simulate_solids || !example.use_collidable_advection) return;
    FOR_EACH_LOCAL_GRID(grid_index){
        if(example.fluids_parameters.water){
            T_FAST_LEVELSET& levelset=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset;
            int sign=1;//if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
            example.advection_collidable(grid_index)->cell_valid_points_current=example.advection_collidable(grid_index)->cell_valid_points_next;
            example.advection_collidable(grid_index)->Average_To_Invalidated_Cells(LOCAL_GRID_ACCESS(grid_index),sign*example.fluids_parameters.collidable_phi_replacement_value,levelset.phi);}}
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Revalidate_Fluid_Velocity()
{
    if(!example.simulate_solids || !example.use_collidable_advection) return;
    FOR_EACH_LOCAL_GRID(grid_index){
        example.advection_collidable(grid_index)->Average_To_Invalidated_Face(example.rigid_grid_collection.Rigid_Grid(grid_index).grid,FACE_VELOCITIES(grid_index));}
}
//#####################################################################
// Function Adjust_Particles_After_Updating_Grid
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Adjust_Particles_After_Updating_Grid(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,RIGID_GRID<T_GRID>& rigid_grid,const T dt,const T time)
{
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
    for(NODE_ITERATOR iterator(rigid_grid.grid);iterator.Valid();iterator.Next()){
        TV_INT block=iterator.Node_Index();
        if(particles(block)){
            typename T_ARRAYS_PARTICLE_LEVELSET_PARTICLES::ELEMENT cell_particles=particles(block);
            while(cell_particles){
                for(int k=1;k<=cell_particles->array_collection->Size();k++)
                    cell_particles->X(k)=rigid_grid.Current_Object_Space_Point(cell_particles->X(k));
                cell_particles=cell_particles->next;}}}
}
//#####################################################################
// Function Adjust_Removed_Particles_After_Updating_Grid
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Adjust_Removed_Particles_After_Updating_Grid(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,RIGID_GRID<T_GRID>& rigid_grid,const T dt,const T time)
{
    for(NODE_ITERATOR iterator(rigid_grid.grid);iterator.Valid();iterator.Next()){
        TV_INT block=iterator.Node_Index();
        if(particles(block)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& cell_particles=*particles(block);
            for(int k=cell_particles.array_collection->Size();k>=1;k--){
                cell_particles.X(k)=rigid_grid.Current_Object_Space_Point(cell_particles.X(k));
                cell_particles.V(k)=rigid_grid.Current_Object_Space_Vector(cell_particles.V(k));}}}
}
//#####################################################################
// Function Adjust_Face_Scalar_Field_After_Updating_Grid
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Adjust_Face_Scalar_Field_After_Updating_Grid(T_FACE_ARRAYS_SCALAR& u, T_FACE_ARRAYS_SCALAR u_ghost, RIGID_GRID<T_GRID>& rigid_grid, const T dt, const T time)
{
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
    for(FACE_ITERATOR iterator(rigid_grid.grid,1);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
        X=rigid_grid.Previous_Object_Space_Point(X);
        TV u_tmp;
        for(int axis_local=1;axis_local<=TV::dimension;axis_local++)
            u_tmp(axis_local)=interpolation.Clamped_To_Array_Face_Component(axis_local,rigid_grid.grid,u_ghost,X);
        u_tmp=rigid_grid.Current_Object_Space_Vector(u_tmp);
        u.Component(axis)(face)=u_tmp(axis);}
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Extrapolate_Velocity_Across_Interface(T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_ARRAYS_SCALAR& phi_ghost,const T ghost,const T band_width,const T damping,const TV& air_speed,const T_FACE_ARRAYS_BOOL_DIMENSION* face_neighbors_visible,const T_FACE_ARRAYS_BOOL* fixed_faces_input)
{
    int current_grid=example.chimera_grid->Find_Current_Grid_Local_Index(grid);
    T_FACE_ARRAYS_BOOL& above_water_face_mask=above_water_face_masks(current_grid);
    above_water_face_mask.Resize(grid.Domain_Indices((int)ghost));above_water_face_mask.Fill(true);
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        T_GRID face_grid=grid.Get_Face_Grid(axis);T_ARRAYS_SCALAR phi_face(face_grid.Domain_Indices((int)ghost));T_ARRAYS_BASE& face_velocity=face_velocities.Component(axis);
        T_ARRAYS_BOOL fixed_face;T_ARRAYS_BASE_BOOL& above_water_face=above_water_face_mask.Component(axis);
        if(fixed_faces_input) fixed_face=fixed_faces_input->Component(axis); else fixed_face=T_ARRAYS_BOOL(face_grid.Domain_Indices((int)ghost));
        for(FACE_ITERATOR iterator(grid,(int)ghost,T_GRID::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Face_Index();
            int number_valid_neighbor_cell=0;
            if(grid.Domain_Indices((int)ghost).Lazy_Inside(iterator.First_Cell_Index())){phi_face(index)+=phi_ghost(iterator.First_Cell_Index());++number_valid_neighbor_cell;}
            if(grid.Domain_Indices((int)ghost).Lazy_Inside(iterator.Second_Cell_Index())){phi_face(index)+=phi_ghost(iterator.Second_Cell_Index());++number_valid_neighbor_cell;}
            if(number_valid_neighbor_cell>1) phi_face(index)/=number_valid_neighbor_cell;
            if((grid.Domain_Indices((int)ghost).Lazy_Inside(iterator.First_Cell_Index()) && phi_ghost(iterator.First_Cell_Index())<=0) ||
                (grid.Domain_Indices((int)ghost).Lazy_Inside(iterator.Second_Cell_Index()) && phi_ghost(iterator.Second_Cell_Index())<=0)){fixed_face(index)=true;above_water_face(index)=false;}}
        T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,(int)ghost);extrapolate.Set_Band_Width((T)band_width);extrapolate.Set_Custom_Seed_Done(&fixed_face);
        if(face_neighbors_visible) extrapolate.Set_Collision_Aware_Extrapolation(face_neighbors_visible->Component(axis));
        extrapolate.Extrapolate(0,false);

/*        if(damping) for(FACE_ITERATOR iterator(grid,0,T_GRID::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Face_Index();
            if(!fixed_face(index) && phi_face(index)<delta) face_velocity(index)=(1-damping)*face_velocity(index)+damping*air_speed[axis];}*/
    }
}
//#####################################################################
// Function Integrate_Fluid_Non_Advection_Forces
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Integrate_Fluid_Non_Advection_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const RIGID_GRID<T_GRID>& rigid_grid,const T dt,const int substep)
{
    //FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters=example.fluids_parameters;
    //int number_of_regions=fluids_parameters.number_of_regions;
    //PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    //PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>* particle_levelset_evolution_multiple=fluids_parameters.particle_levelset_evolution_multiple;
    //INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>* incompressible_multiphase=fluids_parameters.incompressible_multiphase;
    //INCOMPRESSIBLE_UNIFORM<T_GRID>* incompressible=fluids_parameters.incompressible;
    //EULER_UNIFORM<T_GRID>* euler=fluids_parameters.euler;
    //SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>* euler_solid_fluid_coupling_utilities=fluids_parameters.euler_solid_fluid_coupling_utilities;

    LOG::SCOPE integration_scope("INTEGRATE NON ADVECTION FORCES","integrate non advection forces");

    /*if(fluids_parameters.use_body_force){
        if(Simulate_Incompressible_Fluids()) fluids_parameters.callbacks->Get_Body_Force(fluids_parameters.incompressible->force,dt,time);
        if(fluids_parameters.compressible) fluids_parameters.callbacks->Get_Body_Force(fluids_parameters.euler->force,dt,time);}*/

    LOG::Time("updating velocity (explicit part without convection)"); // TODO: fix vorticity confinement for thin shells
    //if(fluids_parameters.analytic_test){time=particle_levelset_evolution.time;example.Get_Analytic_Velocities(particle_levelset_evolution.time+dt);}
    //else{
    /*if(fluids_parameters.compressible){
            euler->Save_State(euler->U_save,euler->euler_projection.face_velocities_save,euler->need_to_remove_added_internal_energy_save);
            euler->Advance_One_Time_Step_Forces(dt,time);
            euler->Fill_Ghost_Cells(dt,time,example.fluids_parameters.number_of_ghost_cells);
            if(!euler->timesplit || !euler->thinshell) euler_solid_fluid_coupling_utilities->Fill_Solid_Cells(); // TODO(kwatra): see if can get rid of this, since one-sided interpolation is used in slip
            euler->Get_Dirichlet_Boundary_Conditions(dt,time);
            if(euler->timesplit) euler->euler_projection.Get_Pressure(euler->euler_projection.p_advected);}*/
    //if(Simulate_Incompressible_Fluids()){
    //if(fluids_parameters.fluid_affects_solid) incompressible->projection.Set_Up_For_Projection(face_velocities);
    // TODO(kwatra): Check if SPH case is handled properly.
    //Write_Substep("before viscosity",substep,1);
    //if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(fluids_parameters.projection))
    //    coupled_evolution->Apply_Viscosity(face_velocities,dt,time);
    //Write_Substep("after viscosity",substep,1);
    //if(number_of_regions>=2)
    //    incompressible_multiphase->Advance_One_Time_Step_Forces(face_velocities,dt,time,fluids_parameters.implicit_viscosity,&particle_levelset_evolution_multiple->phis,
    //        &fluids_parameters.pseudo_dirichlet_regions,fluids_parameters.number_of_ghost_cells);
    // else if(!fluids_parameters.sph) 
    //Update gravity
    //if(example.fluids_parameters.water){
        TV gravity_direction_acceleration;gravity_direction_acceleration+=dt*example.fluids_parameters.gravity*example.fluids_parameters.gravity_direction;
        TV object_space_acceleration=rigid_grid.Object_Space_Vector(gravity_direction_acceleration);
        for(typename GRID<TV>::FACE_ITERATOR iterator(rigid_grid.grid);iterator.Valid();iterator.Next())
            face_velocities.Component(iterator.Axis())(iterator.Face_Index())+=object_space_acceleration(iterator.Axis());
    //incompressible->Advance_One_Time_Step_Forces(face_velocities,dt,time,fluids_parameters.implicit_viscosity,&particle_levelset_evolution.phi,fluids_parameters.number_of_ghost_cells);//}
    //fluids_parameters.Blend_In_External_Velocity(face_velocities,dt,time);}
}
//#####################################################################
// Function Solid_Position_Update
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Solid_Position_Update(const T time,const T dt)
{
    LOG::SCOPE scope("solids position update");
    example.Set_External_Positions(example.solid_body_collection.rigid_body_collection.rigid_body_particle.X,example.solid_body_collection.rigid_body_collection.rigid_body_particle.rotation,time+dt);
}
//#####################################################################
// Function Solid_Velocity_Update
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Solid_Velocity_Update(const T time,const T dt)
{
    LOG::SCOPE scope("solids velocity update");

    example.Set_External_Velocities(example.solid_body_collection.rigid_body_collection.rigid_body_particle.V,example.solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity,time+dt,time+dt);
}
//#####################################################################
// Function Advance_Fluid_One_Time_Step_Implicit_Part
//#####################################################################
template<class T_GRID> class LAPLACE_CALLBACKS_INCOMPRESSIBLE:public LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
public:
    SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution;
    
    LAPLACE_CALLBACKS_INCOMPRESSIBLE(SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution_input):
        evolution(evolution_input){}
    
    bool Get_Neumann_Boundary_Condition(const TV& location,const TV& normal,T& normal_gradient,const T time)
    {return evolution.example.Get_Source_Velocity(location,normal,normal_gradient,time,0);}
    bool Get_Dirichlet_Boundary_Condition(const TV& location,T& value,const T time)
    {return evolution.example.Get_Dirichlet_Boundary_Condition(location,value,time);}
    T Get_Density(const TV& location)
    {return evolution.example.Get_Density(location);}
};

template<class T_GRID> class LAPLACE_CALLBACKS_VISCOSITY:public LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
public:
    SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution;
    
    TV direction;
    bool define_slip_boundary;
    RANGE<TV> background_domain;

    LAPLACE_CALLBACKS_VISCOSITY(SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution_input):
        evolution(evolution_input)
    {
        define_slip_boundary=evolution.example.define_slip_boundary;
        background_domain=evolution.example.chimera_grid->initial_grids(1).domain;
    }
    
    bool Get_Neumann_Boundary_Condition(const TV& location,const TV& normal,T& normal_gradient,const T time)
    {
        if(define_slip_boundary && (location(2)>=background_domain.max_corner(2) || location(2)<=background_domain.min_corner(2)))
            return evolution.example.Get_Source_Velocity(location,normal,normal_gradient,time,0);

        normal_gradient=(T)0;
        return evolution.example.Get_Dirichlet_Boundary_Condition(location,normal_gradient,time);
    }
    bool Get_Dirichlet_Boundary_Condition(const TV& location,T& value,const T time)
    {
        if(define_slip_boundary && (location(2)>background_domain.max_corner(2) || location(2)<background_domain.min_corner(2)))
            return false;

        return evolution.example.Get_Source_Velocity(location,direction,value,time,0);
    }
    T Get_Density(const TV& location)
    {return evolution.example.Get_Density(location)/evolution.example.fluids_parameters.viscosity;}
};

template<class T_GRID> class LAPLACE_CALLBACKS_DIFFUSION_FACE:public LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
public:
    SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution;
    
    LAPLACE_CALLBACKS_DIFFUSION_FACE(SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution_input):
        evolution(evolution_input){}
    
    bool Get_Neumann_Boundary_Condition(const TV& location,const TV& normal,T& normal_gradient,const T time) //these are dirichlet conditions for faces
    {
        TV vector;
        for(int axis=1;axis<=TV::dimension;axis++){
            evolution.example.test_number=evolution.example.test_numbers(axis);
            if(!evolution.example.Get_Dirichlet_Boundary_Condition(location,vector(axis),time))
                return false;}
        //LOG::cout << "dirichlet vector " << location << " " << vector << std::endl;
        normal_gradient=TV::Dot_Product(vector,normal);
        return true;
    }
    bool Get_Dirichlet_Boundary_Condition(const TV& location,T& value,const T time)
    {return evolution.example.Get_Dirichlet_Boundary_Condition(location,value,time);}
    T Get_Density(const TV& location)
    {return evolution.example.Get_Density(location);}
};

#define WRITE_CELL_ARRAY(X,S) {                                         \
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++) if(laplace_grid.Local_Grid(grid_index)) \
            T_ARRAYS_SCALAR::Exchange_Arrays(example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->pressure,X(grid_index)); \
        Write_Substep(S,0,1);                                           \
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++) if(laplace_grid.Local_Grid(grid_index)) \
            T_ARRAYS_SCALAR::Exchange_Arrays(example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->pressure,X(grid_index));}

#define WRITE_CELL_ARRAY_POINTER(X,S) WRITE_CELL_ARRAY(*X,S)

#define WRITE_CELL_VECTOR(L,X,S) {                      \
        ARRAY<T_ARRAYS_SCALAR> debug_x(n_global_grids); \
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++) if(laplace_grid.Local_Grid(grid_index)){ \
            debug_x(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1)); \
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next()){ \
                int matrix_cell_index=L.Matrix_Cell_Index(grid_index,iterator.Cell_Index()); \
                debug_x(grid_index)(iterator.Cell_Index())=matrix_cell_index?X(matrix_cell_index):0;}} \
        WRITE_CELL_ARRAY(debug_x,S);}

#define WRITE_FACE_ARRAY(X,S) {                                         \
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++) if(laplace_grid.Local_Grid(grid_index)) \
            T_FACE_ARRAYS_SCALAR::Exchange_Arrays(example.solve_heat_equation?example.face_temperatures(grid_index):example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->face_velocities,X(grid_index)); \
        Write_Substep(S,0,1);                                           \
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++) if(laplace_grid.Local_Grid(grid_index)) \
            T_FACE_ARRAYS_SCALAR::Exchange_Arrays(example.solve_heat_equation?example.face_temperatures(grid_index):example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->face_velocities,X(grid_index));}

#define WRITE_FACE_VECTOR(L,X,S) {                      \
        ARRAY<T_FACE_ARRAYS_SCALAR> debug_x(n_global_grids); \
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++) if(laplace_grid.Local_Grid(grid_index)){ \
            debug_x(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1)); \
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next()){ \
                int matrix_face_index=L.Matrix_Face_Index(grid_index,iterator.Full_Index()); \
                debug_x(grid_index)(iterator.Full_Index())=matrix_face_index?X(matrix_face_index):0;}} \
        WRITE_FACE_ARRAY(debug_x,S);}

template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Advance_Fluid_One_Time_Step_Implicit_Part(const typename T_GRID::SCALAR time,const typename T_GRID::SCALAR dt)
{
    LOG::SCOPE scope("Advance_Fluid_One_Time_Step_Implicit_Part");
    
    laplace_grid.Construct_Mesh();
    Write_Substep("after constructing laplace grid",0,2);
    LAPLACE_CALLBACKS_INCOMPRESSIBLE<T_GRID> incompressible_callbacks(*this);
    LAPLACE_CALLBACKS_VISCOSITY<T_GRID> viscosity_callbacks(*this);
    Write_Substep("after constructing laplace incompressible",0,2);
    
    if(example.solve_poisson_equation){
        laplace_incompressible.Construct_System(incompressible_callbacks,time,dt,false);
        Build_Poisson_System(laplace_incompressible,incompressible_callbacks,time,dt,0,x_incompressible,rhs_incompressible);
        Solve_System(laplace_incompressible.system_matrix,rhs_incompressible,x_incompressible,example.fluids_parameters.incompressible_tolerance,example.fluids_parameters.incompressible_iterations,&laplace_incompressible.partition);
        ARRAY<T_ARRAYS_SCALAR*> pressures_pointer(n_global_grids);
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
            pressures_pointer(grid_index)=laplace_grid.Local_Grid(grid_index)?&example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->pressure:0;
        Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_incompressible,incompressible_callbacks,time+dt,pressures_pointer,pressures_boundary,&x_incompressible);}
    else if(example.solve_heat_equation){
        VECTOR_ND<T> rhs_diffusion,x_diffusion;
        //LOG::cout << "diffusion system matrix " << laplace_diffusion.system_matrix << std::endl;
        
        if(example.solve_heat_equation_on_faces){
            Write_Substep("before face heat solve",0,2);
            LAPLACE_CALLBACKS_DIFFUSION_FACE<T_GRID> callbacks_face(*this);
            ARRAY<T_FACE_ARRAYS_SCALAR*> face_temperature_pointer(n_global_grids);
            for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                face_temperature_pointer(grid_index)=laplace_grid.Local_Grid(grid_index)?&example.face_temperatures(grid_index):0;
            
            VECTOR<ARRAY<T_ARRAYS_SCALAR>,TV::dimension> cell_temperature_components,cell_temperature_components_updated;
            VECTOR<ARRAY<T_ARRAYS_SCALAR*>,TV::dimension> cell_temperature_components_pointer,cell_temperature_components_updated_pointer;
            VECTOR<ARRAY<ARRAY<T> >,TV::dimension> cell_temperature_boundary_components,cell_temperature_boundary_components_updated;
            for(int axis=1;axis<=TV::dimension;axis++){
                cell_temperature_components(axis).Resize(n_global_grids);
                cell_temperature_components_updated(axis).Resize(n_global_grids);
                cell_temperature_components_pointer(axis).Resize(n_global_grids);
                cell_temperature_components_updated_pointer(axis).Resize(n_global_grids);
                for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
                    cell_temperature_components_pointer(axis)(grid_index)=&cell_temperature_components(axis)(grid_index);
                    cell_temperature_components_updated_pointer(axis)(grid_index)=&cell_temperature_components_updated(axis)(grid_index);
                    if(laplace_grid.Local_Grid(grid_index)){
                        cell_temperature_components(axis)(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1));
                        cell_temperature_components_updated(axis)(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1));
                        for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                            T value=0;
                            if(callbacks_face.Get_Neumann_Boundary_Condition(laplace_grid.Frame(grid_index)*iterator.Location(),laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(iterator.Axis())),value,time))
                                (*face_temperature_pointer(grid_index))(iterator.Full_Index())=value;}}}}

            Write_Substep("after loading boundary conditions",0,2);            

            //solve system, build it if coupled too
            if(example.solve_heat_equation_on_faces_coupled){
                //we need to construct this just to get the face indices, if we want large system we'll need conserve memory and do something else
                laplace_diffusion.Construct_System(callbacks_face,time,example.fluids_parameters.use_trapezoid_rule?((T).5*dt):dt,true);
                
                if(!face_temperature_packed.Size())
                    face_temperature_packed.Resize(laplace_diffusion.n_matrix_faces);
                cell_temperature_packed.Resize(laplace_diffusion.n_matrix_cells*TV::dimension);
                if(time==0)
                    example.Initialize_Velocities();
                //Pack_And_Compute_Voronoi_Face_Values(laplace_diffusion,face_temperature_pointer,face_temperature_packed,time,dt);
                
                averaging.Build_Face_To_Cell_Matrix(laplace_diffusion,callbacks_face,time);
                averaging.Build_Face_To_Cell_Vector(laplace_diffusion,callbacks_face,time);
                averaging.face_to_cell_matrix.Times(face_temperature_packed,cell_temperature_packed);
                cell_temperature_packed+=averaging.face_to_cell_vector;

                for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                    if(laplace_grid.Local_Grid(grid_index))
                        for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                            int matrix_cell_index=laplace_diffusion.Matrix_Cell_Index(grid_index,iterator.Cell_Index());
                            if(matrix_cell_index)
                                for(int axis=1;axis<=TV::dimension;axis++)
                                    cell_temperature_components(axis)(grid_index)(iterator.Cell_Index())=cell_temperature_packed(laplace_diffusion.n_matrix_cells*(axis-1)+matrix_cell_index);}

                VECTOR<VECTOR_ND<T>,TV::dimension> rhs_diffusion_components,x_diffusion_components;
                for(int axis=1;axis<=TV::dimension;axis++){//we need this step to set the unsolved values to the interpolated values in order for flip to difference correctly
                    example.test_number=example.test_numbers(axis); //the callbacks should be fixed to do the right thing for neumann/dirichlet boundary conditions
                    laplace_diffusion_components(axis)->Construct_System(incompressible_callbacks,time,example.fluids_parameters.use_trapezoid_rule?((T).5*dt):dt,false);
                    Unpack_And_Interpolate_And_Inject_Cell_Values(*laplace_diffusion_components(axis),incompressible_callbacks,time,cell_temperature_components_pointer(axis),cell_temperature_boundary_components(axis),0);
                    Build_Diffusion_System(*laplace_diffusion_components(axis),incompressible_callbacks,time,dt,cell_temperature_components_pointer(axis),cell_temperature_boundary_components(axis),x_diffusion_components(axis),rhs_diffusion_components(axis),true);}

                for(int axis=1;axis<=TV::dimension;axis++)                
                    WRITE_CELL_ARRAY(cell_temperature_components(axis),STRING_UTILITIES::string_sprintf("cell vector component %d",axis));

                //construct the matrix
                ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > blocks(TV::dimension);
                for(int axis=1;axis<=TV::dimension;axis++) blocks(axis)=laplace_diffusion_components(axis)->system_matrix;
                SPARSE_MATRIX_FLAT_MXN<T> block_diagonal_negative_cell_laplacians=SPARSE_MATRIX_FLAT_MXN<T>::Concatenate_Matrices(blocks,true);
                SPARSE_MATRIX_FLAT_MXN<T> averaging_transpose;
                averaging.face_to_cell_matrix.Transpose(averaging_transpose);
                SPARSE_MATRIX_FLAT_MXN<T> system_matrix=averaging_transpose*block_diagonal_negative_cell_laplacians*averaging.face_to_cell_matrix;
                system_matrix*=example.fluids_parameters.use_trapezoid_rule?(T).5*dt:dt;
                for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                    if(laplace_grid.Local_Grid(grid_index))
                        for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                            D_FACE_INDEX face_index=iterator.Full_Index();
                            int matrix_face_index=laplace_diffusion.Matrix_Face_Index(grid_index,face_index);
                            if(matrix_face_index)
                                system_matrix.Add_Element(matrix_face_index,matrix_face_index,laplace_grid.Face_Size(grid_index,face_index));}
                for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
                    int matrix_face_index=laplace_diffusion.Matrix_Face_Index(voronoi_face_index);
                    if(matrix_face_index)
                        system_matrix.Add_Element(matrix_face_index,matrix_face_index,laplace_grid.Face_Size(voronoi_face_index));}
                SPARSE_MATRIX_FLAT_NXN<T> A=system_matrix.Create_NXN_Matrix();
                
                //OCTAVE_OUTPUT<T> oo("system.mat");
                //oo.Write("A",A);
                
                VECTOR_ND<T> cell_explicit_rhs(laplace_diffusion.n_matrix_cells*TV::dimension);
                averaging.Build_Face_To_Cell_Vector(laplace_diffusion,callbacks_face,time+dt);
                cell_explicit_rhs=averaging.face_to_cell_vector;
                VECTOR_ND<T> cell_rhs(laplace_diffusion.n_matrix_cells*TV::dimension);
                block_diagonal_negative_cell_laplacians.Times(cell_explicit_rhs,cell_rhs);
                cell_rhs*=example.fluids_parameters.use_trapezoid_rule?((T)-.5*dt):(T)(-dt);
                for(int axis=1;axis<=TV::dimension;axis++)
                    cell_rhs.Add_Subvector((axis-1)*laplace_diffusion.n_matrix_cells+1,rhs_diffusion_components(axis));
                VECTOR_ND<T> rhs(laplace_diffusion.n_matrix_faces);
                averaging_transpose.Times(cell_rhs,rhs);
                
                VECTOR_ND<T> volume_weighted_face_temperature(laplace_diffusion.n_matrix_faces);
                laplace_diffusion.Multiply_By_Dual_Cell_Size(face_temperature_packed,volume_weighted_face_temperature);
                rhs+=volume_weighted_face_temperature;
                //WRITE_FACE_VECTOR(laplace_diffusion,rhs,"face rhs");
                
                VECTOR_ND<T> x;
                Solve_System(A,rhs,x,example.fluids_parameters.implicit_viscosity_tolerance,example.fluids_parameters.implicit_viscosity_iterations,0);
                //WRITE_FACE_VECTOR(laplace_diffusion,x,"face x");
                
                averaging.face_to_cell_matrix.Times(x,cell_temperature_packed);
                cell_temperature_packed+=averaging.face_to_cell_vector;
                for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                    if(laplace_grid.Local_Grid(grid_index))
                        for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                            int matrix_face_index=laplace_diffusion.Matrix_Face_Index(grid_index,iterator.Full_Index());
                            if(matrix_face_index)
                                (*face_temperature_pointer(grid_index))(iterator.Full_Index())=x(matrix_face_index);}
                face_temperature_packed=x;

                //fill solved voronoi faces here
                for(int axis=1;axis<=TV::dimension;axis++){
                    x_diffusion_components(axis).Resize(laplace_diffusion.n_matrix_cells);
                    cell_temperature_packed.Get_Subvector((axis-1)*laplace_diffusion.n_matrix_cells+1,x_diffusion_components(axis));
                    WRITE_CELL_VECTOR((*laplace_diffusion_components(axis)),x_diffusion_components(axis),STRING_UTILITIES::string_sprintf("cell vector component solved %d",axis));}

                for(int axis=1;axis<=TV::dimension;axis++){
                    example.test_number=example.test_numbers(axis);
                    Unpack_And_Interpolate_And_Inject_Cell_Values(*laplace_diffusion_components(axis),incompressible_callbacks,time+dt,
                        cell_temperature_components_updated_pointer(axis),cell_temperature_boundary_components_updated(axis),
                        &x_diffusion_components(axis));}}
            else{
                for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                    if(laplace_grid.Local_Grid(grid_index)){
                        for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                            TV temperature_vector;
                            for(int axis=1;axis<=TV::dimension;axis++)
                                temperature_vector(axis)=(T).5*(face_temperature_pointer(grid_index)->Component(axis)(iterator.Cell_Index())+face_temperature_pointer(grid_index)->Component(axis)(iterator.Cell_Index()+TV_INT::Axis_Vector(axis)));
                            temperature_vector=laplace_grid.Frame(grid_index).r.Rotate(temperature_vector);
                            for(int axis=1;axis<=TV::dimension;axis++)
                                cell_temperature_components(axis)(grid_index)(iterator.Cell_Index())=temperature_vector(axis);}}

                for(int axis=1;axis<=TV::dimension;axis++){//we need this step to set the unsolved values to the interpolated values in order for flip to difference correctly
                    example.test_number=example.test_numbers(axis); //the callbacks should be fixed to do the right thing for neumann/dirichlet boundary conditions
                    laplace_diffusion.Construct_System(incompressible_callbacks,time,example.fluids_parameters.use_trapezoid_rule?((T).5*dt):dt,true);

                    if(example.use_flip_face_update_consistent)
                        Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_diffusion,incompressible_callbacks,time,cell_temperature_components_pointer(axis),cell_temperature_boundary_components(axis),0);
                    else{
                        LOG::cout << "using inconsistent flip update" << std::endl;
                        Set_Dirichlet_Boundary_Conditions(incompressible_callbacks,time,cell_temperature_components_pointer(axis));
                        Pack_And_Exchange_Boundary_Cell_Values(cell_temperature_components_pointer(axis),cell_temperature_boundary_components(axis));
                        ARRAY<T_ARRAYS_SCALAR*> cell_values_tmp(n_local_grids);
                        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
                            cell_values_tmp(grid_index)=cell_temperature_components_pointer(axis)(laplace_grid.Global_Grid_Index(grid_index));}
                        example.Fill_Ghost_Cells_Chimera(1,cell_values_tmp,cell_values_tmp,false);}
                    Build_Diffusion_System(laplace_diffusion,incompressible_callbacks,time,dt,cell_temperature_components_pointer(axis),cell_temperature_boundary_components(axis),x_diffusion,rhs_diffusion,false);
                    Solve_System(laplace_diffusion.system_matrix,rhs_diffusion,x_diffusion,example.fluids_parameters.implicit_viscosity_tolerance,example.fluids_parameters.implicit_viscosity_iterations,&laplace_diffusion.partition);
                    WRITE_CELL_VECTOR(laplace_diffusion,x_diffusion,"x solution");
                    Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_diffusion,incompressible_callbacks,time+dt,cell_temperature_components_updated_pointer(axis),cell_temperature_boundary_components_updated(axis),&x_diffusion);}}
            
            LOG::cout << "after solving all systems" << std::endl;
            
            for(int axis=1;axis<=TV::dimension;axis++){
                WRITE_CELL_ARRAY_POINTER(cell_temperature_components_pointer(axis),STRING_UTILITIES::string_sprintf("old values %d",axis));
                WRITE_CELL_ARRAY_POINTER(cell_temperature_components_updated_pointer(axis),STRING_UTILITIES::string_sprintf("new values %d",axis));}
            for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                if(laplace_grid.Local_Grid(grid_index))
                    for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next())
                        if(!example.solve_heat_equation_on_faces_coupled || !laplace_diffusion.Matrix_Face_Index(grid_index,iterator.Full_Index())){
                            T value=0;
                            if(callbacks_face.Get_Neumann_Boundary_Condition(laplace_grid.Frame(grid_index)*iterator.Location(),laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(iterator.Axis())),value,time+dt))
                                (*face_temperature_pointer(grid_index))(iterator.Full_Index())=value;
                            else{
                                TV temperature_vector;
                                bool flip=!example.solve_heat_equation_on_faces_coupled && example.use_flip_face_update && (example.use_flip_update_overlap || (laplace_grid.Chimera_Cell(grid_index,iterator.First_Cell_Index()) || laplace_grid.Chimera_Cell(grid_index,iterator.Second_Cell_Index())));
                                for(int axis=1;axis<=TV::dimension;axis++){
                                    temperature_vector(axis)=(T).5*(cell_temperature_components_updated(axis)(grid_index)(iterator.First_Cell_Index())+cell_temperature_components_updated(axis)(grid_index)(iterator.Second_Cell_Index()));
                                    if(flip) temperature_vector(axis)-=(T).5*(cell_temperature_components(axis)(grid_index)(iterator.First_Cell_Index())+cell_temperature_components(axis)(grid_index)(iterator.Second_Cell_Index()));}
                                temperature_vector=laplace_grid.Frame(grid_index).r.Inverse_Rotate(temperature_vector);
                                if(flip) (*face_temperature_pointer(grid_index))(iterator.Full_Index())+=temperature_vector(iterator.Axis());
                                else (*face_temperature_pointer(grid_index))(iterator.Full_Index())=temperature_vector(iterator.Axis());}}
            Write_Substep("updated faces",0,2);}
        else{
            laplace_diffusion.Construct_System(incompressible_callbacks,time,example.fluids_parameters.use_trapezoid_rule?((T).5*dt):dt,true);
            
            ARRAY<T_ARRAYS_SCALAR> temperature(n_global_grids);
            ARRAY<T_ARRAYS_SCALAR*> temperature_pointer(n_global_grids);
            for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
                temperature_pointer(grid_index)=&temperature(grid_index);
                if(laplace_grid.Local_Grid(grid_index))
                    temperature(grid_index)=example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->temperature_container.temperature;}
            ARRAY<ARRAY<T> > temperature_boundary;
            Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_diffusion,incompressible_callbacks,time+dt,temperature_pointer,temperature_boundary,0);
            WRITE_CELL_ARRAY(temperature,STRING_UTILITIES::string_sprintf("old cell values"));
            
            Build_Diffusion_System(laplace_diffusion,incompressible_callbacks,time,dt,temperature_pointer,temperature_boundary,x_diffusion,rhs_diffusion);
            Solve_System(laplace_diffusion.system_matrix,rhs_diffusion,x_diffusion,example.fluids_parameters.implicit_viscosity_tolerance,example.fluids_parameters.implicit_viscosity_iterations,&laplace_diffusion.partition);
            WRITE_CELL_VECTOR(laplace_diffusion,x_diffusion,"x solution");
            
            ARRAY<T_ARRAYS_SCALAR> temperature_updated(n_global_grids);
            ARRAY<T_ARRAYS_SCALAR*> temperature_updated_pointer(n_global_grids);
            for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
                temperature_updated_pointer(grid_index)=&temperature_updated(grid_index);
                if(laplace_grid.Local_Grid(grid_index))
                    temperature_updated(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1));}
            ARRAY<ARRAY<T> > temperature_updated_boundary;
            Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_diffusion,incompressible_callbacks,time+dt,temperature_updated_pointer,temperature_updated_boundary,&x_diffusion);
            WRITE_CELL_ARRAY(temperature_updated,STRING_UTILITIES::string_sprintf("new cell values"));

            for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                if(laplace_grid.Local_Grid(grid_index))
                    for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                        if(!example.use_flip_cell_update || laplace_diffusion.Matrix_Cell_Index(grid_index,iterator.Cell_Index()))
                            example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->temperature_container.temperature(iterator.Cell_Index())=temperature_updated(grid_index)(iterator.Cell_Index());
                        else
                            example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->temperature_container.temperature(iterator.Cell_Index())+=temperature_updated(grid_index)(iterator.Cell_Index())-temperature(grid_index)(iterator.Cell_Index());}
        
            Write_Substep(STRING_UTILITIES::string_sprintf("updated cells"),0,2);}}
    else{
        ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_pointer(n_global_grids);
        ARRAY<T_ARRAYS_SCALAR*> pressures_pointer(n_global_grids);
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
            face_velocities_pointer(grid_index)=laplace_grid.Local_Grid(grid_index)?&example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->face_velocities:0;
            pressures_pointer(grid_index)=laplace_grid.Local_Grid(grid_index)?&example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->pressure:0;}
        laplace_incompressible.Construct_System(incompressible_callbacks,time,dt,false);

        Write_Substep(STRING_UTILITIES::string_sprintf("after building laplace incompressible"),0,2);

        //enforce incompressibility before viscosity solve
        //Pack_And_Compute_Voronoi_Face_Values(laplace_incompressible,face_velocities_pointer,face_velocities_packed,time,dt);
        //Build_Poisson_System(laplace_incompressible,callbacks,time,dt,&face_velocities_packed,x_incompressible,rhs_incompressible);
        //Solve_System(laplace_incompressible.system_matrix,rhs_incompressible,x_incompressible,&laplace_incompressible.partition);
        //Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_incompressible,callbacks,time+dt,pressures_pointer,pressures_boundary,&x_incompressible);
        //Update_Velocities(callbacks,face_velocities_pointer,pressures_pointer,time,dt);
        
        //Pack_And_Compute_Voronoi_Face_Values(laplace_incompressible,face_velocities_pointer,face_velocities_packed,time,dt);
        //Compute_Vorticity(time,dt);
        //Write_Substep(STRING_UTILITIES::string_sprintf("vorticity before viscosity"),0,2);
        
        if(example.fluids_parameters.implicit_viscosity){
            VECTOR<ARRAY<T_ARRAYS_SCALAR*>,TV::dimension> cell_velocity_components_pointer,cell_velocity_components_updated_pointer;
            for(int axis=1;axis<=TV::dimension;axis++){
                cell_velocity_components(axis).Resize(n_global_grids);
                cell_velocity_components_updated(axis).Resize(n_global_grids);
                cell_velocity_components_pointer(axis).Resize(n_global_grids);
                cell_velocity_components_updated_pointer(axis).Resize(n_global_grids);
                for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
                    cell_velocity_components_pointer(axis)(grid_index)=&cell_velocity_components(axis)(grid_index);
                    cell_velocity_components_updated_pointer(axis)(grid_index)=&cell_velocity_components_updated(axis)(grid_index);
                    if(laplace_grid.Local_Grid(grid_index)){
                        cell_velocity_components(axis)(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1));
                        cell_velocity_components_updated(axis)(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1));}}}

            for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                if(laplace_grid.Local_Grid(grid_index)){
                    for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                        TV velocity_vector;
                        for(int axis=1;axis<=TV::dimension;axis++)
                            velocity_vector(axis)=(T).5*(face_velocities_pointer(grid_index)->Component(axis)(iterator.Cell_Index())+face_velocities_pointer(grid_index)->Component(axis)(iterator.Cell_Index()+TV_INT::Axis_Vector(axis)));
                        velocity_vector=laplace_grid.Frame(grid_index).r.Rotate(velocity_vector);
                        for(int axis=1;axis<=TV::dimension;axis++)
                            cell_velocity_components(axis)(grid_index)(iterator.Cell_Index())=velocity_vector(axis);}}
            
            VECTOR_ND<T> rhs_diffusion,x_diffusion;
            
            for(int axis=1;axis<=TV::dimension;axis++){
                viscosity_callbacks.direction=TV::Axis_Vector(axis);
                laplace_diffusion.Construct_System(viscosity_callbacks,time,example.fluids_parameters.use_trapezoid_rule?((T).5*dt):dt,true);
                Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_diffusion,viscosity_callbacks,time,cell_velocity_components_pointer(axis),cell_velocity_boundary_components(axis),0);
                Build_Diffusion_System(laplace_diffusion,viscosity_callbacks,time,dt,cell_velocity_components_pointer(axis),cell_velocity_boundary_components(axis),x_diffusion,rhs_diffusion);
                Solve_System(laplace_diffusion.system_matrix,rhs_diffusion,x_diffusion,example.fluids_parameters.implicit_viscosity_tolerance,example.fluids_parameters.implicit_viscosity_iterations,&laplace_diffusion.partition);
                Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_diffusion,viscosity_callbacks,time+dt,cell_velocity_components_updated_pointer(axis),cell_velocity_boundary_components_updated(axis),&x_diffusion);
                
                WRITE_CELL_ARRAY_POINTER(cell_velocity_components_pointer(axis),STRING_UTILITIES::string_sprintf("old values %d",axis));
                WRITE_CELL_ARRAY_POINTER(cell_velocity_components_updated_pointer(axis),STRING_UTILITIES::string_sprintf("new values %d",axis));}
            
            for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                if(laplace_grid.Local_Grid(grid_index))
                    for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                        T& face_velocity=(*face_velocities_pointer(grid_index))(iterator.Full_Index());
                        if(!incompressible_callbacks.Get_Neumann_Boundary_Condition(laplace_grid.Frame(grid_index)*iterator.Location(),laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(iterator.Axis())),face_velocity,time+dt)){
                            TV velocity_vector;
                            GRID_CELL_INDEX first_cell=GRID_CELL_INDEX(grid_index,iterator.First_Cell_Index()),second_cell=GRID_CELL_INDEX(grid_index,iterator.Second_Cell_Index());
                            if(!laplace_grid.Grid(grid_index).Inside_Domain(first_cell.y))
                                laplace_grid.chimera_grid.Find_Real_Grid_Cell_Index(first_cell);
                            if(!laplace_grid.Grid(grid_index).Inside_Domain(second_cell.y))
                                laplace_grid.chimera_grid.Find_Real_Grid_Cell_Index(second_cell);
                            bool flip=laplace_grid.Chimera_Cell(first_cell.x,first_cell.y) || laplace_grid.Chimera_Cell(second_cell.x,second_cell.y);
                            for(int axis=1;axis<=TV::dimension;axis++){
                                velocity_vector(axis)=(T).5*(cell_velocity_components_updated(axis)(grid_index)(iterator.First_Cell_Index())+cell_velocity_components_updated(axis)(grid_index)(iterator.Second_Cell_Index()));
                                if(flip) velocity_vector(axis)-=(T).5*(cell_velocity_components(axis)(grid_index)(iterator.First_Cell_Index())+cell_velocity_components(axis)(grid_index)(iterator.Second_Cell_Index()));}
                            velocity_vector=laplace_grid.Frame(grid_index).r.Inverse_Rotate(velocity_vector);
                            if(flip) face_velocity+=velocity_vector(iterator.Axis());
                            else face_velocity=velocity_vector(iterator.Axis());
                            
                            if(Is_NaN((*face_velocities_pointer(grid_index))(iterator.Full_Index()))){ 
                                LOG::cout<<"NaN face velocity at grid index "<<grid_index<<", axis "<<iterator.Axis()<<", index "<<iterator.Face_Index()<<std::endl;
                                PHYSBAM_FATAL_ERROR("Face velocities should not be NaN here!");}}}
            Write_Substep(STRING_UTILITIES::string_sprintf("updated viscosity solve %f",example.fluids_parameters.viscosity),0,2);}
        else{
            for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
                if(laplace_grid.Local_Grid(grid_index))
                    for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                        T& face_velocity=(*face_velocities_pointer(grid_index))(iterator.Full_Index());
                        incompressible_callbacks.Get_Neumann_Boundary_Condition(laplace_grid.Frame(grid_index)*iterator.Location(),laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(iterator.Axis())),face_velocity,time+dt);}}

        //Pack_And_Compute_Voronoi_Face_Values(laplace_incompressible,face_velocities_pointer,face_velocities_packed,time,dt);
        //Compute_Vorticity(time,dt);
        //Write_Substep(STRING_UTILITIES::string_sprintf("vorticity after viscosity"),0,2);
        
        // enforce incompressibility after viscosity solve
        Pack_And_Compute_Voronoi_Face_Values(laplace_incompressible,face_velocities_pointer,face_velocities_packed,time,dt);
        voronoi_face_velocities.Resize(laplace_grid.voronoi_faces.Size());
        for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
            int matrix_face_index=laplace_incompressible.Matrix_Face_Index(voronoi_face_index);
            if(matrix_face_index)
                voronoi_face_velocities(voronoi_face_index)=face_velocities_packed(matrix_face_index);}
        laplace_grid.voronoi.voronoi_face_velocities=voronoi_face_velocities;
        //Compute_Vorticity(time,dt);
        
        Write_Substep(STRING_UTILITIES::string_sprintf("voronoi face velocities"),0,2);
        
        Build_Poisson_System(laplace_incompressible,incompressible_callbacks,time,dt,&face_velocities_packed,x_incompressible,rhs_incompressible);
        Solve_System(laplace_incompressible.system_matrix,rhs_incompressible,x_incompressible,example.fluids_parameters.incompressible_tolerance,example.fluids_parameters.incompressible_iterations,&laplace_incompressible.partition);
        //Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_incompressible,incompressible_callbacks,time+dt,pressures_pointer,pressures_boundary,&x_incompressible);
        
        Unpack_Cell_Values(laplace_incompressible,incompressible_callbacks,time+dt,pressures_pointer,&x_incompressible);
        Pack_And_Exchange_Boundary_Cell_Values(pressures_pointer,pressures_boundary);
        ARRAY<T_ARRAYS_SCALAR*> pressures_pointer_local(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++) pressures_pointer_local(grid_index)=pressures_pointer(laplace_grid.Global_Grid_Index(grid_index));
        example.Exchange_Split_Grid_Ghost_Cells(pressures_pointer_local);

        Update_Velocities(incompressible_callbacks,face_velocities_pointer,pressures_pointer,time,dt);}

    laplace_grid.Clean_Memory();
    laplace_incompressible.Clean_Memory();
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    driver->Write_Substep(title,substep,level);
}
//#####################################################################
// Function Pack_And_Compute_Voronoi_Velocities
//#####################################################################
template<class TV_INT> int Neighbor_Axis(const TV_INT& index_1,const TV_INT& index_2)
{
    TV_INT offset=index_2-index_1;
    if(offset.L1_Norm()==1)
        return offset.Dominant_Axis();
    return 0;
}
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Pack_And_Compute_Voronoi_Face_Values(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_values,VECTOR_ND<T>& packed_values,const T time,const T dt)
{
    LOG::SCOPE scope("Pack_And_Compute_Voronoi_Velocities");
    
    packed_values.Resize(laplace.n_matrix_faces);
    //pack regular face values
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        if(laplace_grid.Local_Grid(grid_index))
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                D_FACE_INDEX face_index=iterator.Full_Index();
                int matrix_face_index=laplace.Matrix_Face_Index(grid_index,face_index);
                if(matrix_face_index)
                    packed_values(matrix_face_index)=(*face_values(grid_index))(face_index);}}

    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        int matrix_face_index=laplace.Matrix_Face_Index(voronoi_face_index);
        if(matrix_face_index){
            VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=laplace_grid.voronoi_faces(voronoi_face_index).x;
            int grid_index_1=example.chimera_grid->Parent_Grid(grid_cell_indices(1).x);
            int grid_index_2=example.chimera_grid->Parent_Grid(grid_cell_indices(2).x);
            TV_INT cell_index_1=example.chimera_grid->Unsplit_Cell_Index(grid_cell_indices(1));
            TV_INT cell_index_2=example.chimera_grid->Unsplit_Cell_Index(grid_cell_indices(2));
            int neighbor_axis=Neighbor_Axis(cell_index_1,cell_index_2);
            
            if(grid_index_1==grid_index_2 && neighbor_axis){
                TV_INT index=TV_INT::Componentwise_Max(cell_index_1,cell_index_2);
                int grid_index_local=laplace_grid.Local_Grid(grid_cell_indices(1).x)?grid_cell_indices(1).x:grid_cell_indices(2).x;
                TV_INT index_local=example.chimera_grid->Split_Cell_Index(grid_index_local,index);
                packed_values(matrix_face_index)=face_values(grid_index_local)->Component(neighbor_axis)(index_local)*(index==cell_index_2?1:-1);}
            else
                packed_values(matrix_face_index)=std::numeric_limits<T>::quiet_NaN();}}

    ARRAY<ARRAY<int> > voronoi_face_index_interpolated_from_others(n_global_grids);
    ARRAY<ARRAY<PAIR<TV,TV> > > voronoi_face_positions_interpolated_from_others(n_global_grids);
    ARRAY<ARRAY<T> > voronoi_face_values_interpolated_from_others(n_global_grids);
    ARRAY<ARRAY<PAIR<TV,TV> > > voronoi_face_positions_interpolated_here(n_global_grids);
    ARRAY<ARRAY<T> > voronoi_face_values_interpolated_here(n_global_grids);

    //interpolate voronoi face values from local grids and pack requests for nonlocal grids
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        int matrix_face_index=laplace.Matrix_Face_Index(voronoi_face_index);
        if(matrix_face_index && Is_NaN(packed_values(matrix_face_index))){
            VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=laplace_grid.voronoi_faces(voronoi_face_index).x;
            TV location_1=laplace_grid.Frame(grid_cell_indices(1).x)*laplace_grid.Grid(grid_cell_indices(1).x).X(grid_cell_indices(1).y);
            TV location_2=laplace_grid.Frame(grid_cell_indices(2).x)*laplace_grid.Grid(grid_cell_indices(2).x).X(grid_cell_indices(2).y);
            //TV center=laplace_grid.voronoi_faces_geometry(voronoi_face_index).Center(); //centroid
            TV center=(T).5*(location_1+location_2);
            TV normal=(location_2-location_1).Normalized();
            int finest_grid=example.chimera_grid->Find_Finest_Grid(center,false);
            if(!example.chimera_grid->use_mpi || laplace_grid.Local_Grid(finest_grid)){
                TV interpolated_vector=laplace_grid.Frame(finest_grid).r.Rotate(interpolation.Clamped_To_Array_Face(laplace_grid.Grid(finest_grid),*face_values(finest_grid),laplace_grid.Frame(finest_grid).Inverse_Times(center)));
                //LOG::cout << "interpolated face value " << grid_cell_indices(1).x << " " << grid_cell_indices(1).y << " " << grid_cell_indices(2).x << " " << grid_cell_indices(2).y << " " << interpolated_vector << std::endl;
                packed_values(matrix_face_index)=TV::Dot_Product(normal,interpolated_vector);}
            else{
                voronoi_face_index_interpolated_from_others(finest_grid).Append(voronoi_face_index);
                voronoi_face_positions_interpolated_from_others(finest_grid).Append(PAIR<TV,TV>(center,normal));}}}
    
    //exchange requests, interpolate requested locations locally, exchange response values, pack
    if(example.chimera_grid->use_mpi && n_local_grids){
        example.chimera_grid->Exchange_Voronoi_Face_Positions_For_Interpolation(voronoi_face_positions_interpolated_from_others,voronoi_face_positions_interpolated_here);
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
            voronoi_face_values_interpolated_here(grid_index).Resize(voronoi_face_positions_interpolated_here(grid_index).Size());
            for(int face_index=1;face_index<=voronoi_face_positions_interpolated_here(grid_index).Size();face_index++){
                TV& center=voronoi_face_positions_interpolated_here(grid_index)(face_index).x;
                TV& normal=voronoi_face_positions_interpolated_here(grid_index)(face_index).y;
                int finest_grid=example.chimera_grid->Find_Finest_Grid(center,false);
                TV interpolated_vector=laplace_grid.Frame(finest_grid).r.Rotate(interpolation.Clamped_To_Array_Face(laplace_grid.Grid(finest_grid),*face_values(finest_grid),laplace_grid.Frame(finest_grid).Inverse()*center));
                voronoi_face_values_interpolated_here(grid_index)(face_index)=TV::Dot_Product(normal,interpolated_vector);}}
        example.chimera_grid->Exchange_Voronoi_Face_Interpolated_Values(voronoi_face_values_interpolated_here,voronoi_face_values_interpolated_from_others);
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
            for(int face_index=1;face_index<=voronoi_face_index_interpolated_from_others(grid_index).Size();face_index++){
                int voronoi_face_index=voronoi_face_index_interpolated_from_others(grid_index)(face_index);
                int matrix_face_index=laplace.Matrix_Face_Index(voronoi_face_index);
                packed_values(matrix_face_index)=voronoi_face_values_interpolated_from_others(grid_index)(face_index);}}}
}
//#####################################################################
// Function Build_Diffusion_System
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::Build_Diffusion_System(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const T dt,ARRAY<T_ARRAYS_SCALAR*>& u,ARRAY<ARRAY<T> >& u_boundary,VECTOR_ND<T>& x,VECTOR_ND<T>& rhs,const bool compute_explicit_forces_only)
{
    //dirichlet should be set in objects, however we allow slip for now with neumann
    //these should be set to 0 at the domain boundaries to prevent so that velocities are mirrored in the ghost regions
    x.Resize(laplace.n_matrix_cells);
    rhs.Resize(laplace.n_matrix_cells);
    rhs.Fill(0);

    //build right hand side
    VECTOR_ND<T> temperature_time_n;
    if(example.fluids_parameters.use_trapezoid_rule)
        temperature_time_n.Resize(laplace.n_matrix_cells);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    if(!compute_explicit_forces_only)
                        rhs(matrix_cell_index)=laplace_grid.Cell_Size(grid_index,cell_index)*(*u(grid_index))(cell_index);
                    if(example.fluids_parameters.use_trapezoid_rule)
                        temperature_time_n(matrix_cell_index)=(*u(grid_index))(cell_index);}}
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Boundary_Grid(grid_index))
            for(int boundary_cell_index=1;boundary_cell_index<=laplace_grid.boundary_cell_indices(grid_index).Size();boundary_cell_index++){
                TV_INT cell_index=laplace_grid.boundary_cell_indices(grid_index)(boundary_cell_index);
                int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    if(example.fluids_parameters.use_trapezoid_rule)
                        temperature_time_n(matrix_cell_index)=u_boundary(grid_index)(boundary_cell_index);}}
    
    VECTOR_ND<T> dirichlet_divergence(laplace.n_matrix_cells);
    VECTOR_ND<T> neumann_divergence(laplace.n_matrix_cells);
    laplace.Compute_Dirichlet_Divergence(callbacks,time+dt,dirichlet_divergence);
    laplace.Compute_Neumann_Divergence(callbacks,time+dt,neumann_divergence);
    if(example.fluids_parameters.use_trapezoid_rule){
        VECTOR_ND<T> dirichlet_divergence_time_n(laplace.n_matrix_cells);
        VECTOR_ND<T> neumann_divergence_time_n(laplace.n_matrix_cells);
        laplace.Compute_Dirichlet_Divergence(callbacks,time,dirichlet_divergence_time_n);
        laplace.Compute_Neumann_Divergence(callbacks,time,neumann_divergence_time_n);
        
        VECTOR_ND<T> gradient_time_n(laplace.n_matrix_faces);
        VECTOR_ND<T> explicit_time_n(laplace.n_matrix_cells);
        laplace.Multiply_By_Gradient(temperature_time_n,gradient_time_n);
        laplace.Multiply_By_Inverse_Mass(gradient_time_n,gradient_time_n);
        laplace.Multiply_By_Volume_Weighted_Divergence(gradient_time_n,explicit_time_n);
        
        rhs+=dt*(T).5*(dirichlet_divergence+neumann_divergence+dirichlet_divergence_time_n+neumann_divergence_time_n+explicit_time_n);}
    else
        rhs+=dt*(dirichlet_divergence+neumann_divergence);
    
    WRITE_CELL_VECTOR(laplace,rhs,"rhs diffusion");
}
//#####################################################################
// Function Build_Poisson_System
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Build_Poisson_System(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const T dt,VECTOR_ND<T>* packed_values,VECTOR_ND<T>& x,VECTOR_ND<T>& rhs)
{
    LOG::SCOPE scope("Build_System");

    x.Resize(laplace.n_matrix_cells);
    rhs.Resize(laplace.n_matrix_cells);
    rhs.Fill(0);

    //build right hand side
    if(example.solve_poisson_equation){
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
            if(laplace_grid.Local_Grid(grid_index))
                for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                    TV_INT cell_index=iterator.Cell_Index();
                    int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                    if(matrix_cell_index){
                        GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                        TV location=laplace_grid.Frame(grid_index)*iterator.Location();
                        rhs(matrix_cell_index)=-laplace_grid.Cell_Size(grid_index,cell_index)*example.Get_Analytic_Laplacian(location);}}}
    else{
        //VECTOR_ND<T> rhs_velocity(n_matrix_cells);
        //divergence_matrix.Times(face_velocities_packed,rhs_velocity);
        
        /*TV_INT ci;
        ci(1)=10;
        ci(2)=13;
        const ARRAY<int>& incident=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(GRID_CELL_INDEX(1,ci));
        for(int i=1;i<=incident.Size();i++)
            LOG::cout << "incident face value " << packed_values->operator()(laplace.Matrix_Face_Index(incident(i))) << std::endl;*/

        laplace.Multiply_By_Volume_Weighted_Divergence(*packed_values,rhs);
    
        /*for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
            if(laplace_grid.Local_Grid(grid_index)){
                TV dx=laplace_grid.Grid(grid_index).DX();
                int local_grid_index=example.chimera_grid->global_to_local_grid_index_map(grid_index);
                for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                    TV_INT cell_index=iterator.Cell_Index();
                    int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                    if(matrix_cell_index){
                        T divergence=0;
                        for(int axis=1;axis<=TV::dimension;axis++)
                            for(int side=1;side<=2;side++)
                                divergence+=(2*side-3)*FACE_VELOCITIES(local_grid_index).Component(axis)(cell_index+(side-1)*TV_INT::Axis_Vector(axis))/dx(axis);
                        //if(grid_index==2 && cell_index(1)==1 && cell_index(2)==24)
                            //LOG::cout << "divergence " << divergence << " " << laplace_grid.Cell_Size(grid_index,cell_index) << " " << rhs(matrix_cell_index) << std::endl;
                        rhs(matrix_cell_index)=laplace_grid.Cell_Size(grid_index,cell_index)*divergence;////laplace_grid.Grid(grid_index).Cell_Size()*divergence;
                        }}}*/
        
        VECTOR_ND<T> divergence(rhs.Size());
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
            if(laplace_grid.Local_Grid(grid_index)){
                for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                    TV_INT cell_index=iterator.Cell_Index();
                    int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                    if(matrix_cell_index)
                        divergence(matrix_cell_index)=rhs(matrix_cell_index)/laplace_grid.Cell_Size(grid_index,cell_index);}}
        WRITE_CELL_VECTOR(laplace_incompressible,divergence,"divergence");

        rhs*=-(T)1.0/dt;}

    //add boundary conditions
    VECTOR_ND<T> dirichlet_divergence(laplace.n_matrix_cells);
    VECTOR_ND<T> neumann_divergence(laplace.n_matrix_cells);
    laplace.Compute_Dirichlet_Divergence(callbacks,time,dirichlet_divergence);
    laplace.Compute_Neumann_Divergence(callbacks,time,neumann_divergence);
    
    //positive since the system matrix is the negative laplacian
    rhs+=dirichlet_divergence+(example.solve_poisson_equation?(T)1:-(T)1.0/dt)*neumann_divergence;
    
    WRITE_CELL_VECTOR(laplace_incompressible,rhs,"rhs incompressible");
}
//#####################################################################
// Function Solve_System
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Solve_System(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& rhs,VECTOR_ND<T>& x,T tolerance,int iterations,SPARSE_MATRIX_PARTITION* partition)
{
    LOG::SCOPE scope("Solve_System");

    //ZERO SOLUTION VECTOR
    //x_pcg.Resize(n_matrix_cells+n_kinematic_faces);
    int n_variables=A.n;
    x.Resize(n_variables);
    x.Fill(0);
    
    //COMPUTE TOLERANCE AS PROPORTIONAL TO THE DIVERGENCE OF THE ENTIRE DOMAIN ??? WORK THIS OUT

    //SOLVE USING MPI OR SEQUENTIAL PCG
    if(!example.chimera_grid->use_mpi){
        VECTOR_ND<T> q(n_variables),s(n_variables),r(n_variables),k(n_variables),z(n_variables);
        PCG_SPARSE<T> pcg;
        if(example.solve_heat_equation_on_faces_coupled) pcg.Use_Conjugate_Gradient();
        else pcg.Use_Incomplete_Cholesky();
        pcg.maximum_iterations=iterations;
        pcg.cg_restart_iterations=100;
        pcg.show_results=true;
        pcg.show_residual=example.fluids_parameters.show_residual;
        pcg.Solve(A,x,rhs,q,s,r,k,z,tolerance);}
    else{
        assert(partition);
        example.chimera_grid->Parallel_Solve(A,x,rhs,tolerance,*partition,laplace_grid.is_boundary_grid,iterations,true);}
}
//#####################################################################
// Function Unpack_Values
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Unpack_Cell_Values(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,ARRAY<T_ARRAYS_SCALAR*>& cell_values,VECTOR_ND<T>* packed_values)
{
    //UNPACK SOLVED AND BOUNDARY CONDITION PRESSURES, ZERO OUT OTHERS
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(laplace_grid.Chimera_Cell(grid_index,cell_index)){
                    int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                    if(matrix_cell_index){
                        if(packed_values) (*cell_values(grid_index))(cell_index)=(*packed_values)(matrix_cell_index);}}
                else
                    (*cell_values(grid_index))(cell_index)=std::numeric_limits<T>::quiet_NaN();}
    //LOG::cout << "unpacked values" << std::endl;

    Set_Dirichlet_Boundary_Conditions(callbacks,time,cell_values);
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Set_Dirichlet_Boundary_Conditions(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,ARRAY<T_ARRAYS_SCALAR*>& cell_values)
{
    LOG::SCOPE scope("Set_Dirichlet_Boundary_Conditions");
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next()){
                TV location=laplace_grid.Frame(grid_index)*iterator.Location();
                T value;
                if(callbacks.Get_Dirichlet_Boundary_Condition(location,value,time))
                    (*cell_values(grid_index))(iterator.Cell_Index())=value;}
}
//#####################################################################
// Function Pack_And_Exchange_Boundary_Cell_Values
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Pack_And_Exchange_Boundary_Cell_Values(ARRAY<T_ARRAYS_SCALAR*>& cell_values,ARRAY<ARRAY<T> >& boundary_cell_values)
{
    LOG::SCOPE scope("Pack_And_Exchange_Boundary_Cell_Values");
    //EXCHANGE ALL BOUNDARY PRESSURES WITH NEIGHBORING GRIDS AND SET IN PRESSURES
    boundary_cell_values.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        boundary_cell_values(grid_index).Resize(laplace_grid.boundary_cell_indices(grid_index).Size());
        if(laplace_grid.Local_Grid(grid_index))
            for(int i=1;i<=laplace_grid.boundary_cell_indices(grid_index).Size();i++)
                boundary_cell_values(grid_index)(i)=(*cell_values(grid_index))(laplace_grid.boundary_cell_indices(grid_index)(i));}
    example.chimera_grid->Exchange_Boundary_Scalar_Values(boundary_cell_values,laplace_grid.is_boundary_grid);
}
//#####################################################################
// Function Interpolate_Cell_Values_From_Delaunay_Simplices
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Interpolate_Cell_Values_From_Delaunay_Simplices(ARRAY<T_ARRAYS_SCALAR*>& cell_values,ARRAY<ARRAY<T> >& boundary_cell_values)
{
    LOG::SCOPE scope("Interpolate_Cell_Values_From_Delaunay_Simplices");
    CELL_LOOKUP_CHIMERA<T_GRID> cell_lookup(laplace_grid,cell_values,boundary_cell_values);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(Is_NaN((*cell_values(grid_index))(cell_index))){
                    TV location=laplace_grid.Frame(grid_index)*iterator.Location();
                    if(!laplace_grid.Chimera_Cell(grid_index,cell_index))
                        (*cell_values(grid_index))(cell_index)=laplace_grid.interpolation.From_Cell_Centers(cell_lookup,location);}}
}
//#####################################################################
// Function Unpack_And_Interpolate_And_Inject_Pressures
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Unpack_And_Interpolate_And_Inject_Cell_Values(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,ARRAY<T_ARRAYS_SCALAR*>& cell_values,ARRAY<ARRAY<T> >& boundary_cell_values,VECTOR_ND<T>* packed_values)
{
    LOG::SCOPE scope("Unpack_And_Interpolate_And_Inject_Cell_Values");

    Unpack_Cell_Values(laplace,callbacks,time,cell_values,packed_values);
    Pack_And_Exchange_Boundary_Cell_Values(cell_values,boundary_cell_values);
    Interpolate_Cell_Values_From_Delaunay_Simplices(cell_values,boundary_cell_values);
    Write_Substep("after interpolating from delaunay simplices",0,5);

    //this extra copy is needed because cell values is of size n_global_grids while the coupling and ghost cell functions take arrays containing only local grid arrays
    ARRAY<T_ARRAYS_SCALAR*> cell_values_tmp(n_local_grids);
    for(int grid_index=1;grid_index<=n_local_grids;grid_index++) cell_values_tmp(grid_index)=cell_values(laplace_grid.Global_Grid_Index(grid_index));
    example.Coupling_Overlap_Regions_Cell(1,cell_values_tmp,true);

    Write_Substep("after coupling overlapped regions and filling ghost cells",0,5);

    Set_Dirichlet_Boundary_Conditions(callbacks,time,cell_values);
    LOG::cout << "set dirichlet boundary conditions at the end of Unpack_And_Interpolate_And_Inject_Cell_Values" << std::endl;

    Write_Substep("after cell value interpolation",0,5);
}
//#####################################################################
// Function Update_Velocities
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::Add_Cell_Faces_Least_Squares_Terms(ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,const GRID_CELL_INDEX& grid_cell_index,const TV& x,D2_MATRIX& A,TV2& rhs,HASHTABLE<int>& added_voronoi_faces)
{
    if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
        ARRAY<int> incident_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
        for(int i=1;i<=incident_indices.Size();i++){
            int voronoi_face_index=incident_indices(i);
            if(!added_voronoi_faces.Contains(voronoi_face_index)){
                added_voronoi_faces.Set(voronoi_face_index);
                
                TV location_1,location_2,location,normal;T distance;
                laplace_grid.Voronoi_Face(voronoi_face_index,location_1,location_2,location,normal,distance);
                TV offset=location-x;
                TV2 j;
                for(int axis=1;axis<=TV::dimension;axis++){
                    j(axis)=normal(axis);
                    for(int axis2=1;axis2<=TV::dimension;axis2++)
                        j(TV::dimension+(axis-1)*TV::dimension+axis2)=normal(axis2)*offset(axis);}

                A+=D2_MATRIX::Outer_Product(j,j);
                rhs+=voronoi_face_velocities(voronoi_face_index)*j;}}}

    T_GRID& grid=laplace_grid.Grid(grid_cell_index.x);
    const FRAME<TV>& frame=laplace_grid.Frame(grid_cell_index.x);
    for(int axis=1;axis<=TV::dimension;axis++)
        for(int side=1;side<=2;side++){
            D_FACE_INDEX face_index(axis,grid_cell_index.y+(side-1)*TV_INT::Axis_Vector(axis));
            if(laplace_grid.Chimera_Face(grid_cell_index.x,face_index)){
                TV location=frame*grid.Face(axis,face_index.index);
                TV normal=frame.r.Rotate(TV::Axis_Vector(axis));
                TV offset=location-x;
                TV2 j;
                for(int axis=1;axis<=TV::dimension;axis++){
                    j(axis)=normal(axis);
                    for(int axis2=1;axis2<=TV::dimension;axis2++)
                        j(TV::dimension+(axis-1)*TV::dimension+axis2)=normal(axis2)*offset(axis);}

                A+=D2_MATRIX::Outer_Product(j,j);
                rhs+=face_velocities(grid_cell_index.x)->operator()(face_index)*j;}}
}
template<class T_GRID> typename T_GRID::VECTOR_T SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::Compute_Cell_Center_Velocity(ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,const int grid_index,const TV_INT& cell_index)
{
    //LOG::cout << "interpolating to cell " << grid_index << " " << cell_index << std::endl;

    GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
    
    TV2 rhs;
    D2_MATRIX A;
    HASHTABLE<int> added_voronoi_faces;
    TV cell_location=laplace_grid.Frame(grid_index)*laplace_grid.Grid(grid_index).X(cell_index);
    
    if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
        ARRAY<int> incident_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
        for(int i=1;i<=incident_indices.Size();i++){
            int voronoi_face_index=incident_indices(i);
            const VECTOR<GRID_CELL_INDEX,2>& cell_indices=laplace_grid.voronoi_faces(voronoi_face_index).x;
            Add_Cell_Faces_Least_Squares_Terms(face_velocities,(cell_indices(1)==grid_cell_index)?cell_indices(2):cell_indices(1),cell_location,A,rhs,added_voronoi_faces);}}
    
    for(int axis=1;axis<=TV::dimension;axis++)
        for(int side=1;side<=2;side++){
            D_FACE_INDEX face_index(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis));
            if(laplace_grid.Chimera_Face(grid_index,face_index))
                Add_Cell_Faces_Least_Squares_Terms(face_velocities,GRID_CELL_INDEX(grid_index,cell_index+(2*side-3)*TV_INT::Axis_Vector(axis)),cell_location,A,rhs,added_voronoi_faces);}
    
    //T dx=laplace_grid.Grid(grid_index).Minimum_Edge_Length();
    T diagonal_term=0;//1*dx*dx;
    for(int axis=1;axis<=TV::dimension;axis++)
        for(int axis2=1;axis<=TV::dimension;axis++){
            int i=TV::dimension+(axis-1)*TV::dimension+axis2;
            A(i,i)+=diagonal_term;}
    
    TV2 x=A.Solve_Linear_System(rhs);

    TV velocity;
    x.Get_Subvector(1,velocity);
    return velocity;
}

template<class T_GRID> typename T_GRID::VECTOR_T SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::Compute_Cell_Center_Velocity_Local(ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,const int grid_index,const TV_INT& cell_index)
{
    //LOG::cout << "interpolating to cell " << grid_index << " " << cell_index << std::endl;

    GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
    
    TV2 rhs;
    D2_MATRIX A;
    T_GRID& grid=laplace_grid.Grid(grid_cell_index.x);
    const FRAME<TV>& frame=laplace_grid.Frame(grid_cell_index.x);
    TV cell_location=frame*grid.X(cell_index);
    
    int n_faces=0;if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
        ARRAY<int> incident_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
        for(int i=1;i<=incident_indices.Size();i++){
            int voronoi_face_index=incident_indices(i);
            
            TV location_1,location_2,location,normal;T distance;
            laplace_grid.Voronoi_Face(voronoi_face_index,location_1,location_2,location,normal,distance);
            TV offset=location-cell_location;
            TV2 j;
            for(int axis=1;axis<=TV::dimension;axis++){
                j(axis)=normal(axis);
                for(int axis2=1;axis2<=TV::dimension;axis2++)
                    j(TV::dimension+(axis-1)*TV::dimension+axis2)=normal(axis2)*offset(axis);}
            
            //LOG::cout << "face " << voronoi_face_index << " " << j << " " << normal << " " << offset << " " << cell_location << " " << location << std::endl;
            
            n_faces++;
            A+=D2_MATRIX::Outer_Product(j,j);
            rhs+=voronoi_face_velocities(voronoi_face_index)*j;}}
    
    for(int axis=1;axis<=TV::dimension;axis++)
        for(int side=1;side<=2;side++){
            D_FACE_INDEX face_index(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis));
            if(laplace_grid.Chimera_Face(grid_index,face_index)){
                TV location=frame*grid.Face(axis,face_index.index);
                TV normal=frame.r.Rotate(TV::Axis_Vector(axis));
                TV offset=location-cell_location;
                TV2 j;
                for(int axis=1;axis<=TV::dimension;axis++){
                    j(axis)=normal(axis);
                    for(int axis2=1;axis2<=TV::dimension;axis2++)
                        j(TV::dimension+(axis-1)*TV::dimension+axis2)=normal(axis2)*offset(axis);}

                //LOG::cout << "face " << face_index << " " << j << " " << normal << " " << offset << std::endl;

                n_faces++;
                A+=D2_MATRIX::Outer_Product(j,j);
                rhs+=face_velocities(grid_cell_index.x)->operator()(face_index)*j;}}
    
    //LOG::cout << "A" << std::endl << A << std::endl;

    //if(n_faces<A.m){
    //T dx=laplace_grid.Grid(grid_index).Minimum_Edge_Length();
        T diagonal_term=(T)1e-5;//dx*dx;
        for(int axis=1;axis<=TV::dimension;axis++)
            for(int axis2=1;axis2<=TV::dimension;axis2++){
                int i=TV::dimension+(axis-1)*TV::dimension+axis2;
                A(i,i)+=diagonal_term;}
//}
    
        //LOG::cout << "A_mod" << std::endl << A << std::endl;
        //LOG::cout << "rhs " << std::endl << rhs << std::endl;

    TV2 x=A.Solve_Linear_System(rhs);

    TV velocity;
    x.Get_Subvector(1,velocity);
    return velocity;
}

template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Update_Velocities(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,ARRAY<T_ARRAYS_SCALAR*>& pressures,const T time,const T dt)
{
    LOG::SCOPE scope("Update_Velocities");
    ARRAY<T_FACE_ARRAYS_SCALAR> face_velocities_update(n_global_grids);
    //ARRAY<T_FACE_ARRAYS_SCALAR*> face_velocities_update_pointer(n_global_grids);
    ARRAY<T> voronoi_face_velocities_update;
    //ARRAY<VECTOR<HASHTABLE<TV_INT,bool>,3> > invalid_faces_indices(n_global_grids);
    
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        //face_velocities_update_pointer(grid_index)=&face_velocities_update(grid_index);
        if(laplace_grid.Local_Grid(grid_index)){
            face_velocities_update(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(),true,false,(T)0);
            TV dx=laplace_grid.Grid(grid_index).DX();
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                D_FACE_INDEX face_index=iterator.Full_Index();
                TV location=laplace_grid.Frame(grid_index)*iterator.Location();TV normal=laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(iterator.Axis()));
                T& face_velocity=(*face_velocities(grid_index))(face_index);
                //T& face_velocity_update=face_velocities_update(grid_index)(face_index);
                if(!callbacks.Get_Neumann_Boundary_Condition(location,normal,face_velocity,time)){
                    TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
                    //face_velocity_update=dt/(example.fluids_parameters.density*dx(iterator.Axis()))*((*pressures(grid_index))(first_cell_index)-(*pressures(grid_index))(second_cell_index));
                    face_velocity+=dt/(example.fluids_parameters.density*dx(iterator.Axis()))*((*pressures(grid_index))(first_cell_index)-(*pressures(grid_index))(second_cell_index));}}}}

    //voronoi_face_velocities.Resize(laplace_grid.voronoi_faces.Size());
    voronoi_face_velocities_update.Resize(laplace_grid.voronoi_faces.Size());
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=laplace_grid.voronoi_faces(voronoi_face_index).x;
        if(!laplace_grid.Local_Grid(grid_cell_indices(1).x) && !laplace_grid.Local_Grid(grid_cell_indices(2).x)) continue;
        TV location_1,location_2,location,normal;T distance;
        laplace_grid.Voronoi_Face(voronoi_face_index,location_1,location_2,location,normal,distance);
        T& face_velocity=voronoi_face_velocities(voronoi_face_index);
        //T& face_velocity_update=voronoi_face_velocities_update(voronoi_face_index);
        if(!callbacks.Get_Neumann_Boundary_Condition(location,normal,face_velocity,time)){
            VECTOR<T,2> pressures;
            for(int i=1;i<=2;i++){
                if(laplace_grid.Local_Grid(grid_cell_indices(i).x)) pressures(i)=example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_cell_indices(i).x))->pressure(grid_cell_indices(i).y);
                else pressures(i)=pressures_boundary(grid_cell_indices(i).x)(laplace_grid.boundary_cell_indices_to_linear_index(grid_cell_indices(i).x).Get(grid_cell_indices(i).y));}
            //face_velocity_update=dt/(example.fluids_parameters.density*distance)*(pressures(1)-pressures(2));
            face_velocity+=dt/(example.fluids_parameters.density*distance)*(pressures(1)-pressures(2));}}
    laplace_grid.voronoi.voronoi_face_velocities=voronoi_face_velocities;

    Write_Substep("after applying pressure differences",0,1);

    ///////////////////////////////////////////////////////////////////////////////debug
    VECTOR_ND<T> divergence(rhs_incompressible.Size());
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            T_GRID& grid=laplace_grid.Grid(grid_index);
            FRAME<TV> frame=laplace_grid.Frame(grid_index);
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=laplace_incompressible.Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                    T d=0;
                    if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                        ARRAY<int> incident_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                        for(int i=1;i<=incident_indices.Size();i++){
                            int voronoi_face_index=incident_indices(i);
                            VECTOR<GRID_CELL_INDEX,2> gcis=laplace_grid.voronoi_faces(voronoi_face_index).x;
                            d+=(gcis(1)==grid_cell_index?-1:1)*laplace_grid.voronoi_faces(voronoi_face_index).y*voronoi_face_velocities_update(voronoi_face_index);}}
                    for(int axis=1;axis<=TV::dimension;axis++)
                        for(int side=1;side<=2;side++){
                            D_FACE_INDEX face_index(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis));
                            if(laplace_grid.Chimera_Face(grid_index,face_index)){
                                d+=(3-2*side)*grid.Face_Size(axis)*face_velocities_update(grid_index)(face_index);}}
                    divergence(matrix_cell_index)=d/laplace_grid.Cell_Size(grid_index,cell_index);}}}
    WRITE_CELL_VECTOR(laplace_incompressible,divergence,"divergence free cells");
    ///////////////////////////////////////////////////////////////////////////////debug

    /*VORONOI_INTERPOLATION_CHIMERA<T_GRID> interpolation(laplace_grid,face_velocities,voronoi_face_velocities);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            TV dx=laplace_grid.Grid(grid_index).DX();
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                D_FACE_INDEX face_index=iterator.Full_Index();
                T& value=(*face_velocities(grid_index))(face_index);
                if(Is_NaN(value)){
                    TV normal=laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(face_index.axis));
                    TV location=laplace_grid.Frame(grid_index)*iterator.Location();
                    value=TV::Dot_Product(normal,interpolation.From_Voronoi_Faces(location));}}}*/

    VECTOR<ARRAY<T_ARRAYS_SCALAR*>,TV::dimension> cell_velocity_components_pointer;
    for(int axis=1;axis<=TV::dimension;axis++){
        cell_velocity_components(axis).Resize(n_global_grids);
        cell_velocity_components_pointer(axis).Resize(n_global_grids);
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
            cell_velocity_components_pointer(axis)(grid_index)=&cell_velocity_components(axis)(grid_index);
            if(laplace_grid.Local_Grid(grid_index))
                cell_velocity_components(axis)(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1));}}

    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            //T_GRID& grid=laplace_grid.Grid(grid_index);
            FRAME<TV> frame=laplace_grid.Frame(grid_index);
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(laplace_grid.Chimera_Cell(grid_index,cell_index)){
                    TV velocity_vector;
                    if(laplace_grid.Boundary_Cell(grid_index,cell_index)){
                        velocity_vector=Compute_Cell_Center_Velocity_Local(face_velocities,grid_index,cell_index);
                    }else{
                        for(int axis=1;axis<=TV::dimension;axis++)
                            velocity_vector(axis)=(T).5*(face_velocities(grid_index)->Component(axis)(iterator.Cell_Index())+face_velocities(grid_index)->Component(axis)(iterator.Cell_Index()+TV_INT::Axis_Vector(axis)));
                        velocity_vector=laplace_grid.Frame(grid_index).r.Rotate(velocity_vector);}

                    for(int axis=1;axis<=TV::dimension;axis++)
                        cell_velocity_components(axis)(grid_index)(iterator.Cell_Index())=velocity_vector(axis);}}}
    
    LAPLACE_CALLBACKS_INCOMPRESSIBLE<T_GRID> incompressible_callbacks(*this);
    for(int axis=1;axis<=TV::dimension;axis++){
        Unpack_And_Interpolate_And_Inject_Cell_Values(laplace_incompressible,incompressible_callbacks,time,cell_velocity_components_pointer(axis),cell_velocity_boundary_components(axis),0);
        WRITE_CELL_ARRAY_POINTER(cell_velocity_components_pointer(axis),STRING_UTILITIES::string_sprintf("new velocities %d",axis));}

    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            TV dx=laplace_grid.Grid(grid_index).DX();
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                D_FACE_INDEX face_index=iterator.Full_Index();

                T& face_velocity=face_velocities(grid_index)->operator()(face_index);
                //T& face_velocity_update=face_velocities_update(grid_index)(face_index);
                if(Is_NaN(face_velocity)){
                    TV normal=laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(iterator.Axis()));
                    TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
                    TV velocity_vector;
                    for(int axis=1;axis<=TV::dimension;axis++)
                        velocity_vector(axis)=(T).5*(cell_velocity_components(axis)(grid_index)(iterator.First_Cell_Index())+cell_velocity_components(axis)(grid_index)(iterator.Second_Cell_Index()));
                    face_velocity=TV::Dot_Product(normal,velocity_vector);
                }}}
    
    //the remaining faces are interpolated through grid coupling
    Write_Substep("after computing unsolved velocities in voronoi cells",0,1);
}
//#####################################################################
// Function Set_Zero_Air_Velocities
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Set_Zero_Air_Velocities()
{
    if(!example.fluids_parameters.water) return;
    int number_of_ghost_cells=example.fluids_parameters.number_of_ghost_cells;
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index),number_of_ghost_cells);iterator.Valid();iterator.Next()){
                D_FACE_INDEX face_index=iterator.Full_Index();
                TV_INT first_cell_index=iterator.First_Cell_Index();
                TV_INT second_cell_index=iterator.Second_Cell_Index();
                
                if(!((laplace_grid.Grid(grid_index).Domain_Indices(number_of_ghost_cells).Lazy_Inside(first_cell_index) && Phi(grid_index,first_cell_index)<=0) ||
                   (laplace_grid.Grid(grid_index).Domain_Indices(number_of_ghost_cells).Lazy_Inside(second_cell_index) && Phi(grid_index,second_cell_index)<=0)))
                    FACE_VELOCITIES(laplace_grid.Local_Grid_Index(grid_index))(face_index)=0;}
    Write_Substep("after set zero air velocities",0,1);
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface_On_All_Local_Grids
//#####################################################################
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::
Extrapolate_Velocity_Across_Interface_On_All_Local_Grids()
{
    /*if(!example.fluids_parameters.water) return;
    int number_of_ghost_cells=example.fluids_parameters.number_of_ghost_cells;

    above_water_face_masks.Resize(n_local_grids);

    FOR_EACH_LOCAL_GRID(grid_index)
        Extrapolate_Velocity_Across_Interface(LOCAL_GRID_ACCESS(grid_index),FACE_VELOCITIES(grid_index),PHI(grid_index),0,(T)number_of_ghost_cells,0,TV(),0);
    Write_Substep("after extrapolate velocity across interface",0,1);

    // fill velocities
    ARRAY<T_FACE_ARRAYS_SCALAR> face_velocities_ghost(n_local_grids);
    FOR_EACH_LOCAL_GRID(grid_index){
        face_velocities_ghost(grid_index).Resize(LOCAL_GRID_ACCESS(grid_index).Domain_Indices(number_of_ghost_cells));}
    example.Fill_Ghost_Cells_Face_Chimera(number_of_ghost_cells,&face_velocities_ghost);
    // fill levelset
    example.Fill_Ghost_Cells_Chimera(number_of_ghost_cells);
    // debug output
    FOR_EACH_LOCAL_GRID(grid_index)
        T_FACE_ARRAYS_SCALAR::Exchange_Arrays(FACE_VELOCITIES(grid_index),face_velocities_ghost(grid_index));
    Write_Substep("these are the ghost velocities",0,1);
    FOR_EACH_LOCAL_GRID(grid_index)
        T_FACE_ARRAYS_SCALAR::Exchange_Arrays(FACE_VELOCITIES(grid_index),face_velocities_ghost(grid_index));

    FOR_EACH_LOCAL_GRID(grid_index){
        Extrapolate_Velocity_Across_Interface(LOCAL_GRID_ACCESS(grid_index),face_velocities_ghost(grid_index),PHI(grid_index),(laplace_grid.Global_Grid_Index(grid_index)==1?0:(T)number_of_ghost_cells),(T)number_of_ghost_cells,0,TV(),0);
        T_FACE_ARRAYS_SCALAR::Copy(face_velocities_ghost(grid_index),FACE_VELOCITIES(grid_index));}
    Write_Substep("after extrapolate velocity across interface with ghost",0,1);*/
}
//#####################################################################
// Function Print_Matrices
//#####################################################################
template<class T_GRID> void Print_Matrices(const SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution)
{
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;

    LOG::cout << "rhs " << evolution.rhs_incompressible << std::endl;
    LOG::cout << "divergence" << std::endl << evolution.divergence_matrix << std::endl;
    LOG::cout << "mass" << std::endl << evolution.inverse_mass_matrix << std::endl;

    LOG::cout << "cells " << evolution.n_matrix_cells << " faces " << evolution.n_matrix_faces << " kinematic faces " << evolution.n_kinematic_faces << std::endl;

    for(int grid_index=1;grid_index<=evolution.n_global_grids;grid_index++)
        if(evolution.laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(evolution.laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next())
            {
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=evolution.laplace_incompressible.Matrix_Cell_Index(grid_index,cell_index);
                LOG::cout << "cell " << grid_index << " " << cell_index << " " << matrix_cell_index << std::endl;
            }

    for(int grid_index=1;grid_index<=evolution.n_global_grids;grid_index++)
        if(evolution.Local_Grid(grid_index))
            for(FACE_ITERATOR iterator(evolution.laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next())
            {
                D_FACE_INDEX face_index=iterator.Full_Index();
                int matrix_face_index=evolution.Matrix_Face_Index(grid_index,face_index);
                LOG::cout << "face " << grid_index << " " << face_index.axis << " " << face_index.index << " " << matrix_face_index << std::endl;
            }

    for(int voronoi_index=1;voronoi_index<=evolution.laplace_grid.voronoi_faces.Size();voronoi_index++)
    {
        int matrix_face_index=evolution.laplace_incompressible.Matrix_Face_Index(voronoi_index);
        LOG::cout << "voronoi face " << evolution.laplace_grid.voronoi_faces(voronoi_index).x(1).x << " " << evolution.laplace_grid.voronoi_faces(voronoi_index).x(1).y << " " << evolution.laplace_grid.voronoi_faces(voronoi_index).x(2).x << " " << evolution.laplace_grid.voronoi_faces(voronoi_index).x(2).y << " " << matrix_face_index << std::endl;
    }

    LOG::cout << "velocities " << evolution.face_velocities_packed << std::endl;
}
//#####################################################################
// Function Compute_Face_Dirichlet_Update
//#####################################################################
template<class T_GRID> class LAPLACE_CALLBACKS_VORTICITY:public LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
public:
    SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution;
    
    LAPLACE_CALLBACKS_VORTICITY(SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution_input):
        evolution(evolution_input){}
    
    bool Get_Neumann_Boundary_Condition(const TV& location,const TV& normal,T& normal_gradient,const T time)
    {return false;}
    bool Get_Dirichlet_Boundary_Condition(const TV& location,T& value,const T time)
    {return false;}
    T Get_Density(const TV& location)
    {return 0;}
};
template<class T_GRID>
void Compute_Vorticity_Helper(SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>& evolution,const typename T_GRID::SCALAR time,const typename T_GRID::SCALAR dt) {}
template<class T>
void Compute_Vorticity_Helper(SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<GRID<VECTOR<T,2> > >& evolution,const T time,const T dt)
{
    //LOG::cout << std::setprecision(16);

    typedef GRID<VECTOR<T,2> > T_GRID;
    typedef typename T_GRID::VECTOR_T TV;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef NEWMARK_EVOLUTION<TV> BASE;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef MATRIX<T,TV::dimension,TV::dimension> MATRIX_TV;

#define EFV(I) evolution.example.incompressible_fluid_containers(I)->face_velocities

    //int n_global_grids=evolution.example.chimera_grid->number_of_global_grids;
    
    /*ARRAY<VECTOR<GRID_CELL_INDEX,TV::dimension+1> >& simplex_indices=evolution.laplace_grid.interpolation.simplex_indices;
    int n_simplices=simplex_indices.Size();
    ARRAY<T> simplex_curl(n_simplices);

    HASHTABLE<VORONOI_FACE_INDICES

    for(int simplex_index=1;simplex_index<=n_simplices;simplex_index++){
        for(int i=1;i<=TV::dimension+1;i++){
            VECTOR<int,2> edge(simplex_indices(simplex_index)(i),simplex_indices(simplex_index)(i%(TV::dimension+1)+1));
            
    }*/

    /*ARRAY<T_FACE_ARRAYS_SCALAR> face_tangential_velocities(n_global_grids);
    for(int grid_index=1;grid_index<=evolution.n_global_grids;grid_index++)
        if(evolution.laplace_grid.Local_Grid(grid_index)){
            face_tangential_velocities(grid_index).Resize(evolution.laplace_grid.Grid(grid_index).Domain_Indices());
            for(FACE_ITERATOR iterator(evolution.laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next())
                if(evolution.laplace_grid.Chimera_Face(grid_index,iterator.Full_Index())){
                    int axis=iterator.Axis();
                    int axis2=3-axis;
                    TV_INT index=iterator.Face_Index();
                    face_tangential_velocities(grid_index)(iterator.Full_Index())=0.25*(EFV(grid_index).Component(axis2)(index)+EFV(grid_index).Component(axis2)(index+TV_INT::Axis_Vector(axis2))+EFV(grid_index).Component(axis2)(index-TV_INT::Axis_Vector(axis))+EFV(grid_index).Component(axis2)(index-TV_INT::Axis_Vector(axis)+TV_INT::Axis_Vector(axis2)));}}*/
    
    /*if(evolution.voronoi_face_velocities.Size()==evolution.laplace_grid.voronoi_faces.Size()){
    evolution.laplace_grid.voronoi.voronoi_face_velocities.Resize(evolution.laplace_grid.voronoi_faces.Size());
    VECTOR_ND<T> voronoi_face_tangential_velocities(evolution.laplace_grid.voronoi_faces.Size());
    for(int voronoi_face_index=1;voronoi_face_index<=evolution.laplace_grid.voronoi_faces.Size();voronoi_face_index++)
        if(evolution.laplace_incompressible.Matrix_Face_Index(voronoi_face_index)){
            VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=evolution.laplace_grid.voronoi_faces(voronoi_face_index).x;
            TV location_1=evolution.laplace_grid.Frame(grid_cell_indices(1).x)*evolution.laplace_grid.Grid(grid_cell_indices(1).x).X(grid_cell_indices(1).y);
            TV location_2=evolution.laplace_grid.Frame(grid_cell_indices(2).x)*evolution.laplace_grid.Grid(grid_cell_indices(2).x).X(grid_cell_indices(2).y);
            TV location=(T).5*(location_1+location_2);
            TV normal=(location_2-location_1).Normalized();
            TV tangent=normal.Perpendicular();
            evolution.laplace_grid.voronoi.voronoi_face_velocities(voronoi_face_index)=normal*evolution.voronoi_face_velocities(voronoi_face_index);}}*/
    /*for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(evolution.laplace_grid.Local_Grid(grid_index)){
            evolution.vorticity(grid_index).Resize(evolution.laplace_grid.Grid(grid_index).Domain_Indices(1));
            for(CELL_ITERATOR iterator(evolution.laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(evolution.laplace_grid.Chimera_Cell(grid_index,cell_index)){
                    T curl=0;
                    for(int axis=1;axis<=TV::dimension;axis++)
                        for(int side=1;side<=TV::dimension;side++){
                            D_FACE_INDEX face_index(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis));
                            if(evolution.laplace_grid.Chimera_Face(grid_index,face_index)){
                                if(grid_index==2 && cell_index(1)==32 && (cell_index(2)==26 || cell_index(2)==25))
                                    LOG::cout << "cartesian face tangent " << face_index.axis << " " << face_index.index << " " << face_tangential_velocities(grid_index)(face_index) << std::endl;
                                curl+=evolution.laplace_grid.Grid(grid_index).DX()(3-axis)*(3-2*axis)*(3-2*side)*face_tangential_velocities(grid_index)(face_index);}}
                    GRID_CELL_INDEX gci(grid_index,cell_index);
                    if(evolution.laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(gci)){
                        const ARRAY<int>& ii=evolution.laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(gci);
                        for(int i=1;i<=ii.Size();i++){
                            if(grid_index==2 && cell_index(1)==32 && (cell_index(2)==26 || cell_index(2)==25))
                                LOG::cout << "voronoi face tangent " << evolution.laplace_grid.voronoi_faces(ii(i)).x(1).x << " " << evolution.laplace_grid.voronoi_faces(ii(i)).x(1).y << " - " << evolution.laplace_grid.voronoi_faces(ii(i)).x(2).x << " " << evolution.laplace_grid.voronoi_faces(ii(i)).x(2).y << " : " <<  voronoi_face_tangential_velocities(ii(i)) << " " << evolution.laplace_grid.voronoi.voronoi_face_velocities(ii(i)) << " " << evolution.laplace_grid.voronoi_faces(ii(i)).y << " " << voronoi_face_velocities(ii(i)) << std::endl;
                            curl+=((evolution.laplace_grid.voronoi_faces(ii(i)).x(1)==gci)?-1:1)*evolution.laplace_grid.voronoi_faces(ii(i)).y*voronoi_face_tangential_velocities(ii(i));}}
                    if(grid_index==2 && cell_index(1)==32 && (cell_index(2)==26 || cell_index(2)==25))
                        LOG::cout << "curl " << curl << " " << evolution.laplace_grid.Cell_Size(grid_index,cell_index) << " " << evolution.laplace_grid.Grid(grid_index).DX() << std::endl;
                        evolution.vorticity(grid_index)(cell_index)=curl/evolution.laplace_grid.Cell_Size(grid_index,cell_index);}}}
    
    LAPLACE_CALLBACKS_VORTICITY<T_GRID> callbacks(evolution);
    ARRAY<T_ARRAYS_SCALAR*> arrays(n_global_grids);
    for(int i=1;i<=n_global_grids;i++)
        arrays(i)=&evolution.vorticity(i);
    ARRAY<ARRAY<T> > vorticity_boundary(n_global_grids);
    evolution.Unpack_And_Interpolate_And_Inject_Cell_Values(evolution.laplace_incompressible,callbacks,time+dt,arrays,vorticity_boundary,0);*/
}
template<class T_GRID> void SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>::Compute_Vorticity(const T time,const T dt)
{
    Compute_Vorticity_Helper<T_GRID>(*this,time,dt);
}
//#####################################################################
template class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<GRID<VECTOR<float,1> > >;
template class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<GRID<VECTOR<float,2> > >;
template class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<GRID<VECTOR<double,1> > >;
template class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<GRID<VECTOR<double,2> > >;
template class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<GRID<VECTOR<double,3> > >;
#endif

