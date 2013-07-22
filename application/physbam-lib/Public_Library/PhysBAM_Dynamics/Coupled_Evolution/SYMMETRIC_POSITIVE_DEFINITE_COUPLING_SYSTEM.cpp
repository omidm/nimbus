//#####################################################################
// Copyright 2009, Nipun Kwatra, Jon Gretarsson, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_BASE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_CUT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/GENERALIZED_FLUID_MASS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_POISSON.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_FORCES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_VISCOUS_FORCES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
namespace PhysBAM{template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);}

namespace PhysBAM{
template<class TV> void Dump_Extra_Velocities(FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<typename TV::SCALAR,1> > &cut,const VECTOR_ND<typename TV::SCALAR>& fluid_velocity_vector)
{PHYSBAM_FATAL_ERROR();}
template<class TV> void Dump_Extra_Velocities(FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<typename TV::SCALAR,2> > &cut,const VECTOR_ND<typename TV::SCALAR>& fluid_velocity_vector)
{cut.Dump_Extra_Velocities(fluid_velocity_vector);}
template<class TV> void Dump_Extra_Velocities(FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<typename TV::SCALAR,3> > &cut,const VECTOR_ND<typename TV::SCALAR>& fluid_velocity_vector)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM(const bool use_preconditioner_input,BACKWARD_EULER_SYSTEM<TV>* solid_system_input,
    IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>& boundary_condition_collection,UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,
    const SOLID_BODY_COLLECTION<TV>& solid_body_collection,
    const INCOMPRESSIBLE_FLUID_CONTAINER<GRID<TV> >& incompressible_fluid_container,const ARRAY<T,TV_INT>& density_input,const ARRAY<TV,TV_INT>& centered_velocity_input,
    const ARRAY<T,TV_INT>& one_over_rho_c_squared_input,const bool leakproof_solve_input,const bool using_slip,const bool using_viscosity)
    :BASE(use_preconditioner_input,true),outside_fluid(boundary_condition_collection.psi_D),inactive_faces(boundary_condition_collection.psi_N),
    density(density_input),centered_velocity(centered_velocity_input),
    one_over_rho_c_squared(one_over_rho_c_squared_input),
    solid_system(solid_system_input),index_map(*new COLLISION_AWARE_INDEX_MAP<TV>(info,boundary_condition_collection)),
    solid_interpolation(0),fluid_gradient(0),fluid_poisson(0),solid_forces(0),fluid_interpolation(0),fluid_mass(0),fluid_viscous_forces(0),fluid_to_solid_interpolation(0),
    leakproof_solve(leakproof_solve_input),dt(0),run_self_tests(false),print_poisson_matrix(false),print_matrix(false),print_rhs(false),
    print_each_matrix(false),print_index_map(false),use_viscous_forces(using_viscosity),surface_tension_coefficient(0),
    solid_node(true),fluid_node(true),use_full_ic(false),mpi_solid_fluid(0),mpi_grid(0),debug_velocity(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
~SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM()
{
    delete &index_map;
    delete fluid_gradient;
    delete fluid_poisson;
    delete solid_forces;
    delete solid_interpolation;
    delete fluid_interpolation;
    delete fluid_mass;
    delete fluid_viscous_forces;
    delete fluid_to_solid_interpolation;
}
//#####################################################################
// Function Resize_Coupled_System_Vector
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Resize_Coupled_System_Vector(VECTOR_T& b) const
{
    b.pressure.Resize(index_map.Number_Cells());
    b.lambda.Resize(fluid_interpolation->Number_Of_Constraints());
    if(!leakproof_solve && solid_node) b.force_coefficients.Resize(solid_forces->Velocity_Dependent_Forces_Size());
    else b.force_coefficients.Resize(FORCE_AGGREGATE_ID());
    if(use_viscous_forces && fluid_node) b.viscous_force_coefficients.Resize(fluid_viscous_forces->Viscous_Forces_Size());
    else b.viscous_force_coefficients.Resize(VISCOUS_FORCE_ID());
}
//#####################################################################
// Function Get_Pressure
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Get_Pressure(const VECTOR_T& V,ARRAY<T,TV_INT>& fluid_pressures) const
{
    index_map.Distribute_Indexed_Cells(V.pressure,fluid_pressures);
    if(dt) fluid_pressures*=1/dt;
}
//#####################################################################
// Function Zero_Coupling_Faces_Values
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Zero_Coupling_Faces_Values(T_FACE_ARRAYS_SCALAR& face_array) const
{
    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(index_map.iterator_info);iterator.Valid();iterator.Next())
        face_array(iterator.Full_Index())=0;
}
//#####################################################################
// Function Apply_Lambda_To_Euler_State
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Apply_Lambda_To_Euler_State(const VECTOR_T& V,const ARRAY<T,COUPLING_CONSTRAINT_ID>& coupled_faces_solid_interpolated_velocity_n,
    const ARRAY<T,COUPLING_CONSTRAINT_ID>& coupled_faces_solid_interpolated_velocity_np1,const ARRAY<T,TV_INT>& cell_volume,ARRAY<TV_DIMENSION,TV_INT>& U) const
{
    VECTOR_ND<T> impulse_at_coupling_faces(index_map.Number_Faces());
    fluid_interpolation->Transpose_Times(V.lambda,impulse_at_coupling_faces);
    int ghost_cells=0;
    typename GRID<TV>::REGION region_type=ghost_cells?GRID<TV>::INTERIOR_REGION:GRID<TV>::WHOLE_REGION;
    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(index_map.iterator_info,ghost_cells,region_type);iterator.Valid();iterator.Next()){
        TV_INT fluid_cell_index=iterator.Real_Cell_Index();
        T solid_interpolated_velocity_average=(coupled_faces_solid_interpolated_velocity_n(COUPLING_CONSTRAINT_ID(iterator.collision_index))+coupled_faces_solid_interpolated_velocity_np1(COUPLING_CONSTRAINT_ID(iterator.collision_index)))*(T).5;
        T impulse=impulse_at_coupling_faces(index_map.indexed_faces.m+index_map.constraint_indices.Get(SIDED_FACE_INDEX<TV::dimension>(iterator.side,iterator.Full_Index())))/cell_volume(iterator.Real_Cell_Index());

        U(fluid_cell_index)(iterator.Axis()+1)+=impulse; // momentum update
        U(fluid_cell_index)(TV::dimension+2)+=impulse*solid_interpolated_velocity_average;} // energy update
}
//#####################################################################
// Function Test_Viscosity
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Test_Viscosity(const VECTOR_T& V,const VECTOR_T& B) const
{
    ARRAY<T,VISCOUS_FORCE_ID> temporary_forces=V.viscous_force_coefficients;
    fluid_viscous_forces->Transpose_Times(V.viscous_force_coefficients,temporary_viscous_velocities);
    fluid_mass->Inverse_Times(temporary_viscous_velocities,temporary_viscous_velocities);
    fluid_viscous_forces->Times_Add(temporary_viscous_velocities,temporary_forces);
    temporary_forces-=B.viscous_force_coefficients;
    std::stringstream ss;
    ss<<"Printing viscous errors"<<std::endl;
    ss<<temporary_forces<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Apply_Velocity_Update
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Apply_Velocity_Update(const VECTOR_T& V,ARRAY<T,FACE_INDEX<TV::dimension> >& fluid_velocity,ARRAY<T,TV_INT>& fluid_pressures,GENERALIZED_VELOCITY<TV>& solid_velocity,
    GENERALIZED_VELOCITY<TV>& force_on_solid,bool want_solid,bool want_fluid) const
{
    GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);
    Gather(V,temporary_faces,temporary_solids_velocity);
    Exchange_Velocities(temporary_faces,temporary_solids_velocity);
    ARRAY<T> dummy_array(index_map.indexed_constraints.m);

    if(want_solid) force_on_solid.Copy(-1/dt,temporary_solids_velocity);

    Inverse_Mass(temporary_faces,temporary_solids_velocity);

    if(want_solid){
        if(!mpi_solid_fluid){
            GENERALIZED_VELOCITY<TV> pressure_force(pressure_impulses,pressure_impulses_twist,solid_system->solid_body_collection);
            solid_interpolation->Transpose_Times(V.lambda,pressure_force);} // Save the pressure force for fracture... TODO(jontg): Requires an additional MPI exchange
        solid_velocity-=temporary_solids_velocity;}
    if(mpi_solid_fluid) mpi_solid_fluid->Distribute_Lists_From_Solid_Node(solid_velocity);

    if(want_fluid){
        temporary_faces=-temporary_faces;
        Add_Dirichlet_Pressures_To_Velocity(fluid_pressures,temporary_faces);
        Add_Surface_Tension(temporary_faces);
        Apply_Massless_Structure_Force_To_Fluid(temporary_faces,FLT_MAX);
        if(fluid_to_solid_interpolation)
            if(FLUID_TO_SOLID_INTERPOLATION_CUT<TV>* cut=dynamic_cast<FLUID_TO_SOLID_INTERPOLATION_CUT<TV>*>(fluid_to_solid_interpolation))
                Dump_Extra_Velocities<TV>(*cut,temporary_faces);
        index_map.Distribute(temporary_faces,fluid_velocity,dummy_array);

        if(run_self_tests) Test_Incompressibility(fluid_velocity,dummy_array);
        index_map.Distribute(V.pressure,fluid_pressures);
        if(dt) fluid_pressures*=1/dt;
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(index_map.grid);it.Valid();it.Next())
            if(!index_map.face_indices(it.Full_Index()))
                fluid_velocity(it.Full_Index())=0;}

    temporary_faces*=(T)0;
    Dump_Substep(temporary_faces,"before pressure");
//    INTERPOLATED_COLOR_MAP<T> color_map;
//    color_map.Initialize_Colors(V.pressure.Min(),V.pressure.Max(),false,false,false);
//    for(int i=1;i<=index_map.indexed_cells.m;i++){
//        Add_Debug_Particle(index_map.grid.X(index_map.indexed_cells(i)),color_map(V.pressure(i)));}
    Dump_Substep(temporary_faces,"pressure");

    if(fluid_to_solid_interpolation){
        index_map.Collect(fluid_velocity,dummy_array,temporary_faces);
        Fill_Extra_Velocities(temporary_faces);
        fluid_to_solid_interpolation->Times(temporary_faces,solid_velocity);
        Dump_Substep(temporary_faces,"after interpolate");}
}
//#####################################################################
// Function Interpolate_Solid_Velocity_To_Coupled_Faces
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Interpolate_Solid_Velocity_To_Coupled_Faces(const GENERALIZED_VELOCITY<TV>& solids_velocity,ARRAY<T,COUPLING_CONSTRAINT_ID>& coupled_faces_solid_interpolated_velocity)
{
    // Jv
    coupled_faces_solid_interpolated_velocity.Resize(solid_interpolation->Number_Of_Constraints());
    solid_interpolation->Times(solids_velocity,coupled_faces_solid_interpolated_velocity);
}
//#####################################################################
// Function Apply_One_Sided_Interpolation_At_Coupling_Faces
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Apply_One_Sided_Interpolation_At_Coupling_Faces(const T_FACE_ARRAYS_BOOL& psi_N_domain_boundary,
    const bool use_one_sided_face_velocty_interpolation,T_FACE_ARRAYS_SCALAR& fluids_velocity)
{
    beta_face.Resize(index_map.grid);
    constrained_beta_face.Resize(index_map.indexed_constraints.m);
    constrained_fluid_velocity.Resize(index_map.indexed_constraints.m);

    for(int i=1;i<=index_map.indexed_faces.m;i++){
        FACE_INDEX<TV::dimension> face_index=index_map.indexed_faces(i);
        TV_INT first_cell_index=face_index.First_Cell_Index(),second_cell_index=face_index.Second_Cell_Index();
        beta_face(face_index)=Inverse((density(first_cell_index)+density(second_cell_index))*(T).5);}

    // Go through coupling faces and use density and velocity only on the fluid side
    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(index_map.iterator_info);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Real_Cell_Index();
        int constraint_number=index_map.constraint_indices.Get(SIDED_FACE_INDEX<TV::dimension>(iterator.side,iterator.Full_Index()));
        constrained_beta_face(constraint_number) = Inverse(density(cell_index));
        if(use_one_sided_face_velocty_interpolation) constrained_fluid_velocity(constraint_number) = centered_velocity(cell_index)[iterator.Axis()];
        else constrained_fluid_velocity(constraint_number) = fluids_velocity(iterator.Full_Index());}
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Compute(int ghost_cells,const T dt_input,const T current_velocity_time,const T_FACE_ARRAYS_BOOL& psi_N_domain_boundary,
    const bool disable_thinshell,const bool use_one_sided_face_velocty_interpolation,
    T_FACE_ARRAYS_SCALAR& fluids_velocity,T mu,bool use_second_order_cut_cell,const T_LEVELSET* levelset_input)
{
    solve_id++;
    if(!solid_interpolation) solid_interpolation=new MATRIX_SOLID_INTERPOLATION<TV>(index_map.iterator_info);
    if(!fluid_gradient) fluid_gradient=new MATRIX_FLUID_GRADIENT<TV>(index_map);
    if(!fluid_poisson) fluid_poisson=new MATRIX_FLUID_POISSON<TV>(index_map,one_over_rho_c_squared);
    if(!solid_forces) solid_forces=new MATRIX_SOLID_FORCES<TV>(solid_system->solid_body_collection);
    if(!fluid_interpolation) fluid_interpolation=new MATRIX_FLUID_INTERPOLATION<TV>(index_map);
    if(!fluid_mass) fluid_mass=new GENERALIZED_FLUID_MASS<TV>(index_map,beta_face,constrained_beta_face);
    if(!fluid_viscous_forces) fluid_viscous_forces=new MATRIX_VISCOUS_FORCES<TV>(index_map);


    // TODO(kwatra): make robust. Also make sure thickness is *same* as that used for Inside() tests for setting dirichlet cells in boundaries.
    // Otherwise multiple neumann faces without dirichlet cells can occur.
    dt=dt_input;
    solid_mass=&solid_system->projection_data.mass;

    if(fluid_node){
        levelset=levelset_input;
        index_map.iterator_info.Initialize_Collision_Aware_Face_Iterator(outside_fluid,inactive_faces,0,disable_thinshell);
        // TODO(jontg): Fluid-Fluid MPI coupling
        index_map.Construct_Indices(ghost_cells);
        Apply_One_Sided_Interpolation_At_Coupling_Faces(psi_N_domain_boundary,use_one_sided_face_velocty_interpolation,fluids_velocity);
        if(index_map.two_phase && !fluid_to_solid_interpolation && use_second_order_cut_cell) fluid_mass->Compute_For_Two_Phase_Pressure_Jump(levelset->phi,density);
        else{
            fluid_mass->Compute();
            if(use_second_order_cut_cell) fluid_mass->Second_Order_Mass_Correction(levelset->phi);}
        fluid_gradient->Compute(psi_N_domain_boundary);
        if(!fluid_to_solid_interpolation) fluid_interpolation->Compute(0);
        if(use_viscous_forces) fluid_viscous_forces->Compute(dt,psi_N_domain_boundary,mu);
        if(fluid_to_solid_interpolation) fluid_to_solid_interpolation->Compute(0);
        fluid_poisson->Compute(fluid_gradient->gradient,fluid_mass->one_over_fluid_mass_at_faces,dt,use_preconditioner);

        if(dt){T one_over_dt_squared=(T)1/(dt*dt);
            one_over_rho_c_squared_flat.Resize(index_map.Number_Cells());
            for(int i=1;i<=index_map.indexed_cells.m;i++){
                one_over_rho_c_squared_flat(i)=one_over_rho_c_squared(index_map.indexed_cells(i))*one_over_dt_squared*index_map.grid.Cell_Size();}}}

    if(!fluid_to_solid_interpolation) solid_interpolation->Compute(0);

    if(!leakproof_solve){
        temporary_velocities.Resize(solid_system->solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        temporary_twists.Resize(solid_system->solid_body_collection.rigid_body_collection.rigid_geometry_collection.particles.array_collection->Size());
        pressure_impulses.Resize(solid_system->solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        pressure_impulses_twist.Resize(solid_system->solid_body_collection.rigid_body_collection.rigid_geometry_collection.particles.array_collection->Size());}
    if(!leakproof_solve && solid_node){
        solid_forces->Compute(dt,current_velocity_time);}

    if(use_preconditioner) Compute_Lambda_Diagonal_Preconditioner();
    if(use_full_ic) Compute_Full_Preconditioner();
    PHYSBAM_ASSERT(mpi_solid_fluid || solid_interpolation->Number_Of_Constraints()==fluid_interpolation->Number_Of_Constraints());
}
//#####################################################################
// Function Add_Dirichlet_Pressures_To_Velocity
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Add_Dirichlet_Pressures_To_Velocity(const ARRAY<T,TV_INT>& pressure,VECTOR_ND<T>& fluid_velocity_vector) const
{
    for(int i=1;i<=fluid_gradient->ghost_gradient.m;i++){const typename MATRIX_FLUID_GRADIENT_BASE<TV>::GHOST_GRADIENT_ENTRY& ge=fluid_gradient->ghost_gradient(i);
        fluid_velocity_vector(ge.face)-=fluid_mass->one_over_fluid_mass_at_faces(ge.face)*ge.weight*pressure(ge.index)*dt;}
}
//#####################################################################
// Function Add_Surface_Tension
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Add_Surface_Tension(VECTOR_ND<T>& fluid_velocity_vector) const
{
    if(!levelset || !surface_tension_coefficient) return;
    GRID<TV>& grid=index_map.grid;
    
    T mx=0,av=0;
    int n=0;
    for(int i=1;i<=fluid_gradient->interface_gradient.m;i++){const typename MATRIX_FLUID_GRADIENT_BASE<TV>::INTERFACE_ENTRY& ie=fluid_gradient->interface_gradient(i);
        FACE_INDEX<TV::m> face_index=index_map.indexed_faces(ie.face);
        TV_INT cell1=face_index.First_Cell_Index(),cell2=face_index.Second_Cell_Index();
        T phi1=levelset->phi(cell1),phi2=levelset->phi(cell2),theta;
        if(phi1!=phi2) theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
        else theta=0.5;// use the boundary location when evaluate across boundary
        TV X1=grid.X(cell1),X2=grid.X(cell2),X=X1+(X2-X1)*theta;
        PHYSBAM_ASSERT(0<=theta && theta<=1);
        T jump=-surface_tension_coefficient*levelset->Compute_Curvature(X);
        mx=std::max(mx,abs(jump/surface_tension_coefficient-100));
        av+=abs(jump/surface_tension_coefficient-100);
        n++;
        fluid_velocity_vector(ie.face)-=fluid_mass->one_over_fluid_mass_at_faces(ie.face)*ie.weight*jump*dt;}
    if(n) {std::stringstream ss;ss<<"asdfasdf "<<mx<<"   "<<av/n<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Set_Up_RHS
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Set_Up_RHS(VECTOR_T& V,VECTOR_T& F,const GENERALIZED_VELOCITY<TV>& solids_velocity_star,const ARRAY<T,FACE_INDEX<TV::dimension> >& fluids_velocity_star,
    const ARRAY<T,TV_INT>& p_advected_over_rho_c_squared_dt,const ARRAY<T,TV_INT>& p_advected,const ARRAY<T,TV_INT>& fluid_pressures)
{
    Resize_Coupled_System_Vector(tolerances);
    Resize_Coupled_System_Vector(F);
    Resize_Coupled_System_Vector(V);
    coupling_faces.Resize(fluid_interpolation->Number_Of_Constraints());
    VECTOR_ND<T> fluid_velocity_vector;
    if(fluid_node){
        temporary_faces.Resize(index_map.Number_Faces());
        temporary_lambdas.Resize(F.lambda.Size());
        pressure.Resize(V.pressure.n);
        temporary_viscous_velocities.Resize(index_map.Number_Faces());
        index_map.Collect(fluids_velocity_star,constrained_fluid_velocity,fluid_velocity_vector);}

    Exchange_Velocities(temporary_faces,const_cast<GENERALIZED_VELOCITY<TV>&>(solids_velocity_star));
    Fill_Extra_Velocities(fluid_velocity_vector);
    Dump_Substep(fluid_velocity_vector,"before surface tension");
    Add_Dirichlet_Pressures_To_Velocity(fluid_pressures,fluid_velocity_vector);
    Add_Surface_Tension(fluid_velocity_vector);
    Dump_Substep(fluid_velocity_vector,"add surface tension");
    Apply_Massless_Structure_Force_To_Fluid(fluid_velocity_vector,FLT_MAX);
    Dump_Substep(fluid_velocity_vector,"after add solid");
    Scatter(fluid_velocity_vector,solids_velocity_star,F);
    if(fluid_node){
        VECTOR_ND<T> p_advected_over_rho_c_squared_dt_vector;
        index_map.Collect(p_advected_over_rho_c_squared_dt,p_advected_over_rho_c_squared_dt_vector);
        F.pressure+=p_advected_over_rho_c_squared_dt_vector*index_map.grid.Cell_Size();}

    Setup_Tolerances(F,fluid_velocity_vector,solids_velocity_star);
    Setup_Initial_Guess(F,V,p_advected);
    Project(F);

    if(print_index_map && !print_matrix) index_map.Print(solve_id);
    if(print_matrix) Print_Matrix(V);
    if(print_rhs) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",F);
    if(print_each_matrix){
        std::stringstream ss;
        ss<<"Solve id: "<<solve_id<<std::endl;
        LOG::filecout(ss.str());
        Print_Each_Matrix(solve_id);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("us-%i.txt",solve_id).c_str()).Write("us",fluid_velocity_vector);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("vs-%i.txt",solve_id).c_str()).Write("vs",solids_velocity_star);}
}
//#####################################################################
// Function Setup_Tolerances
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Setup_Tolerances(const VECTOR_T& F,const VECTOR_ND<T>& fluid_velocity,const GENERALIZED_VELOCITY<TV>& structure_velocity)
{
    GENERALIZED_VELOCITY<TV> solids_velocity_star_projected(temporary_velocities,temporary_twists,solid_system->solid_body_collection);

    T eps=(T)1e-5;
    if(solid_node){
        for(FORCE_AGGREGATE_ID i(1);i<=F.force_coefficients.Size();i++)
            tolerances.force_coefficients(i)=max(eps,(T)1e-2*abs(F.force_coefficients(i)));
        solids_velocity_star_projected=structure_velocity;
        solid_system->Project(solids_velocity_star_projected);} // TODO: Exchange

    if(mpi_solid_fluid) mpi_solid_fluid->Distribute_Lists_From_Solid_Node(solids_velocity_star_projected);
    if(fluid_node){
        fluid_gradient->Collect_Maxabs_Velocity(fluid_velocity,tolerances.pressure);
        T eps_velocity=eps*index_map.grid.Face_Sizes()[1];
        T scaling_factor=(T)1e-5;
        tolerances.pressure*=scaling_factor;
        for(int i=1;i<=index_map.real_cell_indices.m;++i){int index=index_map.real_cell_indices(i);
            if(tolerances.pressure(index)<eps_velocity) tolerances.pressure(index)=eps_velocity;}

        for(COUPLING_CONSTRAINT_ID i(1);i<=F.lambda.Size();i++)
            tolerances.lambda(i)=max((T)1e-3,abs(F.lambda(i)));
    

        solid_interpolation->Times(solids_velocity_star_projected,temporary_lambdas);
        for(COUPLING_CONSTRAINT_ID i(1);i<=temporary_lambdas.Size();i++)
            tolerances.lambda(i)=max(eps,scaling_factor*max(tolerances.lambda(i),abs(temporary_lambdas(i))));

        if(use_viscous_forces){
            for(VISCOUS_FORCE_ID i(1);i<=F.viscous_force_coefficients.Size();i++)
                tolerances.viscous_force_coefficients(i)=max(eps,(T)1e-2*abs(F.viscous_force_coefficients(i)));}}
    tolerances*=(T)1e-5;
}
//#####################################################################
// Function Setup_Initial_Guess
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Setup_Initial_Guess(const VECTOR_T& F,VECTOR_T& V,const ARRAY<T,TV_INT>& p_advected) const
{
    if(fluid_node){
        // Initial guess
        index_map.Collect(p_advected,V.pressure);
        V.pressure*=dt;
        for(int i=1;i<=index_map.indexed_faces.m;i++){
            FACE_INDEX<TV::dimension> face_index=index_map.indexed_faces(i);
            T face_area=index_map.grid.Face_Size(face_index.axis);
            if(!(*index_map.iterator_info.outside_fluid)(face_index.First_Cell_Index())){
                temporary_faces(i)=-p_advected(face_index.First_Cell_Index())*dt*face_area;}
            else if(!(*index_map.iterator_info.outside_fluid)(face_index.Second_Cell_Index())){
                temporary_faces(i)=p_advected(face_index.Second_Cell_Index())*dt*face_area;}
            else temporary_faces(i)=T();}
        fluid_interpolation->Times(temporary_faces,V.lambda);

        // viscosity initial guess
        V.viscous_force_coefficients=F.viscous_force_coefficients;}

    if(solid_node) V.force_coefficients=F.force_coefficients;
}
//#####################################################################
// Function Print_Matrix
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Print_Matrix(const VECTOR_T& vec) const
{
    index_map.Print(solve_id);
    VECTOR_T V(vec),F(vec);
    std::stringstream ss;
    ss<<"Printing matrix: "<<std::endl;
    ss<<V.pressure.n<<" pressures"<<std::endl;
    ss<<V.lambda.Size()<<" lambdas"<<std::endl;
    ss<<V.force_coefficients.Size()<<" force coefficients"<<std::endl;
    ss<<V.viscous_force_coefficients.Size()<<" viscous force coefficients"<<std::endl;
    ss<<"Solve id: "<<solve_id<<std::endl;
    LOG::filecout(ss.str());
    OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("matrix-%i.txt",solve_id).c_str()).Write("M",*this,V,F);
    if(use_preconditioner) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("precond-%i.txt",solve_id).c_str()).Write_Preconditioner("P",*this,V,F);
}
//#####################################################################
// Function Print_Matrix
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Print_Each_Matrix(int n) const
{
     fluid_gradient->Print_Each_Matrix(n);
     if(!leakproof_solve){
         GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);
         solid_forces->Print_Each_Matrix(n,temporary_solids_velocity,sqrt(dt));
         solid_interpolation->Print_Each_Matrix(n,temporary_solids_velocity);
         if(fluid_to_solid_interpolation) fluid_to_solid_interpolation->Print_Each_Matrix(n,fluid_gradient->gradient.m,temporary_solids_velocity);

         OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("Mi-%i.txt",solve_id).c_str());
         oo.Begin_Sparse_Matrix("Mi",temporary_solids_velocity.Raw_Size(),temporary_solids_velocity.Raw_Size());

         for(int i=1;i<=solid_mass->one_over_mass.Size();i++)
             for(int j=1;j<=TV::m;j++)
                 oo.Append_Sparse_Diagonal_Block(solid_mass->one_over_mass(i));
         for(int i=1;i<=solid_mass->world_space_rigid_mass_inverse.Size();i++){
             for(int j=1;j<=TV::m;j++)
                 oo.Append_Sparse_Diagonal_Block(solid_mass->world_space_rigid_mass_inverse(i).mass);
             oo.Append_Sparse_Diagonal_Block(solid_mass->world_space_rigid_mass_inverse(i).inertia_tensor);}

         oo.End_Sparse_Matrix();}
     fluid_interpolation->Print_Each_Matrix(n);
     fluid_mass->Print_Each_Matrix(n);
     fluid_poisson->Print_Each_Matrix(n);
     if(use_viscous_forces) fluid_viscous_forces->Print_Each_Matrix(n);
     if(use_full_ic){
         OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("full_matrix-%i.txt",solve_id).c_str()).Write("FM",full_matrix);
         if(full_matrix.C) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("full_precon-%i.txt",solve_id).c_str()).Write("FP",*full_matrix.C);}
}
//#####################################################################
// Function Gather
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Gather(const KRYLOV_VECTOR_BASE<T>& bV,VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);
    if(fluid_to_solid_interpolation) return Massless_Gather(V,fluid_velocity,structure_velocity);
    if(!leakproof_solve && solid_node){
        solid_forces->Transpose_Times(V.force_coefficients,structure_velocity);
        structure_velocity*=sqrt(dt);}
    if(!leakproof_solve){
        if(mpi_solid_fluid) {
            if(fluid_node) solid_interpolation->Transpose_Times(V.lambda,structure_velocity);
            mpi_solid_fluid->Aggregate_Lists_To_Solid_Node(structure_velocity);}
        else if(solid_node) solid_interpolation->Transpose_Times_Add(V.lambda,structure_velocity);}
    if(fluid_node){
        fluid_interpolation->Transpose_Times(V.lambda,fluid_velocity);
        fluid_velocity=-fluid_velocity;
        fluid_gradient->Times_Add(V.pressure,fluid_velocity);
        if(use_viscous_forces) fluid_viscous_forces->Transpose_Times_Add(V.viscous_force_coefficients,fluid_velocity);}
}
//#####################################################################
// Function Scatter
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Scatter(const VECTOR_ND<T>& fluid_velocity,const GENERALIZED_VELOCITY<TV>& structure_velocity,KRYLOV_VECTOR_BASE<T>& bF) const
{
    VECTOR_T& F=debug_cast<VECTOR_T&>(bF);
    if(fluid_to_solid_interpolation) return Massless_Scatter(fluid_velocity,structure_velocity,F);
    if(fluid_node){
        fluid_interpolation->Times(fluid_velocity,F.lambda);
        F.lambda=-F.lambda;
        fluid_gradient->Transpose_Times(fluid_velocity,F.pressure);
        if(use_viscous_forces) fluid_viscous_forces->Times(fluid_velocity,F.viscous_force_coefficients);
        if(!leakproof_solve) solid_interpolation->Times_Add(structure_velocity,F.lambda);}
    if(solid_node && !leakproof_solve){
        solid_forces->Times(structure_velocity,F.force_coefficients);
        F.force_coefficients*=sqrt(dt);}
}
//#####################################################################
// Function Inverse_Mass
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Inverse_Mass(VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const
{
    if(!leakproof_solve && solid_node){
        solid_mass->Inverse_Multiply(structure_velocity,structure_velocity,true);
        solid_system->Project(structure_velocity);}
    if(fluid_node) fluid_mass->Inverse_Times(fluid_velocity,fluid_velocity);
}
//#####################################################################
// Function Massless_Gather
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Massless_Gather(const VECTOR_T& V,VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const
{
    GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);
    fluid_gradient->Times(V.pressure,fluid_velocity);
    if(use_viscous_forces) fluid_viscous_forces->Transpose_Times_Add(V.viscous_force_coefficients,fluid_velocity);
    solid_forces->Transpose_Times(V.force_coefficients,temporary_solids_velocity);
    temporary_solids_velocity*=sqrt(dt);
    fluid_to_solid_interpolation->Transpose_Times_Add(temporary_solids_velocity,fluid_velocity);
    structure_velocity*=(T)0;
}
//#####################################################################
// Function Massless_Scatter
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Massless_Scatter(const VECTOR_ND<T>& fluid_velocity,const GENERALIZED_VELOCITY<TV>& structure_velocity,VECTOR_T& F) const
{
    GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);
    fluid_gradient->Transpose_Times(fluid_velocity,F.pressure);
    if(use_viscous_forces) fluid_viscous_forces->Times(fluid_velocity,F.viscous_force_coefficients);
    fluid_to_solid_interpolation->Times(fluid_velocity,temporary_solids_velocity);
    temporary_solids_velocity*=sqrt(dt);
    solid_forces->Times(temporary_solids_velocity,F.force_coefficients);
    ARRAYS_COMPUTATIONS::Fill(F.lambda,(T)0);
}
//#####################################################################
// Function Apply_Massless_Structure_Force_To_Fluid
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Apply_Massless_Structure_Force_To_Fluid(VECTOR_ND<T>& fluid_velocity,T time) const
{
    if(!fluid_to_solid_interpolation) return;
    GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);
    VECTOR_ND<T> temp(temporary_faces.n);
    ARRAYS_COMPUTATIONS::Fill(temporary_solids_velocity.V.array,TV());
    ARRAYS_COMPUTATIONS::Fill(temporary_solids_velocity.rigid_V.array,TWIST<TV>());
    solid_forces->solid_body_collection.Add_Velocity_Independent_Forces(temporary_solids_velocity.V.array,temporary_solids_velocity.rigid_V.array,time);
    if(print_each_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("f-%i.txt",solve_id).c_str()).Write("f",temporary_solids_velocity);
    fluid_to_solid_interpolation->Transpose_Times(temporary_solids_velocity,temp);
    Dump_Substep(temp,temporary_solids_velocity,"after interp");
    fluid_mass->Inverse_Times(temp,temp);
    temp*=dt;
    fluid_velocity+=temp;
    Dump_Substep(fluid_velocity,"after add");
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);VECTOR_T& F=debug_cast<VECTOR_T&>(bF);
    GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);

    Gather(V,temporary_faces,temporary_solids_velocity);
    Inverse_Mass(temporary_faces,temporary_solids_velocity);

    Exchange_Velocities(temporary_faces,temporary_solids_velocity);
    Scatter(temporary_faces,temporary_solids_velocity,F);

    if(!leakproof_solve) F.force_coefficients+=V.force_coefficients;
    if(use_viscous_forces) F.viscous_force_coefficients+=V.viscous_force_coefficients;

    for(int i=1;i<=one_over_rho_c_squared_flat.n;i++) F.pressure(i)+=one_over_rho_c_squared_flat(i)*V.pressure(i);
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& bV) const
{
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Project() const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bR) const  // solve MR=V
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);VECTOR_T& R=debug_cast<VECTOR_T&>(bR);
    if(use_full_ic){
        for(int i=1;i<=full_precondition_in.n;i++)
            full_precondition_in(i)=const_cast<VECTOR_T&>(V).Raw_Get(i);

        full_matrix.C->Solve_Forward_Substitution(full_precondition_in,full_precondition_out,true);
        full_matrix.C->Solve_Backward_Substitution(full_precondition_out,full_precondition_in,false,true);

        for(int i=1;i<=full_precondition_in.n;i++)
            R.Raw_Get(i)=full_precondition_in(i);}
    else if(fluid_node){
        R.Copy(1,V);
        if(fluid_poisson->dt==0) Exchange_Pressure(R.pressure);
        fluid_poisson->Apply_Preconditioner(R.pressure);

        for(COUPLING_CONSTRAINT_ID i(1);i<=lambda_diagonal_preconditioner.m;i++)
            R.lambda(i)*=lambda_diagonal_preconditioner(i);}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& bV1,const KRYLOV_VECTOR_BASE<T>& bV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(bV1);const VECTOR_T& V2=debug_cast<const VECTOR_T&>(bV2);
    double inner_product_pressure=Dot_Product_Double_Precision(V1.pressure,V2.pressure);
    double inner_product_lambda=ARRAYS_COMPUTATIONS::Dot_Product_Double_Precision(V1.lambda,V2.lambda);
    double inner_product_force_coefficients=0;
    if(!leakproof_solve) inner_product_force_coefficients=ARRAYS_COMPUTATIONS::Dot_Product_Double_Precision(V1.force_coefficients,V2.force_coefficients);
    double inner_product_viscous_force_coefficients=ARRAYS_COMPUTATIONS::Dot_Product_Double_Precision(V1.viscous_force_coefficients,V2.viscous_force_coefficients);
    //LOG::cout<<"Inner products: p: "<<inner_product_pressure<<" lambda: "<<inner_product_lambda<<" solid: "<<inner_product_force_coefficients<<std::endl;
    double ret=inner_product_pressure+inner_product_lambda+inner_product_force_coefficients+inner_product_viscous_force_coefficients;
    T tmp=0;
    if(mpi_solid_fluid){mpi_solid_fluid->Reduce_Add((T)ret,tmp);return tmp;}
    else return ret;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(bR);
    bool pressure_converged=true;
    for(int i=1;i<=index_map.real_cell_indices.m;++i){int index=index_map.real_cell_indices(i);
        if(abs(R.pressure(index))>tolerances.pressure(index)) pressure_converged=false;}
    bool lambda_converged=true;
    for(COUPLING_CONSTRAINT_ID i(1);i<=tolerances.lambda.Size();i++)
        if(abs(R.lambda(i))>tolerances.lambda(i)) lambda_converged=false;
    bool forces_converged=true;
    for(FORCE_AGGREGATE_ID i(1);i<=tolerances.force_coefficients.Size();i++)
        if(abs(R.force_coefficients(i))>tolerances.force_coefficients(i)) forces_converged=false;
    bool viscous_forces_converged=true;
    for(VISCOUS_FORCE_ID i(1);i<=tolerances.viscous_force_coefficients.Size();i++)
        if(abs(R.viscous_force_coefficients(i))>tolerances.viscous_force_coefficients(i)) viscous_forces_converged=false;
    T ret=(pressure_converged && lambda_converged && forces_converged && viscous_forces_converged)?0:FLT_MAX;
    if(mpi_solid_fluid){ret=mpi_solid_fluid->Reduce_Max(ret);return ret;}
    else return ret;
}
//#####################################################################
// Function Linf_Norm
//#####################################################################
template<class TV> typename TV::SCALAR SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Linf_Norm(const VECTOR_T& R) const
{
    T fluid_convergence_norm=(T)0;
    TV_INT max_pressure_error_cell=TV_INT();
    for(int i=1;i<=index_map.real_cell_indices.m;i++){
        T pressure_error=abs(R.pressure(index_map.real_cell_indices(i)));
        if(pressure_error>fluid_convergence_norm){
            max_pressure_error_cell=index_map.indexed_cells(index_map.real_cell_indices(i));
            fluid_convergence_norm=pressure_error;}}

    T coupling_convergence_norm=0;
    COUPLING_CONSTRAINT_ID max_lambda_error=COUPLING_CONSTRAINT_ID(0);
    for(COUPLING_CONSTRAINT_ID i(1);i<=R.lambda.Size();i++){
        T lambda_error=abs(R.lambda(i));
        if(lambda_error>coupling_convergence_norm){
            max_lambda_error=i;
            coupling_convergence_norm=lambda_error;}}

    T forces_convergence_norm=(T)0;
    for(FORCE_AGGREGATE_ID i(1);i<=R.force_coefficients.Size();i++) forces_convergence_norm=max(forces_convergence_norm,abs(R.force_coefficients(i)));

    T viscous_forces_convergence_norm=(T)0;
    for(VISCOUS_FORCE_ID i(1);i<=R.viscous_force_coefficients.Size();i++) viscous_forces_convergence_norm=max(viscous_forces_convergence_norm,abs(R.viscous_force_coefficients(i)));
    T convergence_norm=max(fluid_convergence_norm,coupling_convergence_norm,forces_convergence_norm,viscous_forces_convergence_norm);

    std::stringstream ss;
    ss<<"Max pressure error at cell "<<max_pressure_error_cell<<" = "<<fluid_convergence_norm<<std::endl;
    ss<<"Max lambda error at constraint "<<max_lambda_error<<" = "<<coupling_convergence_norm<<std::endl;
    ss<<"Max forces error = "<<forces_convergence_norm<<std::endl;
    ss<<"Max viscous forces error = "<<viscous_forces_convergence_norm<<std::endl;
    LOG::filecout(ss.str());

    return mpi_solid_fluid?mpi_solid_fluid->Reduce_Max(convergence_norm):convergence_norm;
}
//#####################################################################
// Function Residual_Linf_Norm
//#####################################################################
template<class TV> typename TV::SCALAR SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Residual_Linf_Norm(const VECTOR_T& x,const VECTOR_T& rhs) const
{
    VECTOR_T residual;
    Resize_Coupled_System_Vector(residual);
    Multiply(x,residual);
    residual-=rhs;
    return Linf_Norm(residual);
}
//#####################################################################
// Function Apply
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Apply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Test_Matrix() const
{
    fluid_gradient->Test_Matrix();
    fluid_poisson->Test_Matrix(print_poisson_matrix);
    fluid_interpolation->Test_Matrix();
    if(!leakproof_solve && solid_node) solid_forces->Test_Matrix();
    int number=solid_system->solid_body_collection.deformable_body_collection.deformable_geometry.particles.array_collection->number;
    int rigid_number=solid_system->solid_body_collection.rigid_body_collection.rigid_geometry_collection.particles.array_collection->number;
    solid_interpolation->Test_Matrix(number,rigid_number);
    if(use_viscous_forces) fluid_viscous_forces->Test_Matrix();

    VECTOR_T a;
    Resize_Coupled_System_Vector(a);
    VECTOR_T b(a),c(a);
    Test_System(a,b,c);
}
//#####################################################################
// Function Test_Incompressibility
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Test_Incompressibility(const ARRAY<T,FACE_INDEX<TV::dimension> >& fluid_velocity,const ARRAY<T>& constrained_fluids_velocity) const
{
    VECTOR_T a;
    Resize_Coupled_System_Vector(a);
    if(!a.pressure.n) return;
    VECTOR_ND<T> fluid_velocity_vector;
    index_map.Collect(fluid_velocity,constrained_fluids_velocity,fluid_velocity_vector);
    fluid_gradient->Transpose_Times(fluid_velocity_vector,a.pressure);
    T mn=a.pressure.Min(),mx=a.pressure.Max(),vol=index_map.grid.Cell_Size();
    std::stringstream ss;
    ss<<"divergence: ["<<mn/vol<<" "<<mx/vol<<"]     vol weighted: ["<<mn<<" "<<mx<<"]"<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Set_Coupling_Faces
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Set_Coupling_Faces(const int ghost_cells,T_FACE_ARRAYS_BOOL& psi_N) const
{
    typename GRID<TV>::REGION region_type=ghost_cells?GRID<TV>::INTERIOR_REGION:GRID<TV>::WHOLE_REGION;
    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(index_map.iterator_info,ghost_cells,region_type);iterator.Valid();iterator.Next())
        psi_N(iterator.axis,iterator.index)=true;
}
//#####################################################################
// Function Show_Constraints
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Show_Constraints(T_FACE_ARRAYS_BOOL& psi_N) const
{
    T_FACE_ARRAYS_BOOL store_psi_N=psi_N;
    psi_N.Fill(false);
    int ghost_cells=0;
    typename GRID<TV>::REGION region_type=ghost_cells?GRID<TV>::INTERIOR_REGION:GRID<TV>::WHOLE_REGION;
    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(index_map.iterator_info,ghost_cells,region_type);iterator.Valid();iterator.Next())
        psi_N(iterator.axis,iterator.index)=true;
    psi_N=store_psi_N;
}
//#####################################################################
// Function Check_Constraints
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Check_Constraints(const GENERALIZED_VELOCITY<TV>& solids_rhs,const ARRAY<T,FACE_INDEX<TV::dimension> >& fluids_rhs,const ARRAY<T>& constrained_rhs) const
{
    index_map.Collect(fluids_rhs,constrained_rhs,temporary_faces);
    fluid_interpolation->Times(temporary_faces,coupling_faces);
    coupling_faces=-coupling_faces;
    solid_interpolation->Times_Add(solids_rhs,coupling_faces);
    std::stringstream ss;
    ss<<coupling_faces<<std::endl;

    VECTOR_ND<T> fluid_velocity_vector;
    index_map.Collect(fluids_rhs,constrained_rhs,fluid_velocity_vector);
    fluid_gradient->Transpose_Times(fluid_velocity_vector,pressure);
    ss<<"Incompressibility constraints"<<std::endl;
    ss<<pressure/index_map.grid.Cell_Size()<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Compute_Lambda_Diagonal_Preconditioner
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Compute_Lambda_Diagonal_Preconditioner()
{
    if(!fluid_node) return;
    lambda_diagonal_preconditioner.Resize(solid_interpolation->Number_Of_Constraints());
    ARRAYS_COMPUTATIONS::Fill(lambda_diagonal_preconditioner,(T)0);
    solid_interpolation->Add_Diagonal(lambda_diagonal_preconditioner,*solid_mass);
    fluid_interpolation->Add_Diagonal(lambda_diagonal_preconditioner,*fluid_mass);
    for(COUPLING_CONSTRAINT_ID i(1);i<=lambda_diagonal_preconditioner.m;i++)
        lambda_diagonal_preconditioner(i)=lambda_diagonal_preconditioner(i)?1/lambda_diagonal_preconditioner(i):1;
}
//#####################################################################
// Function Compute_Scatter_Matrix
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Compute_Scatter_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& gather_matrix)
{
    SYSTEM_MATRIX_HELPER<T> matrix_helper;
    GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);
    
    int size_p=index_map.Number_Cells(),size_lambda=Value(fluid_interpolation->Number_Of_Constraints()),size_force=leakproof_solve?0:Value(solid_forces->Velocity_Dependent_Forces_Size());
    int size_viscous=use_viscous_forces?Value(fluid_viscous_forces->Viscous_Forces_Size()):0;
    int offset_lambda=size_p,offset_force=offset_lambda+size_lambda,offset_viscous=offset_force+size_force;
    int size_fluids=index_map.Number_Faces(),size_solids_actual=temporary_solids_velocity.Raw_Size(),size_solids=fluid_to_solid_interpolation?0:size_solids_actual;
    solid_interpolation->Store_Maps(temporary_solids_velocity);

    if(fluid_to_solid_interpolation){
        fluid_to_solid_interpolation->Store_Maps(temporary_solids_velocity);
        matrix_helper.Add_Matrix(*fluid_gradient,true);
        SPARSE_MATRIX_FLAT_MXN<T> C,H,CH;
        SYSTEM_MATRIX_HELPER<T>::Base_To_Matrix(size_force,size_solids_actual,*solid_forces,C);
        SYSTEM_MATRIX_HELPER<T>::Base_To_Matrix(size_solids_actual,size_fluids,*fluid_to_solid_interpolation,H);
        CH=C*H;
        CH*=sqrt(dt);
        matrix_helper.Add_Matrix(CH,false,offset_force,0);
        if(use_viscous_forces) matrix_helper.Add_Matrix(*fluid_viscous_forces,false,offset_viscous,0);}
    else{
        if(fluid_node){
            matrix_helper.Add_Matrix(*fluid_gradient,true);
            matrix_helper.Add_Matrix(*fluid_interpolation,false,offset_lambda,0);
            matrix_helper.Scale(-1);
            if(use_viscous_forces) matrix_helper.Add_Matrix(*fluid_viscous_forces,false,offset_viscous,0);
            if(!leakproof_solve) matrix_helper.Add_Matrix(*solid_interpolation,false,offset_lambda,size_fluids);}
        if(solid_node && !leakproof_solve){
            matrix_helper.Add_Matrix(*solid_forces,false,offset_force,size_fluids);
            matrix_helper.Scale(sqrt(dt));}}

    matrix_helper.Compact();
    matrix_helper.Set_Matrix(offset_viscous+size_viscous,size_fluids+size_solids,gather_matrix);

    if(print_each_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("SC-%i.txt",solve_id).c_str()).Write("SC",gather_matrix);
}
//#####################################################################
// Function Compute_Inverse_Mass_Matrix
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Compute_Inverse_Mass_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& inverse_mass)
{
    GENERALIZED_VELOCITY<TV> temporary_solids_velocity(temporary_velocities,temporary_twists,solid_system->solid_body_collection);
    SYSTEM_MATRIX_HELPER<T> matrix_helper;
    int size_fluids=index_map.Number_Faces(),size_solids=fluid_to_solid_interpolation?0:temporary_solids_velocity.Raw_Size();

    if(!leakproof_solve && solid_node && !fluid_to_solid_interpolation){
        matrix_helper.New_Block();
        int k=1;
        for(int i=1;i<=solid_mass->one_over_mass.Size();i++)
            for(int j=1;j<=TV::m;j++,k++)
                matrix_helper.data.Append(TRIPLE<int,int,T>(k,k,solid_mass->one_over_mass(i)));
        for(int i=1;i<=solid_mass->world_space_rigid_mass_inverse.Size();i++){
            for(int j=1;j<=TV::m;j++,k++)
                matrix_helper.data.Append(TRIPLE<int,int,T>(k,k,solid_mass->world_space_rigid_mass_inverse(i).mass));
            for(int j=1;j<=TV::SPIN::m;j++)
                for(int n=1;n<=TV::SPIN::m;n++)
                    matrix_helper.data.Append(TRIPLE<int,int,T>(k+j-1,k+n-1,solid_mass->world_space_rigid_mass_inverse(i).inertia_tensor(j,n)));
            k+=TV::SPIN::m;}
        matrix_helper.Shift(size_fluids,size_fluids);
        /*solid_system->Project(structure_velocity);*/}

    if(fluid_node) matrix_helper.Add_Matrix(*fluid_mass);

    {int mx=0,my=0;for(int i=1;i<=matrix_helper.data.m;i++) mx=std::max(mx,matrix_helper.data(i).x),my=std::max(my,matrix_helper.data(i).y);
    std::stringstream ss;ss<<"max: "<<mx<<"  "<<my<<"  vs "<<size_fluids+size_solids<<"  "<<size_fluids+size_solids<<std::endl;LOG::filecout(ss.str());}

    matrix_helper.Compact();
    matrix_helper.Set_Matrix(size_fluids+size_solids,size_fluids+size_solids,inverse_mass);
}
//#####################################################################
// Function Compute_Full_Preconditioner
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Compute_Full_Preconditioner()
{
    if(mpi_solid_fluid) PHYSBAM_FATAL_ERROR("Not implemented for MPI two-way coupled sims.");
    SPARSE_MATRIX_FLAT_MXN<T> scatter_matrix,inverse_mass,sc_im;
    Compute_Scatter_Matrix(scatter_matrix);
    Compute_Inverse_Mass_Matrix(inverse_mass);

    sc_im=scatter_matrix*inverse_mass;
    full_matrix=sc_im.Times_Transpose(scatter_matrix).Create_NXN_Matrix();

    int size_p=index_map.Number_Cells(),size_lambda=Value(fluid_interpolation->Number_Of_Constraints()),size_force=leakproof_solve?0:Value(solid_forces->Velocity_Dependent_Forces_Size());
    int size_viscous=use_viscous_forces?Value(fluid_viscous_forces->Viscous_Forces_Size()):0;
    int offset_lambda=size_p,offset_force=offset_lambda+size_lambda,offset_viscous=offset_force+size_force;

    if(!leakproof_solve) for(int i=1;i<=size_force;i++) full_matrix(offset_force+i,offset_force+i)+=1;
    if(use_viscous_forces) for(int i=1;i<=size_viscous;i++) full_matrix(offset_viscous+i,offset_viscous+i)+=1;
    for(int i=1;i<=one_over_rho_c_squared_flat.n;i++) full_matrix(i,i)+=one_over_rho_c_squared_flat(i);

    full_matrix.Construct_Incomplete_Cholesky_Factorization();
    full_precondition_in.Resize(full_matrix.n);
    full_precondition_out.Resize(full_matrix.n);
}
//#####################################################################
// Function Exchange_Pressure
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Exchange_Pressure(VECTOR_ND<T>& pressure) const
{
    // TODO: This is very inefficient.
    if(mpi_grid){
        PHYSBAM_FATAL_ERROR("Implement this fully before using.");
        index_map.Distribute(pressure,pressure_on_grid);
        mpi_grid->Exchange_Boundary_Cell_Data(pressure_on_grid,1,true);
        index_map.Collect_Indexed_Cells(pressure_on_grid,pressure);}
}
//#####################################################################
// Function Exchange_Velocities
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Exchange_Velocities(VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& structure_velocity) const
{
    // TODO: Exchange enough information to do:
    //    fluid node: fluid_viscous_forces->Times
    if(mpi_solid_fluid) mpi_solid_fluid->Distribute_Lists_From_Solid_Node(structure_velocity);
    if(fluid_node && mpi_grid){
        fluid_velocity_on_grid.Resize(index_map.grid.Domain_Indices(0));fluid_velocity_on_grid.Fill(0);
        index_map.Distribute_Boundary_Faces(fluid_velocity,fluid_velocity_on_grid);
        mpi_grid->Sum_Common_Face_Data(fluid_velocity_on_grid);
        index_map.Collect_Boundary_Faces(fluid_velocity_on_grid,fluid_velocity);}
}
//#####################################################################
// Function Set_MPI
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Set_MPI(MPI_SOLID_FLUID<TV>& mpi_solid_fluid_input,MPI_UNIFORM_GRID<GRID<TV> >& mpi_grid_input)
{
    mpi_grid=&mpi_grid_input;
    mpi_solid_fluid=&mpi_solid_fluid_input;
    solid_node=mpi_solid_fluid->Solid_Node();
    fluid_node=mpi_solid_fluid->Fluid_Node();
}
//#####################################################################
// Function Dump_Substep
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Dump_Substep(const VECTOR_ND<T>& fluid_velocity,const char* name,int substep,int level) const
{
    ARRAY<T,FACE_INDEX<TV::m> > tmp(debug_velocity->Domain_Indices());
    ARRAY<T> tmp_constrained(index_map.indexed_constraints.m);
    index_map.Distribute(fluid_velocity,tmp,tmp_constrained);
    debug_velocity->Exchange_Arrays(*debug_velocity,tmp);
    if(fluid_to_solid_interpolation)
        if(FLUID_TO_SOLID_INTERPOLATION_CUT<TV>* cut=dynamic_cast<FLUID_TO_SOLID_INTERPOLATION_CUT<TV>*>(fluid_to_solid_interpolation))
            Dump_Extra_Velocities<TV>(*cut,fluid_velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(name,substep,level);
    debug_velocity->Exchange_Arrays(*debug_velocity,tmp);
}
//#####################################################################
// Function Dump_Substep
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Dump_Substep(const VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity,const char* name,int substep,int level) const
{
    debug_generalized_velocity->Exchange(solid_velocity);
    Dump_Substep(fluid_velocity,name,substep,level);
    debug_generalized_velocity->Exchange(solid_velocity);
}
//#####################################################################
// Function Fill_Extra_Velocities
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Fill_Extra_Velocities(VECTOR_ND<T>& fluid_velocity_vector) const
{
    if(FLUID_TO_SOLID_INTERPOLATION_CUT<TV>* interp=dynamic_cast<FLUID_TO_SOLID_INTERPOLATION_CUT<TV>*>(fluid_to_solid_interpolation))
        interp->Fill_Extra_Velocities(fluid_velocity_vector);
    Dump_Substep(fluid_velocity_vector,"fill extra velocities");
}
template<> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<float,1> >::Fill_Extra_Velocities(VECTOR_ND<T>& fluid_velocity_vector) const {PHYSBAM_FATAL_ERROR();}
template<> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<double,1> >::Fill_Extra_Velocities(VECTOR_ND<T>& fluid_velocity_vector) const {PHYSBAM_FATAL_ERROR();}
template<> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<float,3> >::Fill_Extra_Velocities(VECTOR_ND<T>& fluid_velocity_vector) const {PHYSBAM_FATAL_ERROR();}
template<> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<double,3> >::Fill_Extra_Velocities(VECTOR_ND<T>& fluid_velocity_vector) const {PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Mark_Valid_Faces
//#####################################################################
template<class TV> void SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::
Mark_Valid_Faces(ARRAY<bool,FACE_INDEX<TV::m> >& valid) const
{
    const ARRAY<int>& offsets=fluid_gradient->gradient.offsets;
    valid.Fill(false);
    for(int i=1;i<=index_map.indexed_faces.m-1;i++){
        if(offsets(i)<offsets(i+1)){
            valid(index_map.indexed_faces(i))=true;
            /*Add_Debug_Particle(index_map.grid.Axis_X_Face(index_map.indexed_faces(i)),VECTOR<typename TV::SCALAR,3>(0,1,0));*/}}
    Dump_Substep(temporary_faces,"after mark valid faces");
}
template<class TV> int SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>::solve_id=0;
//#####################################################################
template class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<float,1> >;
template class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<float,2> >;
template class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<double,1> >;
template class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<double,2> >;
template class SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<VECTOR<double,3> >;
#endif
}
