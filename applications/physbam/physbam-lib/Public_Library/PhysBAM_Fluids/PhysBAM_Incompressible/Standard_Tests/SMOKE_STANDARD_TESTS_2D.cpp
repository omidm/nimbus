//#####################################################################
// Copyright 2006-2007, Jon Gretarsson, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_INTERSECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;
template<class T_GRID> SMOKE_STANDARD_TESTS_2D<T_GRID>::
SMOKE_STANDARD_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container,
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :example(example),fluids_parameters(fluids_parameters),incompressible_fluid_container(incompressible_fluid_container),rigid_body_collection(rigid_body_collection)
{
    
    // set up the standard fluid environment
    // TODO: *REALLY* need to pick sensible constants and settings
    example.frame_rate=24;
    fluids_parameters.cfl=(T).9;
    fluids_parameters.domain_walls=VECTOR<VECTOR<bool,2>,T_GRID::dimension>::Constant_Vector(VECTOR<bool,2>::Constant_Vector(false));
    fluids_parameters.domain_walls(2)(1)=true;
    fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).04;        
    fluids_parameters.kolmogorov=(T)0;fluids_parameters.gravity=0;
    rho=1;fluids_parameters.rho_bottom=1;fluids_parameters.rho_top=(T).65;
    fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=0;
    fluids_parameters.temperature_container.Set_Cooling_Constant(0);fluids_parameters.temperature_products=(T)3000;
    fluids_parameters.write_velocity=true;example.write_frame_title=true;fluids_parameters.write_debug_data=true;
    fluids_parameters.use_body_force=true;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T_GRID> void SMOKE_STANDARD_TESTS_2D<T_GRID>::
Initialize(const int test_number_input,const int resolution,const T angle_fraction)
{
    test_number=test_number_input;
    std::stringstream ss;ss<<"Running Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;LOG::filecout(ss.str());

    // set up the domain
    int cells=resolution;
    if(test_number==1 || test_number==2 || test_number==3){
        example.first_frame=0;example.last_frame=3840;
        grid.Initialize(10*cells+1,15*cells+1,0,1,0,1.5);}
    else if(test_number==4){
        fluids_parameters.domain_walls=VECTOR<VECTOR<bool,2>,T_GRID::dimension>::Constant_Vector(VECTOR<bool,2>::Constant_Vector(true));
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.gravity=(T)0;
        fluids_parameters.use_poisson=true;  // TODO This is a hack to tell Projection to use poisson rather than laplace
        // TODO: need to call use_variable_beta on the poisson solver if we ACTUALLY want to use variable beta
        fluids_parameters.second_order_cut_cell_method=false;
        example.first_frame=0;example.last_frame=1000;
        grid.Initialize(10*cells+1,10*cells+1,(T)-1,(T)1,(T)-1,(T)1);}
    else{LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

    example.output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests_Smoke/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));

    // set up the source domain
    if(test_number==1 || test_number==2 || test_number==3){
        source=BOX<TV>((T).45,(T).55,(T)0,(T).1);
        world_to_source=MATRIX<T,3>::Identity_Matrix();
        source_velocity=VECTOR<T,2>((T)0,(T)0.5);}

    //set up example-specific parameters
    if(test_number==3){explosion_divergence=20;explosion_end_time=3;fluids_parameters.use_non_zero_divergence=true;}

    rotation_angle=angle_fraction?(T)pi/angle_fraction:(T)0;
    rotation_frame=FRAME<TV>(ROTATION<TV>::From_Angle(rotation_angle));
}
template<class T_GRID> SMOKE_STANDARD_TESTS_2D<T_GRID>::
~SMOKE_STANDARD_TESTS_2D()
{
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,2> SMOKE_STANDARD_TESTS_2D<T_GRID>::
Initial_Velocity(const VECTOR<T,2>& X) const
{
    if(test_number==4) return TV(-10*cos(X.x)*sin(X.y),10*sin(X.x)*cos(X.y));
    return VECTOR<T,2>();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_GRID> void SMOKE_STANDARD_TESTS_2D<T_GRID>::
Initialize_Bodies()
{
    if(test_number==2){
        int id=rigid_body_collection.Add_Rigid_Body(example.stream_type,example.data_directory+"/Rigid_Bodies_2D/circle",(T).125,true,true,false);
        rigid_body_collection.rigid_body_particle.X(id)=VECTOR<T,2>((T).5,(T).75);
        rigid_body_collection.Rigid_Body(id).is_static=true;
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);}
    else if(test_number==4){
        oriented_box=rigid_body_collection.Add_Rigid_Body(example.stream_type,example.data_directory+"/Rigid_Bodies_2D/square",(T).8/sqrt((T)2),true,true,false);
        rigid_body_collection.rigid_body_particle.rotation(oriented_box)=ROTATION<TV>::From_Angle(rotation_angle);
        rigid_body_collection.Rigid_Body(oriented_box).is_static=true;

        // Calculate the volume fraction for our new divergence face weights (see Normal stuff)
        IMPLICIT_OBJECT_INTERSECTOR<TV> implicit_object(*rigid_body_collection.Rigid_Body(oriented_box).implicit_object);
        T_GRID& grid=*fluids_parameters.grid;
        divergence_face_weights.Resize(grid,1);beta_face.Resize(grid,1);
        T one_over_cell_size=(T)1/grid.Cell_Size();

        // TODO This appears to be broken for non-zero angles
        for(typename T_GRID::FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();VECTOR<int,2> face_index=iterator.Face_Index();RANGE<TV> dual_cell=iterator.Dual_Cell();
            T normal_weight=(T)1-sqr(rigid_body_collection.Rigid_Body(oriented_box).implicit_object->Extended_Normal(iterator.Location())(axis));
            T volume_fraction=implicit_object.Negative_Material_In_Box(dual_cell)*one_over_cell_size;
            if(volume_fraction<=1e-3) beta_face(axis,face_index)=0;
            else beta_face(axis,face_index)=volume_fraction;
            if(volume_fraction<=1-1e-3){
                if(normal_weight<=1e-3) divergence_face_weights(axis,face_index)=0;
                else divergence_face_weights(axis,face_index)=normal_weight;}
            else 
                divergence_face_weights(axis,face_index)=1;
            
            if(volume_fraction<=1e-3) beta_face(axis,face_index)=0;
            else beta_face(axis,face_index)=volume_fraction;}
    }

    PHYSBAM_DEBUG_WRITE_SUBSTEP("After calling Initialize Bodies",0,0);
}
//#####################################################################
// Function Get_Divergence
//#####################################################################
template<class T_GRID> void SMOKE_STANDARD_TESTS_2D<T_GRID>::
Get_Divergence(ARRAY<T,VECTOR<int,2> >& divergence,const T dt,const T time)
{
    if(test_number==3){
        T expansion=explosion_divergence*sin(time)/exp(time);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) 
            if(source.Lazy_Inside(iterator.Location())) divergence(iterator.Cell_Index())=expansion;}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class T_GRID> void SMOKE_STANDARD_TESTS_2D<T_GRID>::
Get_Object_Velocities(PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection,const T dt,const T time)
{
    if(test_number==4){
         T_GRID& grid=*fluids_parameters.grid;
         projection.poisson->Set_Variable_beta(true);
         T_FACE_ARRAYS_SCALAR::Copy(beta_face,projection.poisson->beta_face);
         projection.poisson->Use_Weighted_Divergence();
         T_FACE_ARRAYS_SCALAR::Copy(divergence_face_weights,projection.poisson->divergence_face_weights);

         // We're changing how to calculate Neumann Boundary Conditions
         for(typename T_GRID::FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
             int axis=iterator.Axis();VECTOR<int,2> face_index=iterator.Face_Index();
             if(divergence_face_weights(axis,face_index)<=1e-3 || beta_face(axis,face_index)<=1e-3){
                 projection.elliptic_solver->psi_N(axis,face_index)=true;incompressible_fluid_container.face_velocities(axis,face_index)=0;}}}
}
template class SMOKE_STANDARD_TESTS_2D<GRID<VECTOR<float,2> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SMOKE_STANDARD_TESTS_2D<GRID<VECTOR<double,2> > >;
#endif
