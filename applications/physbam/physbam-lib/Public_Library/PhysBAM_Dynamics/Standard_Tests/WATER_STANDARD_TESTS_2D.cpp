//#####################################################################
// Copyright 2005-2008, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_2D
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> WATER_STANDARD_TESTS_2D<T_GRID>::
WATER_STANDARD_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS<T_GRID>& fluids_parameters,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :example(example),fluids_parameters(fluids_parameters),rigid_body_collection(rigid_body_collection_input),inaccurate_union(*fluids_parameters.grid),
    use_inaccurate_body_collisions(false),use_variable_density_for_sph(false),use_two_way_coupling_for_sph(false),convert_sph_particles_to_fluid(false),use_analytic_divergence(false),
    use_analytic_divergence_for_expansion_only(false),adjust_cell_weights_on_neumann_boundaries(false),enforce_density_near_interface(true),flip_ratio(1),target_particles_per_unit_volume(1),
    neumann_boundary_slip_multiplier(0),ballistic_particles_as_percentage_of_target((T).03),particle_targeting_time(1),test_number(0),sphere(0)
{
}
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Initialize(const int test_number_input,const int resolution)
{
    test_number=test_number_input;
    std::stringstream ss;ss<<"Running Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;LOG::filecout(ss.str());

    // set up the standard fluid environment
    example.frame_rate=24;
    example.restart=false;example.restart_frame=0;
    fluids_parameters.domain_walls(1)(1)=true;fluids_parameters.domain_walls(1)(2)=true;fluids_parameters.domain_walls(2)(1)=true;fluids_parameters.domain_walls(2)(2)=false;
    fluids_parameters.number_particles_per_cell=16;
    fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
    fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_debug_data=true;
    fluids_parameters.write_ghost_values=true;fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.delete_fluid_inside_objects=true;
    fluids_parameters.incompressible_iterations=40;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.store_particle_ids=true;
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_vorticity_confinement_fuel=false;

    // set up the domain
    int cells=1*resolution;
    if(test_number==1){
        example.first_frame=0;example.last_frame=10;
        grid.Initialize(15*cells+1,10*cells+1,0,(T)1.5,0,1);}
    else if(test_number==2){
        example.first_frame=0;example.last_frame=100;
        grid.Initialize(10*cells+1,15*cells+1,0,1,0,(T)1.5);}
    else if(test_number==3){
        example.first_frame=0;example.last_frame=150;
        grid.Initialize(10*cells+1,15*cells+1,0,1,0,(T)1.5);}
    else if(test_number==4){
        example.first_frame=0;example.last_frame=200;
        grid.Initialize(15*cells+1,10*cells+1,0,(T)1.5,0,1);
        motion_curve.Add_Control_Point(0,TV((T)1.25,(T).55));
        motion_curve.Add_Control_Point((T).2,TV((T).8,(T).1)); // .03 was old
        motion_curve.Add_Control_Point(3,TV((T).8,(T).1));}
    else if(test_number==5){
        example.first_frame=0;example.last_frame=150;
        grid.Initialize(10*cells+1,25*cells+1,0,(T).1,0,(T).25);}
    else if(test_number==6){
        example.first_frame=0;example.last_frame=50;
        grid.Initialize(25*cells+1,10*cells+1,0,1,0,(T).4);}
    else if(test_number==7){
        fluids_parameters.gravity=(T)0;
        example.first_frame=0;example.last_frame=100;
        grid.Initialize(10*cells+1,15*cells+1,0,1,0,(T)1.5);}
    else if(test_number==8){
        example.first_frame=0;example.last_frame=1000;
        grid.Initialize(10*cells+1,15*cells+1,0,1,0,(T)1.5);}
    else if(test_number==9){
        example.first_frame=0;example.last_frame=1000;
        grid.Initialize(16*cells+1,6*cells+1,0,4,0,(T)1.5);}
    else if(test_number==10){
        example.first_frame=0;example.last_frame=1000;
        grid.Initialize(6*cells+1,12*cells+1,0,(T)1.5,0,3);}
    else if(test_number==11){
        fluids_parameters.gravity=(T)0;
        example.first_frame=0;example.last_frame=100;
        grid.Initialize(4*cells+1,6*cells+1,0,1,0,(T)1.5);}
    else if(test_number==12){
        fluids_parameters.viscosity=(T)1000;fluids_parameters.implicit_viscosity=true;
        example.first_frame=0;example.last_frame=100;
        grid.Initialize(10*cells+1,15*cells+1,0,1,0,(T)1.5);}
    else if(test_number==13){
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=true;
        fluids_parameters.use_explicit_part_of_implicit_viscosity=true;
        example.first_frame=0;example.last_frame=100;
        grid.Initialize(10*cells,15*cells,0,1,0,(T)1.5,true);}
    else if(test_number==14){
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=true;
        fluids_parameters.use_explicit_part_of_implicit_viscosity=true;
        example.first_frame=0;example.last_frame=100;
        grid.Initialize(10*cells,15*cells,0,1,0,(T)1.5,true);}
    else if(test_number==15){
        fluids_parameters.analytic_test=true;
        //fluids_parameters.viscosity=(T)1000;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=false;
        //fluids_parameters.use_explicit_part_of_implicit_viscosity=false;
        example.first_frame=0;example.last_frame=100;
        grid.Initialize(2*cells,3*cells,0,1,0,(T)1.5,true);}
    else if(test_number==20){ //MIKE example
        example.first_frame=0;example.last_frame=2000;
        //grid.Initialize(15*cells+1,10*cells+1,10*cells+1,0,(T)1.5,0,1,0,1);
        grid.Initialize(201,81,0,(T)2,0,(T).8);
        fluids_parameters.cfl=(T)1.9;
        fluids_parameters.incompressible_iterations=100;
        //fluids_parameters.cg_restart_iterations=20;
        //fluids_parameters.evolution_solver_type=krylov_solver_cr;
        motion_curve.Add_Control_Point(0,TV((T)1,(T).1));
        motion_curve.Add_Control_Point((T).075,TV((T)1,(T).1));
        motion_curve.Add_Control_Point(3,TV((T)1,(T).55));}
    else if(test_number==21){
        example.first_frame=0;example.last_frame=1000;
        fluids_parameters.cfl=(T)0.9;
        fluids_parameters.incompressible_iterations=100;
        grid.Initialize(201,81,0,(T)1,0,(T).8);}
    else{
        LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

    example.output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));

    // set up sources for each test case
    if(test_number==3){
        TV domain_center=grid.domain.Center();domain_center.y=(T)1;
        sources.Resize(2);source_velocity.Resize(2);world_to_source.Resize(2);
        sources(1)=sources(2)=BOX<TV>((T).10787,(T).2714532,-(T).10787,(T).10787);
        ARRAY<MATRIX<T,3> > source_to_world(2);
        MATRIX<T,3> rotation=MATRIX<T,3>::Rotation_Matrix_Z_Axis((T)pi);
        MATRIX<T,3> translation=MATRIX<T,3>::Translation_Matrix(domain_center);
        source_to_world(1)=translation;
        source_to_world(2)=translation*rotation;
        for(int i=1;i<=2;i++){source_velocity(i)=source_to_world(i).Extract_Rotation()*TV((T).5,0);world_to_source(i)=source_to_world(i).Inverse();}}
    if(test_number==5){
        world_to_source.Append(MATRIX<T,3>::Rotation_Matrix_Z_Axis((T)pi/(T)4)*MATRIX<T,3>::Translation_Matrix(TV((T)-.03,(T)-.23)));
        sources.Append(BOX<TV>(-(T).010,(T).010,-(T).010,(T).010));
        source_velocity.Append((T).3*TV(-1,-1).Normalized());}
    if(test_number==8){
        TV domain_center=grid.domain.Center();domain_center.y=(T)1;
        sources.Resize(1);source_velocity.Resize(1);world_to_source.Resize(1);
        sources(1)=BOX<TV>(0,(T).1,0,(T).2);
        ARRAY<MATRIX<T,3> > source_to_world(1);
        MATRIX<T,3> translation=MATRIX<T,3>::Translation_Matrix(domain_center);
        source_to_world(1)=translation;world_to_source(1)=source_to_world(1).Inverse();
        source_velocity(1)=TV(0,-.5);}

    // set example-specific parameters
    fluids_parameters.object_friction=(test_number==4||test_number==20)?(T)1:0;
    if(test_number==8){
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        use_two_way_coupling_for_sph=false;
        convert_sph_particles_to_fluid=false;
        use_variable_density_for_sph=true;}
    else if(test_number==9){ 
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        convert_sph_particles_to_fluid=false;
        use_analytic_divergence=true;
        use_two_way_coupling_for_sph=true;
        use_variable_density_for_sph=true;}
    else if(test_number==10){ 
        flip_ratio=1;
        target_particles_per_unit_volume=fluids_parameters.number_particles_per_cell/grid.Cell_Size();
        neumann_boundary_slip_multiplier=0;
        adjust_cell_weights_on_neumann_boundaries=true;
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        convert_sph_particles_to_fluid=false;
        use_analytic_divergence=true;
        use_analytic_divergence_for_expansion_only=false;
        particle_targeting_time=.25;
        use_two_way_coupling_for_sph=false;
        use_variable_density_for_sph=true;
        enforce_density_near_interface=true;}

    if(test_number==4||test_number==5||test_number==6||test_number==20){fluids_parameters.solid_affects_fluid=true;fluids_parameters.fluid_affects_solid=false;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> WATER_STANDARD_TESTS_2D<T_GRID>::
~WATER_STANDARD_TESTS_2D()
{}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Initialize_Advection(const bool always_use_objects)
{
    PHYSBAM_ASSERT(!always_use_objects,"fluids_parameters.solid_affects_fluids overrides this functionality");  // Deprecated parameter
    if(always_use_objects||test_number==4||test_number==5||test_number==6||test_number==20) fluids_parameters.Use_Fluid_Coupling_Defaults();
    else fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,2> WATER_STANDARD_TESTS_2D<T_GRID>::
Initial_Velocity(const TV& X) const
{
    if(test_number==7 || test_number==11) return TV(0,1);
    if(test_number==12 || test_number==13) return TV(-1,0);
    return TV();
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Get_Variable_Viscosity(ARRAY<T,VECTOR<int,2> >& variable_viscosity,const T time) const
{
    if(test_number==13){
        for(typename GRID<TV>::CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.2)?(T)0:(T)1000;}
    else if(test_number==14){
        for(typename GRID<TV>::CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.5)?(T)0:(T)1000;}
    else if(test_number==15){
        for(typename GRID<TV>::CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.5)?(T)0:(T)1000;}
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR WATER_STANDARD_TESTS_2D<T_GRID>::
Initial_Phi(const TV& X) const
{
    T phi=1;
    if(test_number==1) phi=X.y-(T).32424;
    else if(test_number==2){
        static SPHERE<TV> circle((TV((T).5,(T).75)),(T).2);
        phi=min(circle.Signed_Distance(X),X.y-(T).412134);}
    else if(test_number==4||test_number==20) phi=X.y-(T).400235234;
    else if(test_number==6){
        static SPHERE<TV> circle((TV((T).7,(T).14)),(T).1);
        phi=circle.Signed_Distance(X);}
    else if(test_number==7){
        static SPHERE<TV> circle((TV((T).5,(T).5)),(T).2);
        phi=circle.Signed_Distance(X);}
    else if(test_number==9){
        static BOX<TV> box((TV((T)3.4,1)),(TV((T)3.6,(T)1.3)));
        phi=min(box.Signed_Distance(X),X.y-(T).4);}
    else if(test_number==11){
        static TV normal((T)1,(T)1),location((T).5,(T).75);
        static LINE_2D<T> line(normal,location);
        phi=line.Signed_Distance(X);}
    else if(test_number==12 || test_number==13){
        static SPHERE<TV> circle((TV((T).5,(T).75)),(T).2);
        phi=circle.Signed_Distance(X);}
    else if(test_number==14){
        static BOX<TV> box((TV((T).25,(T).75)),(TV((T)0,(T).5)));
        phi=box.Signed_Distance(X);}
    else if(test_number==15){
        static SPHERE<TV> circle((TV((T).5,(T)1)),(T).2);
        phi=circle.Signed_Distance(X);}
    else if(test_number==21) phi=X.y-(T).32424;
    for(int s=1;s<=sources.m;s++) phi=min(phi,sources(s).Signed_Distance(world_to_source(s).Homogeneous_Times(X)));
    return phi;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Initialize_Bodies()
{
    if(test_number==4){
        sphere=rigid_body_collection.Add_Rigid_Body(example.stream_type,example.data_directory+"/Rigid_Bodies_2D/circle",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particle.X(sphere)=TV((T)1.25,(T).55);
        rigid_body_collection.rigid_body_particle.kinematic(sphere)=true;}
    else if (test_number==20){
        sphere=rigid_body_collection.Add_Rigid_Body(example.stream_type,example.data_directory+"/Rigid_Bodies_2D/circle",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particle.X(sphere)=TV((T).8,(T).1);
        rigid_body_collection.rigid_body_particle.kinematic(sphere)=true;}
    if(use_inaccurate_body_collisions){
        inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection.rigid_geometry_collection);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);}
    else fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);
}
//#####################################################################
// Function Update_Sources
//#####################################################################
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Update_Sources(const T time)
{
    if(test_number==3)
        if(time>4) sources.Clean_Memory();
    if(test_number==5)
        if(time>(T)1.6) sources.Clean_Memory();
    if(test_number==8)
        if(time>2) sources.Clean_Memory();
    if(test_number==10){ //TODO : move this to a better callback function
        FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >&>(fluids_parameters);
        SPH_EVOLUTION_UNIFORM<GRID<TV> >& sph_evolution=*fluids_parameters_uniform.sph_evolution;
        if(time>12){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume/2;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}
        else if(time>8){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}
        else if(time>4){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume*2;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}}
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if((test_number==4||test_number==20) && id==sphere) frame.t=motion_curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class T_GRID> bool WATER_STANDARD_TESTS_2D<T_GRID>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if((test_number==4||test_number==20) && id==sphere){twist.linear=motion_curve.Derivative(time);return true;}
    return false;
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Initialize_SPH_Particles()
{
    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >&>(fluids_parameters);
    PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& particle_levelset=fluids_parameters_uniform.particle_levelset_evolution->particle_levelset;
    SPH_EVOLUTION_UNIFORM<T_GRID>& sph_evolution=*fluids_parameters_uniform.sph_evolution;
    T_GRID& grid=fluids_parameters_uniform.particle_levelset_evolution->grid;

    sph_evolution.use_two_way_coupling=use_two_way_coupling_for_sph;
    sph_evolution.use_variable_density_solve=use_variable_density_for_sph;
    sph_evolution.convert_particles_to_fluid=convert_sph_particles_to_fluid;
    sph_evolution.use_analytic_divergence=use_analytic_divergence;
    sph_evolution.use_analytic_divergence_for_expansion_only=use_analytic_divergence_for_expansion_only;
    sph_evolution.flip_ratio=flip_ratio;
    sph_evolution.neumann_boundary_slip_multiplier=neumann_boundary_slip_multiplier;
    sph_evolution.adjust_cell_weights_on_neumann_boundaries=adjust_cell_weights_on_neumann_boundaries;
    sph_evolution.enforce_density_near_interface=enforce_density_near_interface;
    sph_evolution.particle_targeting_time=particle_targeting_time;
    sph_evolution.target_particles_per_unit_volume=fluids_parameters.number_particles_per_cell/grid.Cell_Size();
    sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;

    int particle_id=0,number_of_sph_particles=0;

    if(test_number==8){
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        BOX<TV> particle_region((TV((T).2,1)),TV((T).3,(T)1.2));
        number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume);
        for(int i=0;i<number_of_sph_particles;i++){
            TV X=random.Get_Uniform_Vector(particle_region);
            TV_INT block=grid.Block_Index(X,3);
            if(!removed_negative_particles(block)) removed_negative_particles(block)=particle_levelset.template_removed_particles.Clone();
            int id=removed_negative_particles(block)->array_collection->Add_Element();
            (*removed_negative_particles(block)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID))(id)=particle_id++;
            removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.Minimum_Edge_Length();}}
    else if(test_number==9){
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        BOX<TV> particle_region((TV((T).4,1)),TV((T).6,(T)1.3));
        T particle_multiplier=3;
        number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume/particle_multiplier);
        for(int region=0;region<3;region++){
            for(int i=0;i<number_of_sph_particles;i++){
                TV X=random.Get_Uniform_Vector(particle_region);
                TV_INT block=grid.Block_Index(X,3);
                if(!removed_negative_particles(block)) removed_negative_particles(block)=particle_levelset.template_removed_particles.Clone();
                int id=removed_negative_particles(block)->array_collection->Add_Element();
                (*removed_negative_particles(block)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID))(id)=particle_id++;
                removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.Minimum_Edge_Length();}
            particle_region+=TV(1,0);number_of_sph_particles=int(particle_multiplier*number_of_sph_particles);}}
    else if(test_number==10){
        sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume;
        sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        BOX<TV> particle_region(TV(grid.domain.min_corner.x,grid.domain.min_corner.y),TV(grid.domain.max_corner.x,grid.domain.max_corner.y/3));
        T particle_multiplier=1;
        number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume/particle_multiplier);
        for(int i=0;i<number_of_sph_particles;i++){
            TV X=random.Get_Uniform_Vector(particle_region);
            TV_INT block=grid.Block_Index(X,3);
            if(!removed_negative_particles(block)) removed_negative_particles(block)=particle_levelset.template_removed_particles.Clone();
            int id=removed_negative_particles(block)->array_collection->Add_Element();
            (*removed_negative_particles(block)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID))(id)=particle_id++;
            removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.Minimum_Edge_Length();}}
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
template<class T_GRID> void WATER_STANDARD_TESTS_2D<T_GRID>::
Limit_Dt(T& dt,const T time)
{
    if(test_number==4||test_number==20){
        TV velocity=rigid_body_collection.rigid_body_particle.V(sphere);
        T rigid_dt_denominator=abs(velocity.x)/grid.dX.x+abs(velocity.y)/grid.dX.y;
        if(rigid_dt_denominator>1e-8) dt=min(dt,1/rigid_dt_denominator);}
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR WATER_STANDARD_TESTS_2D<T_GRID>::
Initial_Phi_Object(const TV& X) const
{
    if(test_number==4||test_number==20) return rigid_body_collection.Rigid_Body(sphere).Implicit_Geometry_Extended_Value(X);
    else if(test_number==5){
        static BOX<TV> glass((T).005,(T).095,(T)-.1,1);
        return -glass.Signed_Distance(X);}
    else if(test_number==6){
        static LINE_2D<T> line((TV((T)-.2,1)).Normalized(),TV());
        return line.Signed_Distance(X);}
    return 1;
}
//#####################################################################
template class WATER_STANDARD_TESTS_2D<GRID<VECTOR<float,2> > >;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::Get_Variable_Viscosity(ARRAY<float,VECTOR<int,2> >&,float) const;
template float WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::Initial_Phi(VECTOR<float,2> const&) const;
template VECTOR<float,2> WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::Initial_Velocity(VECTOR<float,2> const&) const;
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::Initialize_Advection(bool);
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::Limit_Dt(float&,float);
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::Update_Sources(float);
template WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::WATER_STANDARD_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,2> >&,FLUIDS_PARAMETERS<RLE_GRID_2D<float> >&,
    RIGID_BODY_COLLECTION<VECTOR<float,2> >&,int,int);
template WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::~WATER_STANDARD_TESTS_2D();
template float WATER_STANDARD_TESTS_2D<RLE_GRID_2D<float> >::Initial_Phi_Object(VECTOR<float,2> const&) const;
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class WATER_STANDARD_TESTS_2D<GRID<VECTOR<double,2> > >;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::Get_Variable_Viscosity(ARRAY<double,VECTOR<int,2> >&,double) const;
template double WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::Initial_Phi(VECTOR<double,2> const&) const;
template VECTOR<double,2> WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::Initial_Velocity(VECTOR<double,2> const&) const;
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::Initialize_Advection(bool);
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::Limit_Dt(double&,double);
template void WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::Update_Sources(double);
template WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::WATER_STANDARD_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,2> >&,FLUIDS_PARAMETERS<RLE_GRID_2D<double> >&,
    RIGID_BODY_COLLECTION<VECTOR<double,2> >&,int,int);
template WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::~WATER_STANDARD_TESTS_2D();
template double WATER_STANDARD_TESTS_2D<RLE_GRID_2D<double> >::Initial_Phi_Object(VECTOR<double,2> const&) const;
#endif
#endif
