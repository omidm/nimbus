//#####################################################################
// Copyright 2006-2007, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_3D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;
template<class T_GRID> SMOKE_STANDARD_TESTS_3D<T_GRID>::
SMOKE_STANDARD_TESTS_3D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container,RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :example(example),fluids_parameters(fluids_parameters),incompressible_fluid_container(incompressible_fluid_container),rigid_body_collection(rigid_body_collection),test_number(0)
{
}
template<class T_GRID> SMOKE_STANDARD_TESTS_3D<T_GRID>::
~SMOKE_STANDARD_TESTS_3D()
{
}
template<class T_GRID> void SMOKE_STANDARD_TESTS_3D<T_GRID>::
Initialize(const int test_number_input,const int resolution)
{
    test_number=test_number_input;
    std::stringstream ss;ss<<"Running Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;LOG::filecout(ss.str());
    
    // set up the standard fluid environment
    // TODO: *REALLY* need to pick sensible constants and settings
    example.frame_rate=24;
    fluids_parameters.cfl=3;
    fluids_parameters.domain_walls=VECTOR<VECTOR<bool,2>,T_GRID::dimension>::Constant_Vector(VECTOR<bool,2>::Constant_Vector(false));
    fluids_parameters.domain_walls(2)(1)=true;
    fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).15;
    fluids_parameters.kolmogorov=(T)0;fluids_parameters.gravity=(T)0;
    rho=(T)1;fluids_parameters.rho_bottom=(T)1;fluids_parameters.rho_top=(T).65;
    fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=0;
    fluids_parameters.temperature_container.Set_Cooling_Constant(0);fluids_parameters.temperature_products=(T)3000;
    fluids_parameters.write_velocity=true;example.write_frame_title=true;fluids_parameters.write_debug_data=true;
    fluids_parameters.use_body_force=true;
    if(test_number==3 || test_number==4) fluids_parameters.temperature_buoyancy_constant=(T)0.001;
        
    //set up the domain
    int cells=resolution;
    if(test_number==1||test_number==2||test_number==3||test_number==4){
        example.first_frame=0;example.last_frame=3840;
        grid.Initialize(10*cells+1,10*cells+1,25*cells+1,-10,10,-10,10,-5,45);}
    //grid.Initialize(10*cells+1,15*cells+1,10*cells+1,0,1,0,1.5,0,1);}
    else{LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}
    
    example.output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests_Smoke/Test_%d__Resolution_%d_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1),(grid.counts.z-1));
    
    // set up the source domain
    if(test_number==1||test_number==2||test_number==3||test_number==4){
        source=BOX<TV>((T).35,(T).65,(T).35,(T).65,(T)0,(T).1);
        //source=BOX<TV>((T).45,(T).55,(T)0,(T).1,(T).45,(T).55);
        world_to_source=MATRIX<T,4>::Identity_Matrix();
        source_velocity=TV((T)0,(T)0,(T)2);}
    //source_velocity=TV((T)0,(T)0.5,(T)0);}
    
    //set up example-specific parameters
    if(test_number==3 || test_number==4){
        explosion_divergence=100;explosion_end_time=3;fluids_parameters.use_non_zero_divergence=true;}
    if(test_number==4){
        vortex_particle_evolution=new VORTEX_PARTICLE_EVOLUTION_3D<T>();
        vortex_particle_evolution->Initialize(grid);
        vortex_particle_evolution->particle_confinement_parameter=(T)1;
        vortex_particle_evolution->renormalize_vorticity_after_stretching_tilting=true;
        vortex_particle_evolution->force_scaling=(T)1e-3;
        particle_vorticity_minimum_density=(T).9;source_vorticity_magnitude=(T)1000;particle_radius=(T).04;
        fluids_parameters.use_vorticity_confinement=false;}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_GRID> void SMOKE_STANDARD_TESTS_3D<T_GRID>::
Initialize_Bodies()
{
    if(test_number==2){
        int id=rigid_body_collection.Add_Rigid_Body(example.stream_type,example.data_directory+"/Rigid_Bodies/sphere",(T).125,true,true,false);
        rigid_body_collection.rigid_body_particle.X(id)=TV((T).5,(T).75,(T).5);
        rigid_body_collection.Rigid_Body(id).is_static=true;
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);}
}
//#####################################################################
// Function Get_Divergence
//#####################################################################
template<class T_GRID> void SMOKE_STANDARD_TESTS_3D<T_GRID>::
Get_Divergence(ARRAY<T,VECTOR<int,3> >& divergence,const T dt,const T time)
{
    LOG::Time("Getting divergence");
    if(test_number==3 || test_number==4){
        T expansion=explosion_divergence*sin(time)/exp(time);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) 
            if(source.Lazy_Inside(iterator.Location())) divergence(iterator.Cell_Index())=expansion;}
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
template<class T_GRID> void SMOKE_STANDARD_TESTS_3D<T_GRID>::
Get_Body_Force(ARRAY<T,FACE_INDEX<3> >& force,const T dt,const T time)
{
    if(test_number==4){
        ARRAY<T,FACE_INDEX<3> > face_velocities_ghost(*fluids_parameters.grid,fluids_parameters.number_of_ghost_cells,false);
        fluids_parameters.incompressible->boundary->Fill_Ghost_Cells_Face(*fluids_parameters.grid,incompressible_fluid_container.face_velocities,face_velocities_ghost,time,fluids_parameters.number_of_ghost_cells);
        vortex_particle_evolution->Compute_Body_Force(face_velocities_ghost,force,dt,time);
        if(time<=explosion_end_time){
            int add_count=0;VORTICITY_PARTICLES<TV >& vorticity_particles=vortex_particle_evolution->vorticity_particles;
            TV cell_upper=(T).5*fluids_parameters.grid->dX,cell_lower=-cell_upper;
            for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
                if(source.Lazy_Inside(iterator.Location()) && time>(T)1/24 && random.Get_Uniform_Number((T)0,(T)1)<(T).005){
                    std::stringstream ss;ss<<"adding particle now have "<<vorticity_particles.array_collection->Size()+1<<std::endl;LOG::filecout(ss.str());
                    add_count++;int particle_id=vorticity_particles.array_collection->Add_Element(); 
                    vorticity_particles.radius(particle_id)=particle_radius;
                    vorticity_particles.X(particle_id)=iterator.Location()+random.Get_Uniform_Vector(cell_lower,cell_upper);
                    vorticity_particles.vorticity(particle_id)=(T)source_vorticity_magnitude*TV::Cross_Product(TV(0,1,0),(vorticity_particles.X(particle_id)-source.Center()).Normalized()).Normalized();}}
            vortex_particle_evolution->Euler_Step(face_velocities_ghost,dt,time);}
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,3> SMOKE_STANDARD_TESTS_3D<T_GRID>::
Initial_Velocity(const TV& X) const
{
    return TV();
}
template class SMOKE_STANDARD_TESTS_3D<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SMOKE_STANDARD_TESTS_3D<GRID<VECTOR<double,3> > >;
#endif
