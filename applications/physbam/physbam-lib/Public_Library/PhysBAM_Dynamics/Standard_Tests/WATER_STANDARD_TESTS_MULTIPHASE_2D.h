//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_MULTIPHASE_2D
//#####################################################################
// Test descriptions:
//   11. Two drops colliding
//   12. Sphere splashing into a pool of two liquid phases
//   13. 4 phase splash
//   14. Rising air bubble in water
// Also supports a variety of standard resolutions in powers of 2.
//#####################################################################
#ifndef __WATER_STANDARD_TESTS_MULTIPHASE_2D__
#define __WATER_STANDARD_TESTS_MULTIPHASE_2D__

#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_MULTIPHASE.h>
namespace PhysBAM{

template<class T_GRID>
class WATER_STANDARD_TESTS_MULTIPHASE_2D:public WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,WATER_STANDARD_TESTS_2D<T_GRID> >
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<T,2> TV;
public:
    typedef WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,WATER_STANDARD_TESTS_2D<T_GRID> > BASE;
    using BASE::rigid_body_collection;using BASE::fluids_parameters;using BASE::grid;using BASE::example;using BASE::sphere;using BASE::test_number;using BASE::world_to_source;
    using BASE::source_velocity;using BASE::source_region;using BASE::sources;

    WATER_STANDARD_TESTS_MULTIPHASE_2D(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>& example_input,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
        :WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,WATER_STANDARD_TESTS_2D<T_GRID> >(example_input,fluids_parameters_input,incompressible_fluid_container,rigid_body_collection_input)
    {
    }

void Initialize(const int test_number_input,const int resolution,const int restart_frame=-1)
{
    BASE::Initialize(test_number_input,resolution,restart_frame);
    std::stringstream ss;ss<<"Running Multiphase Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;LOG::filecout(ss.str());
    int cells=1*resolution;
    if(test_number==11){
        fluids_parameters.domain_walls[2][2]=true;
        example.first_frame=0;example.last_frame=500;example.frame_rate=250;
        grid.Initialize(10*cells+1,10*cells+1,0,(T).1,0,(T).1);}
    if(test_number==12){
        fluids_parameters.domain_walls[2][2]=false;
        grid.Initialize(20*cells+1,15*cells+1,(T)-.5,(T)1.5,(T)-.5,1);}
    if(test_number==13){
        fluids_parameters.domain_walls[2][2]=true;
        example.first_frame=0;example.last_frame=500;example.frame_rate=100;
        grid.Initialize(10*cells+1,10*cells+1,0,2,0,2);}
    if(test_number==14){
        fluids_parameters.domain_walls[2][2]=true;
        example.first_frame=0;example.last_frame=20;example.frame_rate=400;
        grid.Initialize(20*cells+1,30*cells+1,(T)-.01,(T).01,(T)-.01,(T).02);}
    if(test_number==15){
        fluids_parameters.domain_walls[2][2]=true;
        example.first_frame=0;example.last_frame=500;example.frame_rate=60;
        grid.Initialize(20*cells+1,10*cells+1,0,1,0,(T).5);
        world_to_source.Append(MATRIX<T,3>::Identity_Matrix());
        sources.Append(BOX<TV>(TV((T).75,(T).45),TV((T).875,(T).6)));
        source_velocity.Append(TV(0,(T)-.6));
        source_region.Append(3);
        world_to_source.Append(MATRIX<T,3>::Identity_Matrix());
        sources.Append(BOX<TV>(TV((T).6,(T).45),TV((T).725,(T).6)));
        source_velocity.Append(TV(0,0));
        source_region.Append(4);}
    if(test_number==16){
        example.first_frame=0;example.last_frame=500;example.frame_rate=120;
        grid.Initialize(10*cells+1,11*cells+1,0,1,0,(T)1.1);}
    if(test_number==17){
        example.first_frame=0;example.last_frame=500;example.frame_rate=250;
        grid.Initialize(10*cells+1,10*cells+1,0,(T).1,0,(T).1);}
    example.output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests_Multiphase/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));
    std::stringstream ss1;ss1<<"output directory="<<example.output_directory<<std::endl;LOG::filecout(ss1.str());
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,WATER_STANDARD_TESTS_2D<T_GRID> >::Initialize_Bodies();
    if(test_number==15){
        int ground=rigid_body_collection.Add_Rigid_Body(example.stream_type,example.data_directory+"/Rigid_Bodies_2D/ground",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particle.X(ground)=TV((T).1,0);
        rigid_body_collection.rigid_body_particle.rotation(ground)=ROTATION<TV>::From_Angle((T)pi/8);
        rigid_body_collection.rigid_body_particle.kinematic(ground)=true;}
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);
}
//#####################################################################
};
}
#endif
