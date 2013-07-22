//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_MULTIPHASE.h>
using namespace PhysBAM;
template<class T_GRID,class T_WATER_STANDARD_TESTS> WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
WATER_STANDARD_TESTS_MULTIPHASE(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>& example,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container_input,
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :T_WATER_STANDARD_TESTS(example,fluids_parameters_input,rigid_body_collection_input),fluids_parameters_uniform(fluids_parameters_input),
    incompressible_fluid_container(incompressible_fluid_container_input)
{
}
template<class T_GRID,class T_WATER_STANDARD_TESTS> WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
~WATER_STANDARD_TESTS_MULTIPHASE()
{
}
template<class T_GRID,class T_WATER_STANDARD_TESTS> void WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Initialize(const int test_number_input,const int resolution,const int restart_frame)
{
    BASE::Initialize(Non_Multiphase_Test_Number(test_number_input),resolution);
    test_number=test_number_input;
    fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=false;
    fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
    test_number=test_number_input;
    fluids_parameters.incompressible_iterations=100;

    use_open_wall=false;air_region=-1;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input=static_cast<FLUIDS_PARAMETERS_UNIFORM<T_GRID>&>(fluids_parameters);

    if(test_number<10){
        fluids_parameters_input.densities(1)=1000;fluids_parameters_input.densities(2)=1;
        source_region.Resize(sources.m);ARRAYS_COMPUTATIONS::Fill(source_region,1);
        fluids_parameters_input.pseudo_dirichlet_regions(2)=true;
        fluids_parameters_input.second_order_cut_cell_method=false;}
    else if(test_number==11){
        fluids_parameters_input.densities(1)=1000;fluids_parameters_input.densities(2)=500;fluids_parameters_input.densities(3)=2000;
        fluids_parameters_input.surface_tensions(1,2)=fluids_parameters_input.surface_tensions(2,1)=(T).1;
        fluids_parameters_input.surface_tensions(1,3)=fluids_parameters_input.surface_tensions(3,1)=(T).1;
        fluids_parameters_input.surface_tensions(2,3)=fluids_parameters_input.surface_tensions(3,2)=(T).1;}
    else if(test_number==12){
        fluids_parameters_input.dirichlet_regions(1)=true;use_open_wall=true;air_region=1;
        fluids_parameters_input.densities(1)=1;fluids_parameters_input.densities(2)=800;fluids_parameters_input.densities(3)=1000;fluids_parameters_input.densities(4)=3000;}
    else if(test_number==13){
        fluids_parameters_input.densities(1)=1400;fluids_parameters_input.densities(2)=500;fluids_parameters_input.densities(3)=1000;fluids_parameters_input.densities(4)=1;}
    else if(test_number==14){
        fluids_parameters_input.densities(1)=(T)1.226;fluids_parameters_input.densities(2)=1000; 
        fluids_parameters_input.surface_tensions(1,2)=fluids_parameters_input.surface_tensions(2,1)=(T).0728;
        fluids_parameters_input.viscosities(1)=(T).0000178;
        fluids_parameters_input.viscosities(2)=(T).001137;
        fluids_parameters_input.implicit_viscosity=false;
        fluids_parameters_input.incompressible_iterations=200;
        fluids_parameters_input.implicit_viscosity_iterations=200;}
    else if(test_number==15){
        fluids_parameters_input.densities(1)=1000;
        fluids_parameters_input.viscosities(1)=(T)500;           
        //fluids_parameters_input.viscosities(1)=(T)50;           
        //fluids_parameters_input.use_multiphase_strain(1)=true;
        //fluids_parameters_input.elastic_moduli(1)=20000;
        //fluids_parameters_input.plasticity_alphas(1)=0;

        fluids_parameters_input.densities(2)=1000;
        fluids_parameters_input.viscosities(2)=(T)60;
            
        fluids_parameters_input.densities(3)=1000;
        fluids_parameters_input.viscosities(3)=10;
        // SOURCE SET IN DERIVED EXAMPLE

        fluids_parameters_input.densities(4)=1000;
        fluids_parameters_input.viscosities(4)=0;           
        // SOURCE SET IN DERIVED EXAMPLE
            
        fluids_parameters_input.densities(5)=1;
        fluids_parameters_input.dirichlet_regions(5)=true;//use_open_wall=true;air_region=6;
        //fluids_parameters_input.cfl/=8;

        fluids_parameters_input.implicit_viscosity_iterations=50;
        fluids_parameters_input.implicit_viscosity=true;
        //fluids_parameters_input.cfl=4;
    }
    else if(test_number==16){
        if(restart_frame>=500){
            fluids_parameters_input.densities(1)=500;        
            fluids_parameters_input.use_multiphase_strain(1)=true;
            fluids_parameters_input.elastic_moduli(1)=15000;
            fluids_parameters_input.plasticity_alphas(1)=0;                
            fluids_parameters_input.implicit_viscosity_iterations=50;
            fluids_parameters_input.implicit_viscosity=true;
            fluids_parameters_input.viscosities(1)=(T)200;

            fluids_parameters_input.densities(4)=1000;          
            fluids_parameters_input.use_multiphase_strain(4)=true;
            fluids_parameters_input.elastic_moduli(4)=15000;
            fluids_parameters_input.plasticity_alphas(4)=0;
            fluids_parameters_input.viscosities(4)=(T)200;}
        else if(restart_frame>=296) fluids_parameters_input.densities(1)=500;         
        else if(restart_frame>=68){
            fluids_parameters_input.densities(1)=1500;   
            fluids_parameters_input.implicit_viscosity_iterations=50;
            fluids_parameters_input.implicit_viscosity=false;
            fluids_parameters_input.viscosities(1)=0;           
        }      
        else{
            fluids_parameters_input.densities(1)=500;        
            fluids_parameters_input.use_multiphase_strain(1)=true;
            fluids_parameters_input.elastic_moduli(1)=15000;
            fluids_parameters_input.plasticity_alphas(1)=0;                
            fluids_parameters_input.implicit_viscosity_iterations=50;
            fluids_parameters_input.implicit_viscosity=true;
            fluids_parameters_input.viscosities(1)=(T)200;           
        }
        fluids_parameters_input.densities(2)=1000;
        fluids_parameters_input.densities(3)=1;
        fluids_parameters_input.dirichlet_regions(3)=true;
        fluids_parameters_input.reseeding_frame_rate=10;
        fluids_parameters_input.cfl/=2;
    }
    else if(test_number==17){
        fluids_parameters_input.densities(1)=1000;
        fluids_parameters_input.densities(2)=(T)1.226;
        //fluids_parameters_input.surface_tensions(1,2)=fluids_parameters_input.surface_tensions(2,1)=(T)2;
        //fluids_parameters_input.surface_tensions(1,2)=fluids_parameters_input.surface_tensions(2,1)=(T)2;
        fluids_parameters_input.cfl/=2;
        fluids_parameters_input.reseeding_frame_rate=1;
    }
    else{
        std::stringstream ss;ss<<"Unrecognized example: "<<test_number<<std::endl;LOG::filecout(ss.str());
        PHYSBAM_FATAL_ERROR();}


    if(test_number==15){fluids_parameters_input.solid_affects_fluid=true;fluids_parameters_input.fluid_affects_solid=false;}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class T_GRID,class T_WATER_STANDARD_TESTS> void WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Initialize_Advection(const bool always_use_objects)   
{
    if(always_use_objects||test_number==4||test_number==5||test_number==6||test_number==15) fluids_parameters.Use_Fluid_Coupling_Defaults();
    else fluids_parameters.Use_No_Fluid_Coupling_Defaults();

    if(use_open_wall)
        for(int i=1;i<=Number_Of_Regions(test_number);i++){
            BOUNDARY_PHI_WATER<T_GRID>* boundary=new BOUNDARY_PHI_WATER<T_GRID>();
            boundary->Set_Velocity_Pointer(incompressible_fluid_container.face_velocities);
            if(i==air_region)boundary->sign=-1;
            fluids_parameters.phi_boundary_multiphase(i)=boundary;}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_GRID,class T_WATER_STANDARD_TESTS> void WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Initialize_Bodies()
{
    T_WATER_STANDARD_TESTS::Initialize_Bodies();
    if(test_number==16){
        GRID<VECTOR<T,3> > armadillo_temp_grid;
        ARRAY<T,VECTOR<int,3> > armadillo_temp_phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > armadillo_temp(armadillo_temp_grid,armadillo_temp_phi);
        FILE_UTILITIES::Read_From_File<float>(example.data_directory+"/Rigid_Bodies/armadillo_high_res.phi",armadillo_temp);
        armadillo=new T_LEVELSET(*(new T_GRID(*fluids_parameters.grid)),*(new T_ARRAYS_SCALAR(fluids_parameters.grid->Domain_Indices(1))));
        for(CELL_ITERATOR iterator(*fluids_parameters.grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            VECTOR<T,3> vec=VECTOR<T,3>(iterator.Location());exchange(vec.x,vec.z);
            if(T_GRID::dimension==3){armadillo->phi(cell)=armadillo_temp.Extended_Phi((vec-VECTOR<T,3>((T).5,(T).35,(T).5+(T).07*vec.y))*(T)145)/(T)145;}
            else{
                armadillo->phi(cell)=1;
                for(T i=0;i<=1;i+=(T).01)armadillo->phi(cell)=min(armadillo->phi(cell),armadillo_temp.Extended_Phi((VECTOR<T,3>(i,vec.y,vec.z)-VECTOR<T,3>((T).5,(T).35,(T).5+(T).07*vec.y))*(T)145)/(T)145);}}
    }
}
//#####################################################################
// Function Number_Of_Regions
//#####################################################################
template<class T_GRID,class T_WATER_STANDARD_TESTS> int WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Number_Of_Regions(int test_number)
{
    if(test_number<=5) return 2;
    if(test_number==11) return 3;
    if(test_number==12) return 4;
    if(test_number==13) return 4;
    if(test_number==14) return 2;
    if(test_number==15) return 5;
    if(test_number==16) return 4;
    if(test_number==17) return 2;
    std::stringstream ss;ss<<"Unrecognized example: "<<test_number<<std::endl;LOG::filecout(ss.str());
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Non_Multiphase_Test_Number
//#####################################################################
template<class T_GRID,class T_WATER_STANDARD_TESTS> int WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Non_Multiphase_Test_Number(int test_number)
{
    if(test_number<10) return test_number;
    if(test_number==11) return 1;
    if(test_number==12) return 4;
    if(test_number==13) return 1;
    if(test_number==14) return 1;
    if(test_number==15) return 1;
    if(test_number==16) return 1;
    if(test_number==17) return 1;
    std::stringstream ss;ss<<"Unrecognized example: "<<test_number<<std::endl;LOG::filecout(ss.str());
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
template<class T_GRID,class T_WATER_STANDARD_TESTS> typename T_GRID::VECTOR_T::SCALAR WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Initial_Phi(const int region,const TV& X) const
{
    if(test_number<10){
        if(region==1)return Initial_Phi(X);
        else if(region==2)return -Initial_Phi(X);}
    
    ARRAY<T> phis(50);
    // two drops test
    if(test_number==11){
        TV center1=TV(VECTOR<T,2>((T).055,(T).03)),center2=TV(VECTOR<T,2>((T).05,(T).07));
        if(T_GRID::dimension==3){center1[3]=(T).055;center2[3]=(T).05;}
        T radius=(T).015;
        phis(2)=(X-center1).Magnitude()-radius;
        phis(3)=(X-center2).Magnitude()-radius;
        phis(1)=-min(phis(2),phis(3));}
    // splash
    if(test_number==12){
        phis(2)=abs(X.y-(T).30)-(T).10;
        phis(3)=abs(X.y-(T).10)-(T).10;
        phis(4)=abs(X.y+(T).25)-(T).25;
        phis(1)=-min(phis(2),phis(3),phis(4));}
    // RT
    if(test_number==13){
        VECTOR<double,3> vec=VECTOR<double,3>(VECTOR<T,3>(X*(T)10));
        if(T_GRID::dimension==3) vec[3]=X[3]*10;
        T y_pos=X.y-(T)NOISE<double>::Noise1(vec,5,(T).5)*(T)0.002;
        phis(1)=abs(y_pos-(T).25)-(T).25;
        phis(2)=abs(y_pos-(T).75)-(T).25;
        phis(3)=abs(y_pos-(T)1.25)-(T).25;
        phis(4)=abs(y_pos-(T)1.75)-(T).25;}
    // rising bubble
    if(test_number==14){
        T radius=(T)1/(T)300;
        phis(2)=radius-X.Magnitude(); // center is at 0,0,0
        phis(1)=-phis(2);}
    // incline plane
    if(test_number==15){
        SPHERE<VECTOR<T,3> > sphere1((VECTOR<T,3>((T).4,(T).35,(T).35)),(T).1);if(T_GRID::dimension==2)sphere1.center.z=0;
        SPHERE<VECTOR<T,3> > sphere2((VECTOR<T,3>((T).4,(T).35,(T).65)),(T).1);if(T_GRID::dimension==2){sphere2.center.z=0;sphere2.center.x=(T).15;}
        phis(1)=sphere1.Signed_Distance(VECTOR<T,3>(X));
        phis(2)=sphere2.Signed_Distance(VECTOR<T,3>(X));
        phis(3)=1;
        phis(4)=1;
        phis(5)=-min(phis(1),phis(2),phis(3),phis(4));}
    // viscoelastic armadillo
    if(test_number==16){
        phis(1)=armadillo->Extended_Phi(X);
        phis(2)=max(X.y-(T).3,-phis(1));
        phis(3)=-min(phis(1),phis(2));
        phis(4)=1;}
    // milk crown
    if(test_number==17){
        TV center1=TV(VECTOR<T,2>((T).05,(T).03));if(T_GRID::dimension==3)center1[3]=(T).05;
        T radius=(T).01;
        phis(1)=min((X-center1).Magnitude()-radius,X.y-(T).01);
        phis(2)=-phis(1);}
    return phis(region);
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class T_GRID,class T_WATER_STANDARD_TESTS> typename T_GRID::VECTOR_T WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Initial_Velocity(const TV& X) const
{
    // milk crown
    if(test_number==17){
        TV center1=TV(VECTOR<T,2>((T).05,(T).03));if(T_GRID::dimension==3)center1[3]=(T).05;
        T radius=(T).011;
        if((X-center1).Magnitude()-radius<=0) return -TV::Axis_Vector(2);}
    return TV();
}
//#####################################################################
// Function Update_Sources
//#####################################################################
template<class T_GRID,class T_WATER_STANDARD_TESTS> void WATER_STANDARD_TESTS_MULTIPHASE<T_GRID,T_WATER_STANDARD_TESTS>::
Update_Sources(const T time)
{
    if(test_number==15){
//        if(time>(T)1.15) source_velocity(2)=TV(VECTOR<T,2>(0,(T)-.6));
        if(time>(T)1.15){source_region(1)=5;}
        if(time>(T)1.5){source_region(1)=4;}
    }
}
template class WATER_STANDARD_TESTS_MULTIPHASE<GRID<VECTOR<float,2> >,WATER_STANDARD_TESTS_2D<GRID<VECTOR<float,2> > > >;
template class WATER_STANDARD_TESTS_MULTIPHASE<GRID<VECTOR<float,3> >,WATER_STANDARD_TESTS_3D<GRID<VECTOR<float,3> > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class WATER_STANDARD_TESTS_MULTIPHASE<GRID<VECTOR<double,2> >,WATER_STANDARD_TESTS_2D<GRID<VECTOR<double,2> > > >;
template class WATER_STANDARD_TESTS_MULTIPHASE<GRID<VECTOR<double,3> >,WATER_STANDARD_TESTS_3D<GRID<VECTOR<double,3> > > >;
#endif
