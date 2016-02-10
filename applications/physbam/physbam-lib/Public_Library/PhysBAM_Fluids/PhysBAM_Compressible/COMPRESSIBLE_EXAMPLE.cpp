//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_1D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_G.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_G.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_H.h>
// #include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/COMPRESSIBLE_ADVECTION_ENO_RF.h>
// #include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Solvers/COMPRESSIBLE_ADVECTION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/COMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Read_Write/Particles/READ_WRITE_RIGIDS_PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_EXAMPLE<TV>::
COMPRESSIBLE_EXAMPLE(const STREAM_TYPE& stream_type,const int array_collection_type): 
    BASE(stream_type),mpi_grid(0),resolution(0),solid_affects_fluid(false),eos(new EOS_GAMMA<T>()),boundary(0),boundary_scalar(0),
    conservation_law_solver(0),kinematic_coupling_utilities(0),collision_bodies_affecting_fluid(0),solid_body_collection(*new SOLID_BODY_COLLECTION<TV>(this,array_collection_type))
{
    LOG::Initialize_Logging(false,false,1<<30,true,1);
    Initialize_Geometry_Particle();
    Initialize_Rigids_Particles();
    Initialize_Read_Write_Structures();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMPRESSIBLE_EXAMPLE<TV>::
~COMPRESSIBLE_EXAMPLE()
{
    if(mpi_grid) delete boundary;
    for(int i=1;i<=TV::dimension;i++) delete cfl_eigensystem[i];
    delete conservation_law_solver;
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf("/%d",frame);

    // Restart data
    if(frame==first_frame){FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
        if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);}
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/euler_U",U);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/euler_psi",psi);

    // Detailed output data
    T_ARRAYS_SCALAR rho(grid.Domain_Indices()),velocity(grid.Domain_Indices()),pressure(grid.Domain_Indices()),energy(grid.Domain_Indices()),internal_energy(grid.Domain_Indices());
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)){
            rho(cell_index)=U(cell_index)(1);
            pressure(cell_index)=eos->p(rho(cell_index),EULER<T_GRID>::e(U(cell_index)));
            internal_energy(cell_index)=EULER<T_GRID>::e(U(cell_index));
            energy(cell_index)=EULER<T_GRID>::Get_Total_Energy(U,cell_index);
            velocity(cell_index)=EULER<T_GRID>::Get_Velocity_Component(U(cell_index),1);}}

    if(write_debug_data){T_FACE_ARRAYS_SCALAR density_flux(grid);
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face_index=iterator.Face_Index();
            if(conservation_law_solver->fluxes.Valid_Index(iterator.Full_Index())){int axis=iterator.Axis();
                density_flux.Component(axis)(face_index)=conservation_law_solver->fluxes.Component(axis)(face_index)(1);}}
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/density_flux",density_flux);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/psi_N",psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/mac_velocities",face_velocities);}

    solid_body_collection.Write(stream_type,output_directory,frame,frame,true,true,false,true,false);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/density",rho);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/energy",energy);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/pressure",pressure);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/internal_energy",internal_energy);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/centered_velocities",velocity);
}
//#####################################################################
// Function Clamp_Dt_Adaptively
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Clamp_Dt_Adaptively(T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_BOOL& psi,T_ARRAYS_SCALAR& rho_dt,T_ARRAYS_SCALAR& e_dt,const T& dt,T& clamp_rho,T& clamp_e)
{
    if(adaptive_time_step) for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)){T clamp_rho_cell=clamp_rho*U(cell_index)(1);
            if(rhs(cell_index)(1)>0 && U(cell_index)(1)>=clamp_rho_cell) rho_dt(cell_index)=(U(cell_index)(1)-clamp_rho_cell)/rhs(cell_index)(1);
            else rho_dt(cell_index)=dt;

            T e=EULER<T_GRID>::e(U,cell_index);int d=T_GRID::dimension+2;
            T clamp_e_cell=clamp_e*e,tmp_dt=dt,momentum_flux_sqr=0,momentum_sqr=0,momentum_flux_dot_product=0;
            for(int axis=1;axis<=TV::dimension;axis++){
                momentum_flux_sqr+=rhs(cell_index)(axis+1)*rhs(cell_index)(axis+1);
                momentum_sqr+=U(cell_index)(axis+1)*U(cell_index)(axis+1);
                momentum_flux_dot_product+=U(cell_index)(axis+1)*rhs(cell_index)(axis+1);}

                T rho_np1=U(cell_index)(1)-rho_dt(cell_index)*rhs(cell_index)(1);
                T rho_np1_sqr=rho_np1*rho_np1;
                PHYSBAM_ASSERT(rho_np1>0);
                T E_np1=U(cell_index)(d)-rho_dt(cell_index)*rhs(cell_index)(d);
                T mom_np1_sqr=momentum_sqr-2*rho_dt(cell_index)*momentum_flux_dot_product+rho_dt(cell_index)*rho_dt(cell_index)*momentum_flux_sqr;
                T e_np1=E_np1/rho_np1-(T).5*mom_np1_sqr/rho_np1_sqr;

                T a=2*rhs(cell_index)(1)*rhs(cell_index)(d)-momentum_flux_sqr-2*clamp_e_cell*rhs(cell_index)(1)*rhs(cell_index)(1);
                T c=2*U(cell_index)(d)*U(cell_index)(1)-2*clamp_e_cell*U(cell_index)(1)*U(cell_index)(1)-momentum_sqr;

                LOG::cout.precision(20);
                if(e_np1<clamp_e_cell){
                    T b_over_two=momentum_flux_dot_product-U(cell_index)(1)*rhs(cell_index)(d)-U(cell_index)(d)*rhs(cell_index)(1)+2*clamp_e_cell*U(cell_index)(1)*rhs(cell_index)(1);
                    T b_sqr_over_four=b_over_two*b_over_two,ac=a*c;

                    if(b_sqr_over_four>ac){
                        if((a>=0 && a<1e-16 && b_over_two<1e-16 && b_over_two>=0)||(a<=0 && a>-1e-16 && b_over_two>-1e-16 && b_over_two<=0)) tmp_dt=dt;
                        else if((a>=0 && a<1e-16)||(a<=0 && a>-1e-16)){T tmp_dt1=(-(T).5*c)/b_over_two;tmp_dt=(tmp_dt1>0)?tmp_dt1:dt;}
                        else if((b_over_two>=0 && b_over_two<1e-16)||(b_over_two<=0 && b_over_two>-1e-16)){T tmp_dt1=abs(sqrt(-c/a));tmp_dt=(tmp_dt1>0)?tmp_dt1:dt;}
                        else{T tmp_dt1=(-1*b_over_two-sqrt(b_sqr_over_four-ac))/a;
                            if(tmp_dt1==(T)0.) tmp_dt=(T)1e-16;
                        else{T tmp_dt2=c/(-1*b_over_two-sqrt(b_sqr_over_four-ac));
                            tmp_dt=(tmp_dt1>0)?(tmp_dt2>0)?min(tmp_dt1,tmp_dt2):tmp_dt1:(tmp_dt2>0)?tmp_dt2:dt;}}}
                    else tmp_dt=dt;
                PHYSBAM_ASSERT(tmp_dt>0);
                e_dt(cell_index)=min(dt,tmp_dt);}}}
}
//#####################################################################
// Function Compute_Dt_For_Internal_Energy
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Compute_Dt_For_Internal_Energy(T_GRID& grid,T_FACE_ARRAYS_DIMENSION_SCALAR& fluxes,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,TV_INT& face_index,TV_INT& cell_index,const bool first,T& clamp_e,T& rho_dt,T& e_dt,const int axis)
{
    int d=T_GRID::dimension+2;
    T rho=U_ghost(cell_index)(1)/pow((T)2.,T_GRID::dimension),momentum_sqr=(T)0.,E=U_ghost(cell_index)(d)/pow((T)2.,T_GRID::dimension);
    for(int complete_axis=1;complete_axis<=T_GRID::dimension;complete_axis++) momentum_sqr+=U_ghost(cell_index)(complete_axis+1)*U_ghost(cell_index)(complete_axis+1);
    momentum_sqr/=pow((T)4.,T_GRID::dimension);
    T e=E/rho-(T).5*momentum_sqr/(rho*rho),clamp_e_cell=clamp_e*e/pow((T)2.,T_GRID::dimension);
    PHYSBAM_ASSERT(e>0);
    T density_flux=fluxes(axis,face_index)(1)/grid.DX()(axis),momentum_flux=fluxes(axis,face_index)(axis+1)/grid.DX()(axis),energy_flux=fluxes(axis,face_index)(d)/grid.DX()(axis);

    if(!first){density_flux*=(T)-1.;momentum_flux*=(T)-1.;energy_flux*=(T)-1.;}

    T momentum_flux_sqr=momentum_flux*momentum_flux,momentum_flux_dot_product=U_ghost(cell_index)(axis+1)*momentum_flux/pow((T)2.,T_GRID::dimension);
    T rho_np1=rho-rho_dt*density_flux,momentum_np1_sqr=momentum_sqr-(T)2.*rho_dt*momentum_flux_dot_product+rho_dt*rho_dt*momentum_flux_sqr,E_np1=E-rho_dt*energy_flux;
    PHYSBAM_ASSERT(rho_np1>0);
    if(rho_np1<1e-16){e_dt=(T)0.;return;}

    T e_np1=E_np1/rho_np1-(T).5*momentum_np1_sqr/(rho_np1*rho_np1);

    T a=(T)2.*density_flux*energy_flux-momentum_flux_sqr-(T)2.*clamp_e_cell*density_flux*density_flux;
    T c=(T)2.*E*rho-(T)2.*clamp_e_cell*rho*rho-momentum_sqr;

    if(e_np1<clamp_e_cell){
        T b_over_two=momentum_flux_dot_product-rho*energy_flux-E*density_flux+(T)2.*clamp_e_cell*rho*density_flux;
        T b_sqr_over_four=b_over_two*b_over_two,ac=a*c;
        if(b_sqr_over_four>ac){
            if((a>=0 && a<1e-16 && b_over_two<1e-16 && b_over_two>=0)||(a<=0 && a>-1e-16 && b_over_two>-1e-16 && b_over_two<=0)) e_dt=rho_dt;
            else if((a>=0 && a<1e-16)||(a<=0 && a>-1e-16)){T tmp_dt1=(-(T).5*c)/b_over_two;e_dt=(tmp_dt1>0)?tmp_dt1:rho_dt;}
            else if((b_over_two>=0 && b_over_two<1e-16)||(b_over_two<=0 && b_over_two>-1e-16)){T tmp_dt1=abs(sqrt(-c/a));e_dt=(tmp_dt1>0)?tmp_dt1:rho_dt;}
            else{T tmp_dt1=(-b_over_two-sqrt(b_over_two*b_over_two-ac))/a;
                if(tmp_dt1==0) e_dt=(T)0.;
                else{T tmp_dt2=c/(-b_over_two-sqrt(b_over_two*b_over_two-ac));
                    e_dt=(tmp_dt1>=0)?(tmp_dt2>=0)?min(tmp_dt1,tmp_dt2):tmp_dt1:(tmp_dt2>=0)?tmp_dt2:rho_dt;}}}
        else e_dt=rho_dt;
        PHYSBAM_ASSERT(e_dt>=0);}

    if(e_dt<(T)1e-16) e_dt=(T)0.;
    T tmp_dt=min(rho_dt,e_dt);
    rho_np1=rho-tmp_dt*density_flux,momentum_np1_sqr=momentum_sqr-(T)2.*tmp_dt*momentum_flux_dot_product+tmp_dt*tmp_dt*momentum_flux_sqr,E_np1=E-tmp_dt*energy_flux;
    PHYSBAM_ASSERT(rho_np1>0);
    e_np1=E_np1/rho_np1-(T).5*momentum_np1_sqr/(rho_np1*rho_np1);
    PHYSBAM_ASSERT(e_np1>0);
}
//#####################################################################
// Function Clamp_Fluxes
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Clamp_Fluxes(T_GRID& grid,T_FACE_ARRAYS_DIMENSION_SCALAR& fluxes,T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T& dt,T& clamp_rho,T& clamp_e)
{
    //TODO: FIXME (causes a section type conflict)
    /*if(adaptive_time_step){for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        T alpha=(T)1.,rho_dt=dt,e_dt1=dt,e_dt2=dt;
        TV_INT outgoing_cell_index=(fluxes(axis,face_index)(1)>0)?iterator.First_Cell_Index():iterator.Second_Cell_Index();
        if(psi.Valid_Index(outgoing_cell_index)&&psi(outgoing_cell_index)){
            T clamp_rho_cell=((T)1.-clamp_rho)*U_ghost(outgoing_cell_index)(1)/pow((T)2.,T_GRID::dimension),flux_value=abs(fluxes(axis,face_index)(1));
            if(dt*flux_value/grid.DX()(axis)>clamp_rho_cell && flux_value!=(T)0.) rho_dt=grid.DX()(axis)*clamp_rho_cell/flux_value;}

        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        if(psi.Valid_Index(first_cell_index)&&psi(first_cell_index)) Compute_Dt_For_Internal_Energy(grid,fluxes,U_ghost,face_index,first_cell_index,true,clamp_e,rho_dt,e_dt1,axis);
        if(psi.Valid_Index(second_cell_index)&&psi(second_cell_index)) Compute_Dt_For_Internal_Energy(grid,fluxes,U_ghost,face_index,second_cell_index,false,clamp_e,rho_dt,e_dt2,axis);
        alpha=min(rho_dt,min(e_dt1,e_dt2))/dt;PHYSBAM_ASSERT(alpha<=1.);if(alpha<(T)1.) fluxes(axis,face_index)*=alpha;}

        LOG::cout.precision(20);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            if(psi(cell_index)){rhs(cell_index)=TV_DIMENSION(); 
                for(int axis=1;axis<=T_GRID::dimension;axis++) rhs(cell_index)+=(fluxes(axis,iterator.Second_Face_Index(axis))-fluxes(axis,iterator.First_Face_Index(axis)))/grid.DX()(axis);
                TV_DIMENSION U_np1=U_ghost(cell_index)-dt*rhs(cell_index);
                PHYSBAM_ASSERT(U_np1(1)>0);
                T momentum_sqr=(T)0.;for(int axis=1;axis<=T_GRID::dimension;axis++) momentum_sqr+=U_np1(axis+1)*U_np1(axis+1);
                PHYSBAM_ASSERT((T)2.*U_np1(T_GRID::dimension+2)*U_np1(1)>momentum_sqr);}}}*/
}
//#####################################################################
// Initialize_Grids
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Initialize_Grids()
{
    U.Resize(grid.Domain_Indices(),true,false);
    U_save.Resize(grid.Domain_Indices(),true,false);
    U_ghost.Resize(grid.Domain_Indices(3),true,false);
    face_velocities.Resize(grid.Domain_Indices(3),true,false);
    psi.Resize(grid.Domain_Indices(3),true,false);
    psi_N.Resize(grid.Domain_Indices(3),true,false);
}
//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("/%d",frame);

    // Restart data
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/common/grid",grid);
    if(mpi_grid) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+f+"/euler_U",U);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+f+"/euler_psi",psi);

    solid_body_collection.Read(stream_type,output_directory,frame,frame,true,true,false,true);
}
//#####################################################################
// Limit_Dt
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Limit_Dt(T& dt,const T time)
{
    TV max_lambdas=TV::Constant_Vector((T)0.);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)) for(int k=1;k<=T_GRID::dimension;++k){
            max_lambdas(k)=max(max_lambdas(k),cfl_eigensystem(k)->Maximum_Magnitude_Eigenvalue(U(cell_index)));}}
    TV max_lambda_over_dx=max_lambdas*grid.one_over_dX;

    dt=min(dt,cfl/max_lambda_over_dx.Sum());
}
//#####################################################################
// Log_Parameters
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Log_Parameters() const
{
    BASE::Log_Parameters();
    LOG::SCOPE scope("COMPRESSIBLE_EXAMPLE parameters");
    LOG::cout<<"CFL Number = "<<cfl<<std::endl;
    LOG::cout<<"Global Minimum dt = "<<global_min_dt<<std::endl;
    LOG::cout<<"Resolution = "<<resolution<<std::endl;
    LOG::cout<<"ENO Order = "<<eno_order<<std::endl;
    LOG::cout<<"R-K Order = "<<rk_order<<std::endl;
    LOG::cout<<"Use ENO RF = "<<use_rf<<std::endl;
    LOG::cout<<"Use adaptive time stepping = "<<adaptive_time_step<<std::endl;
    LOG::cout<<"Use flux clamping = "<<clamp_fluxes<<std::endl;
}
//#####################################################################
// Register_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Register_Options()
{
    BASE::Register_Options();
    if(!parse_args) return;

    // Custom Stuff
    parse_args->Add_Double_Argument("-cfl",(T).5,"Compressible Flow CFL number.");
    parse_args->Add_Double_Argument("-g",(T)1e-8,"Global minimum value for dt.");
    parse_args->Add_Integer_Argument("-resolution",100,"Resolution of the fluid mesh (must be divisible by 2!).");
    parse_args->Add_Integer_Argument("-eno_order",3,"Spatial order of accuracy for flux calculation.");
    parse_args->Add_Integer_Argument("-rk_order",3,"Temporal order of accuracy for time stepping.");
    parse_args->Add_Option_Argument("-rf","Use ENO RF.");
    parse_args->Add_Option_Argument("-write_debug_data","Write debug data.");
    parse_args->Add_Option_Argument("-adaptive_time_step","Use adaptive time stepping.");
    parse_args->Add_Option_Argument("-clamp_fluxes","For clamping fluxes along with adaptive time stepping.");
}
//#####################################################################
// Parse_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Parse_Options()
{
    BASE::Parse_Options();
    resolution=parse_args->Get_Integer_Value("-resolution");
    eno_order=parse_args->Get_Integer_Value("-eno_order");
    rk_order=parse_args->Get_Integer_Value("-rk_order");
    cfl=(T)parse_args->Get_Double_Value("-cfl");
    global_min_dt=(T)parse_args->Get_Double_Value("-g");
    use_rf=parse_args->Is_Value_Set("-rf");
    adaptive_time_step=parse_args->Is_Value_Set("-adaptive_time_step");
    write_debug_data=parse_args->Is_Value_Set("-write_debug_data");
    clamp_fluxes=parse_args->Is_Value_Set("-clamp_fluxes");
}
//#####################################################################
// Save_State
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Save_State()
{
    T_ARRAYS_DIMENSION_SCALAR::Copy(U,U_save);
}
//#####################################################################
// Restore_State
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Restore_State()
{
    T_ARRAYS_DIMENSION_SCALAR::Copy(U_save,U);
}
//#####################################################################
// Set_Eigensystems
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Set_Eigensystems()
{
    Set_Eigensystems_Helper(cfl_eigensystem);
}
//#####################################################################
// Function Set_Eigensystems_Helper
//#####################################################################
template<class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,1>& eigensystems)
{
    if(eigensystems[1]) delete eigensystems[1];eigensystems[1]=new EULER_1D_EIGENSYSTEM_F<T>();
}
template<class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2>& eigensystems)
{
    if(eigensystems[1]) delete eigensystems[1];eigensystems[1]=new EULER_2D_EIGENSYSTEM_F<T>();
    if(eigensystems[2]) delete eigensystems[2];eigensystems[2]=new EULER_2D_EIGENSYSTEM_G<T>();
}
template<class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,5> >*,3>& eigensystems)
{
    if(eigensystems[1]) delete eigensystems[1];eigensystems[1]=new EULER_3D_EIGENSYSTEM_F<T>();
    if(eigensystems[2]) delete eigensystems[2];eigensystems[2]=new EULER_3D_EIGENSYSTEM_G<T>();
    if(eigensystems[3]) delete eigensystems[3];eigensystems[3]=new EULER_3D_EIGENSYSTEM_H<T>();
}
//#####################################################################
// Parse_Late_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Parse_Late_Options()
{
    BASE::Parse_Late_Options();
    Set_Eigensystems();

    if(use_rf) conservation_law_solver=new CONSERVATION_ENO_RF<GRID<TV>,TV::dimension+2>();
    else conservation_law_solver=new CONSERVATION_ENO_LLF<GRID<TV>,TV::dimension+2>();
    conservation_law_solver->Set_Order(eno_order);
    conservation_law_solver->adaptive_time_step=adaptive_time_step;
    kinematic_coupling_utilities=new COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>(grid,U,psi,mpi_grid);
    Initialize_Grids();
}
//#####################################################################
template class COMPRESSIBLE_EXAMPLE<VECTOR<float,1> >;
template class COMPRESSIBLE_EXAMPLE<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_EXAMPLE<VECTOR<double,1> >;
template class COMPRESSIBLE_EXAMPLE<VECTOR<double,2> >;
#endif
