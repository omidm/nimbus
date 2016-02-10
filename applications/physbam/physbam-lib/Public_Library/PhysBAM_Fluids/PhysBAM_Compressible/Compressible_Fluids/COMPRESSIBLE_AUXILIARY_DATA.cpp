//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/GRADIENT_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Compressible_Fluids/COMPRESSIBLE_AUXILIARY_DATA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Compressible_Fluids/COMPRESSIBLE_FLUID_COLLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
namespace PhysBAM{
namespace COMPRESSIBLE_AUXILIARY_DATA{
template<class T_GRID,class T_ARRAYS,class T_ARRAYS_BOOL_INPUT,class T_FACE_ARRAYS_DIMENSION_SCALAR>
void Write_Auxiliary_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame,const T_GRID& grid,
    const int number_of_ghost_cells,const T_ARRAYS& U,const T_ARRAYS_BOOL_INPUT& psi,const EOS<typename T_GRID::SCALAR>& eos,const bool write_debug_data,const T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;

    STATIC_ASSERT((IS_SAME<T_ARRAYS,T_ARRAYS_DIMENSION_SCALAR>::value));
    STATIC_ASSERT((IS_SAME<T_ARRAYS_BOOL_INPUT,T_ARRAYS_BOOL>::value));

    RANGE<TV_INT> domain_indices=grid.Domain_Indices(number_of_ghost_cells);
    T_ARRAYS_VECTOR velocity(domain_indices),velocity_plus_c,velocity_minus_c,momentum;
    if(write_debug_data){velocity_plus_c.Resize(domain_indices);velocity_minus_c.Resize(domain_indices);momentum.Resize(domain_indices);}
    T_ARRAYS_SCALAR density(domain_indices),pressure(domain_indices),energy,internal_energy,temperature,
        entropy,enthalpy,speedofsound,machnumber,density_gradient,pressure_gradient;
    if(write_debug_data){energy.Resize(domain_indices);internal_energy.Resize(domain_indices);temperature.Resize(domain_indices);entropy.Resize(domain_indices);
        enthalpy.Resize(domain_indices);speedofsound.Resize(domain_indices);machnumber.Resize(domain_indices);density_gradient.Resize(domain_indices);pressure_gradient.Resize(domain_indices);}

    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(grid.Domain_Indices().Lazy_Inside(cell) && !psi(cell)) continue;
        density(cell)=EULER<T_GRID>::Get_Density(U,cell);
        pressure(cell)=eos.p(density(cell),EULER<T_GRID>::e(U(cell)));
        velocity(cell)=EULER<T_GRID>::Get_Velocity(U,cell);
        if(write_debug_data){
            momentum(cell)=velocity(cell)*density(cell);
            energy(cell)=EULER<T_GRID>::Get_Total_Energy(U,cell);
            entropy(cell)=eos.S(density(cell),eos.e_From_p_And_rho(pressure(cell),density(cell)));
            enthalpy(cell)=EULER<T_GRID>::enthalpy(eos,U(cell));
            internal_energy(cell)=EULER<T_GRID>::e(U(cell));
            temperature(cell)=eos.T(density(cell),internal_energy(cell));
            speedofsound(cell)=eos.c(density(cell),eos.e_From_p_And_rho(pressure(cell),density(cell)));
            if(speedofsound(cell)) machnumber(cell)=velocity(cell).Magnitude()/speedofsound(cell);
            else machnumber(cell)=(T)-1; // output non-physical negative value when machnumber is infinite
            velocity_plus_c(cell)=velocity(cell)+speedofsound(cell)*TV::All_Ones_Vector();
            velocity_minus_c(cell)=velocity(cell)-speedofsound(cell)*TV::All_Ones_Vector();}}

    T_FACE_ARRAYS_SCALAR density_flux(grid),momentum_flux(grid),energy_flux(grid);
    if(write_debug_data && fluxes){
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT face_index=iterator.Face_Index();
            if(fluxes->Valid_Index(iterator.Full_Index())){int axis=iterator.Axis();
                density_flux.Component(axis)(face_index)=fluxes->Component(axis)(face_index)(1);
                momentum_flux.Component(axis)(face_index)=fluxes->Component(axis)(face_index)(2);
                energy_flux.Component(axis)(face_index)=fluxes->Component(axis)(face_index)(T_GRID::dimension+2);}}}

    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(write_debug_data){
        GRADIENT::Compute_Magnitude(grid,number_of_ghost_cells,density,density_gradient);
        GRADIENT::Compute_Magnitude(grid,number_of_ghost_cells,pressure,pressure_gradient);
        
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density_flux",density_flux);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/momentum_flux",momentum_flux);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/energy_flux",energy_flux);

        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/momentum",momentum);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/energy",energy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density_gradient",density_gradient);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure_gradient",pressure_gradient);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/entropy",entropy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/enthalpy",enthalpy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/speedofsound",speedofsound);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/machnumber",machnumber);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/internal_energy",internal_energy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/temperature",temperature);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/velocity_plus_c",velocity_plus_c);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/velocity_minus_c",velocity_minus_c);}

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",pressure);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/centered_velocities",velocity);
}
template<class T_GRID>
void Write_Auxiliary_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame,const COMPRESSIBLE_FLUID_COLLECTION<T_GRID>& compressible_fluid_collection,const bool write_debug_data)
{
    typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;

    const T_GRID& grid=compressible_fluid_collection.grid;
    const T_ARRAYS_DIMENSION_SCALAR& U=compressible_fluid_collection.U;
    const EOS<T>* eos=compressible_fluid_collection.eos;
    const T_FACE_ARRAYS_DIMENSION_SCALAR *fluxes = NULL;

    Write_Auxiliary_Files(stream_type,output_directory,frame,grid,0,U,compressible_fluid_collection.psi,*eos,write_debug_data,fluxes);
}
template void Write_Auxiliary_Files<GRID<VECTOR<float,1> >,ARRAY<VECTOR<float,3>,VECTOR<int,1> >,ARRAY<bool,VECTOR<int,1> >,ARRAY<VECTOR<float,3>,FACE_INDEX<1> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<float,1> > const&,int,ARRAY<VECTOR<float,3>,VECTOR<int,1> > const&,ARRAY<bool,VECTOR<int,1> > const&,EOS<GRID<VECTOR<float,1> >::SCALAR> const&,const bool,const ARRAY<VECTOR<float,3>,FACE_INDEX<1> >*);
template void Write_Auxiliary_Files<GRID<VECTOR<float,2> >,ARRAY<VECTOR<float,4>,VECTOR<int,2> >,ARRAY<bool,VECTOR<int,2> >,ARRAY<VECTOR<float,4>,FACE_INDEX<2> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<float,2> > const&,int,ARRAY<VECTOR<float,4>,VECTOR<int,2> > const&,ARRAY<bool,VECTOR<int,2> > const&,EOS<GRID<VECTOR<float,2> >::SCALAR> const&,const bool,const ARRAY<VECTOR<float,4>,FACE_INDEX<2> >*);
template void Write_Auxiliary_Files<GRID<VECTOR<float,3> >,ARRAY<VECTOR<float,5>,VECTOR<int,3> >,ARRAY<bool,VECTOR<int,3> >,ARRAY<VECTOR<float,5>,FACE_INDEX<3> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<float,3> > const&,int,ARRAY<VECTOR<float,5>,VECTOR<int,3> > const&,ARRAY<bool,VECTOR<int,3> > const&,EOS<GRID<VECTOR<float,3> >::SCALAR> const&,const bool,const ARRAY<VECTOR<float,5>,FACE_INDEX<3> >*);
template void Write_Auxiliary_Files<GRID<VECTOR<float,1> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<float,1> > > const&,const bool);
template void Write_Auxiliary_Files<GRID<VECTOR<float,2> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<float,2> > > const&,const bool);
template void Write_Auxiliary_Files<GRID<VECTOR<float,3> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<float,3> > > const&,const bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Write_Auxiliary_Files<GRID<VECTOR<double,1> >,ARRAY<VECTOR<double,3>,VECTOR<int,1> >,ARRAY<bool,VECTOR<int,1> >,ARRAY<VECTOR<double,3>,FACE_INDEX<1> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<double,1> > const&,int,ARRAY<VECTOR<double,3>,VECTOR<int,1> > const&,ARRAY<bool,VECTOR<int,1> > const&,EOS<GRID<VECTOR<double,1> >::SCALAR> const&,const bool,const ARRAY<VECTOR<double,3>,FACE_INDEX<1> >*);
template void Write_Auxiliary_Files<GRID<VECTOR<double,2> >,ARRAY<VECTOR<double,4>,VECTOR<int,2> >,ARRAY<bool,VECTOR<int,2> >,ARRAY<VECTOR<double,4>,FACE_INDEX<2> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<double,2> > const&,int,ARRAY<VECTOR<double,4>,VECTOR<int,2> > const&,ARRAY<bool,VECTOR<int,2> > const&,EOS<GRID<VECTOR<double,2> >::SCALAR> const&,const bool,const ARRAY<VECTOR<double,4>,FACE_INDEX<2> >*);
template void Write_Auxiliary_Files<GRID<VECTOR<double,3> >,ARRAY<VECTOR<double,5>,VECTOR<int,3> >,ARRAY<bool,VECTOR<int,3> >,ARRAY<VECTOR<double,5>,FACE_INDEX<3> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,GRID<VECTOR<double,3> > const&,int,ARRAY<VECTOR<double,5>,VECTOR<int,3> > const&,ARRAY<bool,VECTOR<int,3> > const&,EOS<GRID<VECTOR<double,3> >::SCALAR> const&,const bool,const ARRAY<VECTOR<double,5>,FACE_INDEX<3> >*);
template void Write_Auxiliary_Files<GRID<VECTOR<double,1> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<double,1> > > const&,const bool);
template void Write_Auxiliary_Files<GRID<VECTOR<double,2> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<double,2> > > const&,const bool);
template void Write_Auxiliary_Files<GRID<VECTOR<double,3> > >(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,int,COMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<double,3> > > const&,const bool);
#endif
}
}
