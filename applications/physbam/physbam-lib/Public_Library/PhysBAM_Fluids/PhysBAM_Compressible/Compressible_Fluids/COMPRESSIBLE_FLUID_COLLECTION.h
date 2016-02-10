//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_FLUID_COLLECTION
//#####################################################################
#ifndef __COMPRESSIBLE_FLUID_COLLECTION__
#define __COMPRESSIBLE_FLUID_COLLECTION__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{

template<class T_GRID>
class COMPRESSIBLE_FLUID_COLLECTION:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    const T_GRID& grid;

    EOS<T>* eos;
    T_ARRAYS_BOOL psi;
    T_ARRAYS_DIMENSION_SCALAR U;
        
    COMPRESSIBLE_FLUID_COLLECTION(const T_GRID& grid_input);
    ~COMPRESSIBLE_FLUID_COLLECTION();

    void Set_Equation_Of_State(EOS<T>* eos_input)
    {eos=eos_input;}
    
    //#####################################################################
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const;
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame);
    void Initialize_Grids();
    //#####################################################################
};
}
#endif
