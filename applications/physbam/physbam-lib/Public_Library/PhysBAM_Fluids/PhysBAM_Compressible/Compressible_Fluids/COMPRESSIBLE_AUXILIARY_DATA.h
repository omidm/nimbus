//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace COMPRESSIBLE_AUXILIARY_DATA
//#####################################################################
#ifndef __COMPRESSIBLE_AUXILIARY_DATA__
#define __COMPRESSIBLE_AUXILIARY_DATA__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Compressible_Fluids/COMPRESSIBLE_FLUID_COLLECTION.h>
#include <string>

namespace PhysBAM{
class STREAM_TYPE;
template<class T_GRID> class COMPRESSIBLE_FLUID_COLLECTION;
template<class T> class EOS;

namespace COMPRESSIBLE_AUXILIARY_DATA{

template<class T_GRID,class T_ARRAYS,class T_ARRAYS_BOOL_INPUT,class T_FACE_ARRAYS_DIMENSION_SCALAR>
void Write_Auxiliary_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame,const T_GRID& grid,
    const int number_of_ghost_cells,const T_ARRAYS& U,const T_ARRAYS_BOOL_INPUT& psi,const EOS<typename T_GRID::SCALAR>& eos,const bool write_debug_data,const T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes);

template<class T_GRID>
void Write_Auxiliary_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame,const COMPRESSIBLE_FLUID_COLLECTION<T_GRID>& compressible_fluid_collection,const bool write_debug_data);

}
}
#endif
