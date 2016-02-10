//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace MPI_UTILITIES
//#####################################################################
#ifndef __MPI_UTILITIES__
#define __MPI_UTILITIES__

#ifdef USE_MPI

#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
namespace PhysBAM{

template<class TV> struct POINT_FACE_REPULSION_PAIR;
template<class TV> struct EDGE_EDGE_REPULSION_PAIR;

namespace MPI_UTILITIES{

template<class T,int d> struct DATATYPE_HELPER<POINT_FACE_REPULSION_PAIR<VECTOR<T,d> > >{static MPI::Datatype Datatype();};
template<class T,int d> struct DATATYPE_HELPER<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,d> > >{static MPI::Datatype Datatype();};

}
}
#endif
#endif
