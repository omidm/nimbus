//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_AWARE_INDEX_MAP
//#####################################################################
#ifndef __COLLISION_AWARE_INDEX_MAP__
#define __COLLISION_AWARE_INDEX_MAP__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/SIDED_FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;
template<class TV> class IMPLICIT_BOUNDARY_CONDITION_COLLECTION;
template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;

template<class TV>
class COLLISION_AWARE_INDEX_MAP:public NONCOPYABLE
{
    enum WORKAROUND {d=TV::dimension};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;

    int last_coupling_cell;

public:
    UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& iterator_info;
    GRID<TV>& grid;
    ARRAY<int,FACE_INDEX<d> > face_indices;
    HASHTABLE<SIDED_FACE_INDEX<d>,int> constraint_indices;
    ARRAY<int,TV_INT> cell_indices;

    ARRAY<FACE_INDEX<d> > indexed_faces;
    ARRAY<SIDED_FACE_INDEX<d> > indexed_constraints;
    ARRAY<TV_INT> indexed_cells;
    ARRAY<int> real_cell_indices,real_cell_indices_reverse_map;
    ARRAY<int> mpi_face_indices;
    int number_extra_cells;
    bool two_phase;

    IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>& boundary_condition_collection;

    COLLISION_AWARE_INDEX_MAP(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>& boundary_condition_collection_input);

    int Number_Faces() const
    {return indexed_faces.m+indexed_constraints.m+number_extra_cells;}

    int Number_Cells() const
    {return two_phase?indexed_cells.m+number_extra_cells:indexed_cells.m;}

//#####################################################################
    void Register_Dirichlet_Cell(const TV_INT& index);
    void Clear(const int ghost_cells);
    void Construct_Indices(const int ghost_cells);
    bool Register_Face_Index(const FACE_INDEX<d>& index);
    bool Register_Constrained_Face_Index(const SIDED_FACE_INDEX<d>& index);
    void Register_Cell_Index(const TV_INT& index,const int ghost_cells);
    void Collect(const ARRAY<T,FACE_INDEX<d> >& faces,const ARRAY<T>& constrained_faces,VECTOR_ND<T>& flattened_faces) const;
    void Collect(const ARRAY<T,TV_INT>& cells,VECTOR_ND<T>& flattened_cells) const;
    void Collect_Indexed_Cells(const ARRAY<T,TV_INT>& cells,VECTOR_ND<T>& flattened_cells) const;
    void Collect_Boundary_Faces(const ARRAY<T,FACE_INDEX<d> >& cells,VECTOR_ND<T>& flattened_cells) const;
    void Distribute(const VECTOR_ND<T>& flattened_faces,ARRAY<T,FACE_INDEX<d> >& faces,ARRAY<T>& constrained_faces) const;
    void Distribute(const VECTOR_ND<T>& flattened_cells,ARRAY<T,TV_INT>& cells) const;
    void Distribute_Indexed_Cells(const VECTOR_ND<T>& flattened_cells,ARRAY<T,TV_INT>& cells) const;
    void Distribute_Boundary_Faces(const VECTOR_ND<T>& flattened_cells,ARRAY<T,FACE_INDEX<d> >& cells) const;
    void Print(int id) const;
//#####################################################################
};
}
#endif
