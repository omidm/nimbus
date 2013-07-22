//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION_CUT
//#####################################################################
#ifndef __FLUID_TO_SOLID_INTERPOLATION_CUT__
#define __FLUID_TO_SOLID_INTERPOLATION_CUT__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_BASE.h>

namespace PhysBAM{

template<class T> class SEGMENTED_CURVE_2D;
template<class TV> class GENERALIZED_FLUID_MASS;
template<class TV> class MATRIX_FLUID_GRADIENT_CUT;
template<class TV>
class FLUID_TO_SOLID_INTERPOLATION_CUT:public FLUID_TO_SOLID_INTERPOLATION_BASE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef FLUID_TO_SOLID_INTERPOLATION_BASE<TV> BASE;
public:
    using BASE::index_map;
    enum FACE_TYPE {unused,inside,outside,cut,coupling};

    struct ENTRY
    {
        TV w;
        int i;
    };

    SEGMENTED_CURVE_2D<T>& curve;

    ARRAY<ARRAY<ENTRY> > entries;

    struct CLIP_ENTRY
    {
        int i;
        T a,b;
    };
    struct CUT_CELL
    {
        ARRAY<CLIP_ENTRY> clipped_segments;
        int face;
        int other_cell;
        SEGMENT_2D<T> segment;
        T inside_volume;
    };
    HASHTABLE<TV_INT,CUT_CELL> cut_cells;
    GENERALIZED_FLUID_MASS<TV>** fluid_mass;
    MATRIX_FLUID_GRADIENT_CUT<TV>* gradient;
    T density;
    T outside_density;
    HASHTABLE<int> unused_faces;
    HASHTABLE<FACE_INDEX<TV::m>,T> face_lengths;
    ARRAY<bool,TV_INT>* outside_fluid;
    const ARRAY<bool,FACE_INDEX<TV::m> >* psi_N;
    bool use_cut_volume;

    FLUID_TO_SOLID_INTERPOLATION_CUT(const COLLISION_AWARE_INDEX_MAP<TV>& map,SEGMENTED_CURVE_2D<T>& curve,T density_input);
    virtual ~FLUID_TO_SOLID_INTERPOLATION_CUT();
//#####################################################################
    virtual void Setup_Before_Compute(ARRAY<bool,TV_INT>& outside_fluid_input,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_input);
    void Compute_Dirichlet_Cells();
    void Remove_Degeneracy();
    void Compute_Gradient();
    void Add_Gradient_Entry(int fi,const FACE_INDEX<TV::m>& f,int side,bool outside);
    void Add_Cut_Gradient_Entry(int fi,const FACE_INDEX<TV::m>& f,int side);
    void Compute_Beta();
    FACE_TYPE Face_Type(int f) const;
    void Compute(const int ghost_cells) PHYSBAM_OVERRIDE;
    void Times_Add(const VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const PHYSBAM_OVERRIDE;
    void Transpose_Times_Add(const GENERALIZED_VELOCITY<TV>& solid_force,VECTOR_ND<T>& fluid_force) const PHYSBAM_OVERRIDE;
    void Print_Each_Matrix(int n,int fluid_faces,GENERALIZED_VELOCITY<TV>& G) const PHYSBAM_OVERRIDE;
    void Cut_Face(FACE_INDEX<TV::m>& f,const TV& normal,const SEGMENT_2D<T>& segment);
    void Fill_Extra_Velocities(VECTOR_ND<T>& fluid_velocity_vector) const;
    void Dump_Extra_Velocities(const VECTOR_ND<T>& fluid_velocity_vector);
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
