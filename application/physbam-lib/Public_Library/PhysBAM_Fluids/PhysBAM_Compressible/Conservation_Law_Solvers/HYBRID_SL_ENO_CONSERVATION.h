//#####################################################################
// Copyright 2010-2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_SL_ENO_CONSERVATION
//##################################################################### 
#ifndef __HYBRID_SL_ENO_CONSERVATION__
#define __HYBRID_SL_ENO_CONSERVATION__   

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
namespace PhysBAM{

template<class T_GRID> class EULER_EIGENSYSTEM;
template<class T,class TV_DIMENSION> class EIGENSYSTEM;
template<class TV> class GRID;

template<class T_GRID,int d>
class HYBRID_SL_ENO_CONSERVATION:public CONSERVATION<T_GRID,d>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::INDEX INDEX;
    typedef VECTOR<bool,2*T_GRID::dimension> TV_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef CONSERVATION<T_GRID,d> BASE;

    const T_FACE_ARRAYS_BOOL& flux_face;
    CONSERVATION<T_GRID,d> *conservation;

public:
    HYBRID_SL_ENO_CONSERVATION(const T_FACE_ARRAYS_BOOL& flux_face_in, CONSERVATION<T_GRID,d> *conservation_law_solver)
        : flux_face(flux_face_in),conservation(conservation_law_solver)
    {conservation->Save_Fluxes();}

    void Set_Order(const int order_input=3)
    {conservation->Set_Order(order_input);}
            
    void Use_Field_By_Field_Alpha()
    {conservation->Use_Field_By_Field_Alpha();}

    void Use_Maximum_Alpha()
    {conservation->Use_Maximum_Alpha();}

    void Set_Use_Exact_Neumann_Face_Location(const bool use_exact_neumann_face_location_input)
    {conservation->Set_Use_Exact_Neumann_Face_Location(use_exact_neumann_face_location_input);}

    void Amplify_Alpha(const T amplification_factor_input=1)
    {conservation->Amplify_Alpha(amplification_factor_input);}

    void Save_Fluxes()
    {conservation->Save_Fluxes();}

    void Scale_Outgoing_Fluxes_To_Clamp_Variable(bool scale_outgoing_fluxes_to_clamp_variable_input,int clamped_variable_index_input,T clamped_value_input)
    {conservation->Scale_Outgoing_Fluxes_To_Clamp_Variable(scale_outgoing_fluxes_to_clamp_variable_input, clamped_variable_index_input, clamped_value_input);}

    void Set_Callbacks(CONSERVATION_CALLBACKS<T_GRID,TV_DIMENSION> *callbacks_input)
    {conservation->Set_Callbacks(callbacks_input);}

    void Set_Custom_Object_Boundary(BOUNDARY_OBJECT<T_GRID,TV_DIMENSION>& object_boundary_input)
    {conservation->Set_Custom_Object_Boundary(object_boundary_input);}

    virtual void Log_Parameters() const
    {BASE::Log_Parameters();conservation->Log_Parameters();}

    virtual void Update_Conservation_Law(T_GRID& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T dt,
        VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems_explicit,const T_FACE_ARRAYS_BOOL& psi_N,
        const T_FACE_ARRAYS_SCALAR& face_velocities,const bool thinshell=false,const TV_BOOL& outflow_boundaries=TV_BOOL(),VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>* eigensystems_auxiliary=0,
        T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary=0);
};
}
#endif
