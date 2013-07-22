//#####################################################################
// Copyright 2012, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MODEL_REDUCTION_SNAPSHOTS
// For model reduction via POD/SVD on a set of simulations
//#####################################################################
#ifndef _MODEL_REDUCTION_SNAPSHOTS_
#define _MODEL_REDUCTION_SNAPSHOTS_

#include <PhysBAM_Tools/Model_Reduction/MODEL_REDUCTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Model_Reduction/STATE_VECTOR.h>
#include <string.h>

namespace PhysBAM{

template<class TV>
class DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS:public MODEL_REDUCTION<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    
    typedef MODEL_REDUCTION<T> BASE;
    using BASE::reduced_basis;

    bool orthonormal_basis;
    MATRIX_MXN<T> reduced_basis_inverse;
public:
    using BASE::Get_Reduced_Basis;

    DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS(bool orthonormal_basis=true)
        :orthonormal_basis(orthonormal_basis)
    {};

    ~DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS()
    {};

//#####################################################################
    void Generate_Reduced_Basis(std::string folder_name,int rom_size);
    void Load_Reduced_Basis_From_Simulation_Folder(std::string folder_name,int rom_size);
    void Convert_Reduced_To_Full(VECTOR_ND<T>& v_reduced, STATE_VECTOR<TV>& v_full);
    void Convert_Full_To_Reduced(VECTOR_ND<T>& v_reduced, STATE_VECTOR<TV>& v_full);
    void Set_Reduced_Basis(MATRIX_MXN<T>& m_in) PHYSBAM_OVERRIDE;
    MATRIX_MXN<T>& Get_Reduced_Basis_Inverse(){return reduced_basis_inverse;};
private:
    void Load_State_From_Simulation_Folder(std::string folder_name,int rom_size);
//##################################################################### 
};
}
#endif
