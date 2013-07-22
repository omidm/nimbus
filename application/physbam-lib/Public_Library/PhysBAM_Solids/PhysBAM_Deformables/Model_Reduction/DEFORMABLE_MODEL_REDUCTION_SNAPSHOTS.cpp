//#####################################################################
// Copyright 2012, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Model_Reduction/DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>

using namespace PhysBAM;
//#####################################################################
// Function Load_State_From_Simulation_Folder
//#####################################################################
template <class TV> void DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<TV>::Load_State_From_Simulation_Folder(std::string folder_name,int rom_size)
{
    const int dimension = TV::dimension;
    const int doubleDimension = dimension*2;

    unsigned last_frame=0;
    FILE_UTILITIES::Read_From_Text_File(folder_name+"/common/last_frame",last_frame);

    const STREAM_TYPE stream_type((T()));
    Initialize_Geometry_Particle(); //register particle read
    DEFORMABLE_GEOMETRY_COLLECTION<TV> deformable_geometry_collection(*(new GEOMETRY_PARTICLES<TV>()));
    deformable_geometry_collection.Read(stream_type,folder_name,0,0,false);
    
    MATRIX_MXN<T> state_matrix;
    //store position and velocity for each particle in a state vector, one for each frame
    state_matrix.Resize(deformable_geometry_collection.particles.X.m*doubleDimension,last_frame);

    //take the SVD of the state matrix
    unsigned frame; 
    for(frame=0;frame<last_frame;frame++){
        deformable_geometry_collection.Read(stream_type,folder_name,frame,frame,false);
        for(int i=1;i<=deformable_geometry_collection.particles.X.m*doubleDimension;i+=doubleDimension)
        {
            state_matrix(i,frame+1)=deformable_geometry_collection.particles.X((i/doubleDimension)+1)(1);
            state_matrix(i+1,frame+1)=deformable_geometry_collection.particles.X((i/doubleDimension)+1)(2);
            if(dimension==2)
            {
                state_matrix(i+2,frame+1)=deformable_geometry_collection.particles.V((i/doubleDimension)+1)(1);
                state_matrix(i+3,frame+1)=deformable_geometry_collection.particles.V((i/doubleDimension)+1)(2);
            }
            else if(dimension==3)
            {
                state_matrix(i+2,frame+1)=deformable_geometry_collection.particles.X((i/doubleDimension)+1)(3);
                state_matrix(i+3,frame+1)=deformable_geometry_collection.particles.V((i/doubleDimension)+1)(1);
                state_matrix(i+4,frame+1)=deformable_geometry_collection.particles.V((i/doubleDimension)+1)(2);
                state_matrix(i+5,frame+1)=deformable_geometry_collection.particles.V((i/doubleDimension)+1)(3);
            }    
        }
    }
    MATRIX_MXN<T> U,E,V;
    U.Resize(deformable_geometry_collection.particles.X.m*doubleDimension,deformable_geometry_collection.particles.X.m*doubleDimension);
    E.Resize(deformable_geometry_collection.particles.X.m*doubleDimension,last_frame);
    V.Resize(last_frame,last_frame);
    state_matrix.Singular_Value_Decomposition(U,E,V);

    //estimate rank
    VECTOR_ND<T> E_diag;E.Get_Diagonal(E_diag);
    T rank_tolerance=(T)1e-12;//what we will consider to be a threshold for a zero singular value
    int rank=0;
    while(E_diag(rank+1)>rank_tolerance){rank++;if(rank==E_diag.Size()) break;}
    
    //extract phi
    MATRIX_MXN<T> phi;phi.Resize(U.Rows(),rank);
    for(int i=1;i<=rank;i++){for(int j=1;j<=U.Rows();j++){phi(j,i)=U(j,i);}}
    
    //find reduced basis
    if(rom_size>rank){PHYSBAM_WARNING("rom_size specified was greater than rank, using rom_size=rank instead");
        rom_size=rank;}
    MATRIX_MXN<T> reduction_matrix;reduction_matrix.Resize(rank,rom_size);
    reduction_matrix.Set_Zero_Matrix();
    for(int i=1;i<=rom_size;i++){reduction_matrix(i,i)=1;}
    reduced_basis=phi*reduction_matrix;

    FILE_UTILITIES::Write_To_Text_File(folder_name+"/common/reduced_basis_"+STRING_UTILITIES::string_sprintf("%d",rom_size),reduced_basis);
}
//#####################################################################
// Function Load_Reduced_Basis_From_Simulation_Folder
//#####################################################################
template <class TV> void DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<TV>::
Load_Reduced_Basis_From_Simulation_Folder(std::string folder_name,int rom_size)
{
    const int doubleDimension = TV::dimension*2;

    const STREAM_TYPE stream_type((T()));
    DEFORMABLE_GEOMETRY_COLLECTION<TV> deformable_geometry_collection(*(new GEOMETRY_PARTICLES<TV>()));
    deformable_geometry_collection.Read(stream_type,folder_name,0,-1,false);

    int rank = 0;
    FILE_UTILITIES::Read_From_Text_File(folder_name+"/common/reduced_basis_rank",rank);

    int numRows = 0;
    FILE_UTILITIES::Read_From_Text_File(folder_name+"/common/reduced_basis_size",numRows);

    //find reduced basis
    if(rom_size>rank){PHYSBAM_WARNING("rom_size specified was greater than rank, using rom_size=rank instead");
        rom_size=rank;}
    
    MATRIX_MXN<T> phi;
    phi.Resize(numRows,rank);
    FILE_UTILITIES::Read_From_Text_File(folder_name+"/common/reduced_basis",phi);

    MATRIX_MXN<T> reduction_matrix;reduction_matrix.Resize(rank,rom_size);
    reduction_matrix.Set_Zero_Matrix();
    for(int i=1;i<=rom_size;i++){reduction_matrix(i,i)=1;}
    reduced_basis.Resize(deformable_geometry_collection.particles.X.m*doubleDimension,rom_size);
    reduced_basis=phi*reduction_matrix;    
}
//#####################################################################
// Function Generate_Reduced_Basis
//#####################################################################
template <class TV> void DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<TV>::
Generate_Reduced_Basis(std::string folder_name,int rom_size)
{
    //function defined in case we want to use multiple folders later
    Load_State_From_Simulation_Folder(folder_name,rom_size);
}
//#####################################################################
// Function Convert_Reduced_To_Full
//#####################################################################
template <class TV> void DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<TV>::
Convert_Reduced_To_Full(VECTOR_ND<T>& v_reduced, STATE_VECTOR<TV>& v_full)
{
    //V is orthogonal here, so W = V
    //might need W for non orthogonal matrices
    v_full=reduced_basis*v_reduced;
}
//#####################################################################
// Function Convert_Full_To_Reduced
//#####################################################################
template <class TV> void DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<TV>::
Convert_Full_To_Reduced(VECTOR_ND<T>& v_reduced, STATE_VECTOR<TV>& v_full)
{
    if(orthonormal_basis)
    {
        //V is orthogonal here, so W = V
        //v_reduced=reduced_basis.Transpose_Times(v_full);
        for(int i=1;i<=reduced_basis.Columns();i++)
        {
            v_reduced(i)=0;
            for(int j=1;j<=reduced_basis.Rows();j++)
            {
                v_reduced(i)+=reduced_basis(j,i)*v_full(j);
            }
        }
    }
    else
    {
        for(int i=1;i<=reduced_basis_inverse.Rows();i++)
        {
            v_reduced(i)=0;
            for(int j=1;j<=reduced_basis_inverse.Columns();j++)
            {
                v_reduced(i)+=reduced_basis_inverse(i,j)*v_full(j);
            }
        }
    }
}
//#####################################################################
// Function Set_Reduced_Basis
//#####################################################################
template <class TV> void DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<TV>::
Set_Reduced_Basis(MATRIX_MXN<T>& m_in)
{
    reduced_basis=m_in;
    if(!orthonormal_basis) reduced_basis_inverse=m_in.Inverse();
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<VECTOR<T,2> >;  \
    template class DEFORMABLE_MODEL_REDUCTION_SNAPSHOTS<VECTOR<T,3> >;
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
