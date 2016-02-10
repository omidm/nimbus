//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ron Fedkiw, Eran Guendelman, Don Hatch, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/FACTORIAL.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_Properties
//#####################################################################
namespace{
template<class T,int d> inline MATRIX<T,1> Scaled_Element_Covariance(const T scaled_element_volume,const MATRIX<T,2,d>& DX) // actually returns only the trace of the covariance matrix
{
    static const SYMMETRIC_MATRIX<T,d> canonical=(T)1+SYMMETRIC_MATRIX<T,d>::Unit_Matrix();
    return MATRIX<T,1>(scaled_element_volume*MATRIX<T,2,d>::Inner_Product(DX*canonical,DX));
}
template<class T,int d> inline SYMMETRIC_MATRIX<T,3> Scaled_Element_Covariance(const T scaled_element_volume,const MATRIX<T,3,d>& DX)
{
    static const SYMMETRIC_MATRIX<T,d> canonical=(T)1+SYMMETRIC_MATRIX<T,d>::Unit_Matrix();
    return SYMMETRIC_MATRIX<T,3>::Conjugate(DX,scaled_element_volume*canonical);
}
template<class T> inline T Inertia_Tensor_From_Covariance(const T covariance_trace)
{
    return covariance_trace;
}
template<class T> inline SYMMETRIC_MATRIX<T,3> Inertia_Tensor_From_Covariance(const SYMMETRIC_MATRIX<T,3>& covariance)
{
    return covariance.Trace()-covariance;
}}
template<class TV,int d> template<bool thin_shell> void MASS_PROPERTIES<TV,d>::
Compute_Properties(MASS_PROPERTIES<TV,d>&,NORMAL_IMPLEMENTATION)
{
    typedef typename TV::SCALAR T;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;

    static const bool filled=!thin_shell;
    if(!object.mesh.elements.m) PHYSBAM_FATAL_ERROR("mesh has no elements");
    if(filled && d!=TV::m-1) PHYSBAM_FATAL_ERROR("only codimension 1 objects can be filled");

    // compute volume and center
    const TV base=object.particles.X(object.mesh.elements(1)[1]);
    T scaled_volume=0; // (d+filled)!*volume
    TV scaled_center_times_volume; // (d+1+filled)!*center*volume
    TV particle_X[d+1];
    for(int t=1;t<=object.mesh.elements.m;t++){const VECTOR<int,d+1>& nodes=object.mesh.elements(t);
        MATRIX<T,TV::m,d+1> DX;
        for(int i=1;i<=nodes.m;i++) particle_X[i-1]=object.particles.X(nodes(i)); // THIS NEEDS TO BE HERE BECAUSE OF A COMPILER BUG IN GCC 4.0.2
        for(int i=1;i<=nodes.m;i++) DX.Column(i)=particle_X[i-1]-base;
        T scaled_element_volume=filled?DX.Parallelepiped_Measure():STRAIN_MEASURE<TV,d>::Ds(object.particles.X,nodes).Parallelepiped_Measure();
        if(t==-1){{std::stringstream ss;ss<<object.particles.X(nodes(1))<<object.particles.X(nodes(2))<<std::endl;LOG::filecout(ss.str());} // THIS NEEDS TO BE HERE BECAUSE OF A COMPILER BUG IN GCC 4.0.1
            {std::stringstream ss;ss<<object.particles.X(nodes(1))<<object.particles.X(nodes(3))<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<object.particles.X(nodes(3))<<base<<std::endl;LOG::filecout(ss.str());}}
        scaled_volume+=scaled_element_volume;
        scaled_center_times_volume+=scaled_element_volume*DX.Column_Sum();}
    volume=(T)1/FACTORIAL<d+filled>::value*scaled_volume;
    if(!volume) PHYSBAM_FATAL_ERROR("zero volume");
    center=base+(T)1/(d+1+filled)*scaled_center_times_volume/scaled_volume;

    // compute inertia tensor (see http://number-none.com/blow/inertia for explanation of filled case)
    T_WORLD_SPACE_INERTIA_TENSOR scaled_covariance=T_WORLD_SPACE_INERTIA_TENSOR(); // (d+2+filled)!*covariance (or trace(covariance) in 2d)
    for(int t=1;t<=object.mesh.elements.m;t++){const VECTOR<int,d+1>& nodes=object.mesh.elements(t);
        MATRIX<T,TV::m,d+1> DX;
        for(int i=1;i<=nodes.m;i++) particle_X[i-1]=object.particles.X(nodes(i)); // THIS NEEDS TO BE HERE BECAUSE OF A COMPILER BUG IN GCC 4.0.2
        for(int i=1;i<=nodes.m;i++) DX.Column(i)=particle_X[i-1]-center;
        T scaled_element_volume=filled?DX.Parallelepiped_Measure():STRAIN_MEASURE<TV,d>::Ds(object.particles.X,nodes).Parallelepiped_Measure();
        scaled_covariance+=Scaled_Element_Covariance(scaled_element_volume,DX);}
    
    inertia_tensor_over_density=Inertia_Tensor_From_Covariance((T)1/FACTORIAL<d+2+filled>::value*scaled_covariance);
    {std::stringstream ss;ss<<"inertia tensor over density "<<inertia_tensor_over_density<<std::endl;LOG::filecout(ss.str());} // THIS NEEDS TO BE HERE BECAUSE OF A COMPILER BUG IN GCC 4.0.2 (Frank wants it to stay std::endl)
}
template<class TV,int d> template<bool thin_shell> void MASS_PROPERTIES<TV,d>::
Compute_Properties(MASS_PROPERTIES<TV,d>&,DUMMY_IMPLEMENTATION )
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> MASS_PROPERTIES<TV,d>::
MASS_PROPERTIES(const T_SIMPLICIAL_OBJECT& object,const bool thin_shell)
    :object(object),mass(1),density(1),use_mass(true)
{
    if(thin_shell) Compute_Properties<true>(*this,ACCESS_IMPLEMENTATION());else Compute_Properties<false>(*this,ACCESS_IMPLEMENTATION());
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T_SIMPLICIAL_OBJECT> typename T_SIMPLICIAL_OBJECT::SCALAR
Shell_Volume_Helper(const T_SIMPLICIAL_OBJECT& object)
{
    typedef typename T_SIMPLICIAL_OBJECT::SCALAR T;typedef typename T_SIMPLICIAL_OBJECT::VECTOR_T TV;
    static const int d=T_SIMPLICIAL_OBJECT::MESH::dimension;
    T scaled_volume=0; // d!*volume
    for(int t=1;t<=object.mesh.elements.m;t++){const VECTOR<int,d+1>& nodes=object.mesh.elements(t);
        MATRIX<T,TV::m,d> Ds=STRAIN_MEASURE<TV,d>::Ds(object.particles.X,nodes);
        scaled_volume+=Ds.Parallelepiped_Measure();}
    return (T)1/FACTORIAL<d>::value*scaled_volume;
}
template<class T> T
Shell_Volume_Helper(const POINT_SIMPLICES_1D<T>& object)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV,int d> typename TV::SCALAR MASS_PROPERTIES<TV,d>::
Thin_Shell_Volume(const T_SIMPLICIAL_OBJECT& object)
{
    if(!object.mesh.elements.m) PHYSBAM_FATAL_ERROR("mesh has no elements");
    return Shell_Volume_Helper(object);
}
//#####################################################################
// Function Transform_To_Object_Frame
//#####################################################################
namespace{
template<class T> ROTATION<VECTOR<T,1> > Diagonalize_Inertia_Tensor(MATRIX<T,0>& object_space_inertia_tensor,const MATRIX<T,0> inertia_tensor)
{
    return ROTATION<VECTOR<T,1> >();
}
template<class T> ROTATION<VECTOR<T,2> > Diagonalize_Inertia_Tensor(MATRIX<T,1>& object_space_inertia_tensor,const MATRIX<T,1>& inertia_tensor)
{
    object_space_inertia_tensor=inertia_tensor;
    return ROTATION<VECTOR<T,2> >();
}
template<class T> ROTATION<VECTOR<T,3> > Diagonalize_Inertia_Tensor(DIAGONAL_MATRIX<T,3>& object_space_inertia_tensor,const SYMMETRIC_MATRIX<T,3>& inertia_tensor)
{
    MATRIX<T,3> rotation;inertia_tensor.Solve_Eigenproblem(object_space_inertia_tensor,rotation);
    return ROTATION<VECTOR<T,3> >(rotation);
}}
template<class TV,int d> void MASS_PROPERTIES<TV,d>::
Transform_To_Object_Frame(FRAME<TV>& frame,T_INERTIA_TENSOR& object_space_inertia_tensor) const
{
    frame.t=center;
    frame.r=Diagonalize_Inertia_Tensor<T>(object_space_inertia_tensor,Inertia_Tensor());
}
//#####################################################################
// Function Transform_To_Object_Frame
//#####################################################################
template<class TV,int d> void MASS_PROPERTIES<TV,d>::
Transform_To_Object_Frame(FRAME<TV>& frame,T_INERTIA_TENSOR& object_space_inertia_tensor,POINT_CLOUD<TV>& point_cloud) const
{
    Transform_To_Object_Frame(frame,object_space_inertia_tensor);
    if(frame.r==ROTATION<TV>()) point_cloud.X-=frame.t;
    else{PHYSBAM_ASSERT(TV::m==3);point_cloud.X=frame.Inverse()*point_cloud.X;}
}
//#####################################################################
template class MASS_PROPERTIES<VECTOR<float,1>,0>;
template class MASS_PROPERTIES<VECTOR<float,2>,1>;
template class MASS_PROPERTIES<VECTOR<float,3>,1>;
template class MASS_PROPERTIES<VECTOR<float,3>,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MASS_PROPERTIES<VECTOR<double,1>,0>;
template class MASS_PROPERTIES<VECTOR<double,2>,1>;
template class MASS_PROPERTIES<VECTOR<double,3>,1>;
template class MASS_PROPERTIES<VECTOR<double,3>,2>;
#endif
