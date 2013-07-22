//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION_PHI
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Images/EPS_FILE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_PHI.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/GENERALIZED_FLUID_MASS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_CUT.h>
using namespace PhysBAM;
//namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
//namespace PhysBAM{template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
FLUID_TO_SOLID_INTERPOLATION_PHI(const COLLISION_AWARE_INDEX_MAP<TV>& map,const ARRAY<T,TV_INT>& phi,SEGMENTED_CURVE_2D<T>& curve_input,T density_input)
    :BASE(map,curve_input,density_input),phi(phi),cut_order(4)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
~FLUID_TO_SOLID_INTERPOLATION_PHI()
{
}
//#####################################################################
// Function Setup_Mesh
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
Setup_Mesh()
{
    ARRAY_VIEW<TV> X=curve.particles.X;
    for(int i=1;i<=curve.mesh.elements.m;i++) X(i)=TV();

    GRID<TV> dual_grid(index_map.grid.Get_Regular_Grid());
    ARRAY<T,TV_INT> dual_phi(index_map.grid.Node_Indices(2));
    PHYSBAM_ASSERT(cut_order>=3 && cut_order<=4);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(dual_grid,1);it.Valid();it.Next())
        dual_phi(it.index)=Node_Average(it.index);

    T mx=0;
    curve.mesh.elements.Remove_All();
    HASHTABLE<TV_INT,ARRAY<int> > used_cells;
    HASHTABLE<FACE_INDEX<TV::m>,TV> HX;
    HASHTABLE<TV_INT,ARRAY<FACE_INDEX<TV::m> > > cut_faces;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(index_map.grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face=it.Full_Index();
        TV_INT a=TV_INT::Axis_Vector(3-it.Axis()),node0=it.index-a,node1=it.index,node2=it.index+a,node3=it.index+2*a;
        T phi1=dual_phi(node1),phi2=dual_phi(node2);
        if((phi1>0)==(phi2>0)) continue;
        T theta=0;
        if(cut_order==3) theta=Compute_Cut3(it.Full_Index(),phi1,phi2);
        else if(cut_order==4) theta=Compute_Cut4(it.Full_Index(),phi1,phi2);

        TV X1=dual_grid.X(node1),X2=dual_grid.X(node2);
        TV XX=(1-theta)*X1+theta*X2;
        mx=std::max(mx,abs((XX-index_map.grid.domain.Center()).Magnitude()-(T)0.01));

        HX.Set(it.Full_Index(),XX);
        cut_faces.Get_Or_Insert(face.First_Cell_Index()).Append(face);
        cut_faces.Get_Or_Insert(face.Second_Cell_Index()).Append(face);}
    std::stringstream ss;
    ss<<"interface linf error "<<mx<<std::endl;
    LOG::filecout(ss.str());


    HASHTABLE_ITERATOR<FACE_INDEX<TV::m>,TV> it(HX);
    if(!it.Valid()) return;
    FACE_INDEX<TV::m> prev_face(it.Key());
    TV_INT node1=prev_face.index,node2=prev_face.index+TV_INT::Axis_Vector(3-prev_face.axis);
    T phi1=dual_phi(node1),phi2=dual_phi(node2);
    int first_side=phi1<phi2?prev_face.axis:3-prev_face.axis;
    int next=0;
    X(++next)=HX.Get(prev_face);
    TV_INT cell=prev_face.Cell_Index(first_side);
    typename BASE::CLIP_ENTRY ce={0,0,1};
    T min_length=index_map.grid.dX.Min()*(T).1;
    for(int i=1;;i++){
        const ARRAY<FACE_INDEX<TV::m> >& array=cut_faces.Get(cell);
        if(array.m!=2 && !index_map.grid.Domain_Indices().Lazy_Inside(cell)) PHYSBAM_FATAL_ERROR("Level set extends outside of domain");
        PHYSBAM_ASSERT(array.m==2);
        FACE_INDEX<TV::m> next_face=(array(1)==prev_face)?array(2):array(1);
        if(next_face==it.Key()) break;
        TV next_X=HX.Get(next_face);
        if((next_X-X(next)).Magnitude()>min_length){
            X(++next)=next_X;
            ce.i=curve.mesh.elements.Append(VECTOR<int,2>(next-1,next));
            cut_cells.Get_Or_Insert(cell).clipped_segments.Append(ce);}
        else {std::stringstream ss1;ss1<<"PRUNE SEGMENT  "<<(next_X-X(next)).Magnitude()<<"  "<<min_length<<"   "<<next<<std::endl;LOG::filecout(ss.str());}
        cell=(cell==next_face.First_Cell_Index())?next_face.Second_Cell_Index():next_face.First_Cell_Index();
        prev_face=next_face;}
    if((X(1)-X(next)).Magnitude()>min_length){
        ce.i=curve.mesh.elements.Append(VECTOR<int,2>(next,1));
        cut_cells.Get_Or_Insert(cell).clipped_segments.Append(ce);}
    else{
        curve.mesh.elements.Last().y=1;
        std::stringstream ss1;
        ss1<<"PRUNE SEGMENT  "<<(X(1)-X(next)).Magnitude()<<"  "<<min_length<<"   "<<next<<std::endl;LOG::filecout(ss.str());}

    Remove_Degeneracy();
}
//#####################################################################
// Function Linear_Average
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
Linear_Average(const FACE_INDEX<TV::m>& face)
{
    TV_INT cell=face.index;
    T phi2=phi(cell);
    cell(face.axis)--;
    T phi1=phi(cell);
    cell(face.axis)--;
    T phi0=phi(cell);
    cell(face.axis)+=3;
    T phi3=phi(cell);
    return ((T)9/16)*(phi1+phi2)-((T)1/16)*(phi0+phi3);
}
//#####################################################################
// Function Node_Average
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
Node_Average(const TV_INT& cell)
{
    T c=phi(cell)+phi(cell-TV_INT(1,0))+phi(cell-TV_INT(1,1))+phi(cell-TV_INT(0,1));
    T h=phi(cell-TV_INT(2,0))+phi(cell-TV_INT(2,1))+phi(cell-TV_INT(-1,0))+phi(cell-TV_INT(-1,1));
    T v=phi(cell-TV_INT(0,2))+phi(cell-TV_INT(1,2))+phi(cell-TV_INT(0,-1))+phi(cell-TV_INT(1,-1));
    return ((T)5/16)*c-((T)1/32)*(h+v);
}
//#####################################################################
// Function Compute_Cut3
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
Compute_Cut3(const FACE_INDEX<TV::m>& face,T phi0,T phi1)
{
    T phia=Linear_Average(face);
    QUADRATIC<T> q(0,0,0);
    T mxabs=maxabs(phi0,phia,phi1);
    phi0/=mxabs;
    phia/=mxabs;
    phi1/=mxabs;
    q.Coefficients_From_Interpolation(0,phi0,(T).5,phia,1,phi1);
    return ITERATIVE_SOLVER<T>().Bisection_Secant_Root(q,0,1);
}
//#####################################################################
// Function Compute_Cut4
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
Compute_Cut4(const FACE_INDEX<TV::m>& face,T phi0,T phi1)
{
    FACE_INDEX<TV::m> tface=face;
    tface.index(3-tface.axis)--;
    T phia=Linear_Average(tface);
    tface.index(3-tface.axis)+=2;
    T phib=Linear_Average(tface);
    CUBIC<T> c(0,0,0,0);
    T mxabs=maxabs(phia,phi0,phi1,phib);
    phia/=mxabs;
    phi0/=mxabs;
    phi1/=mxabs;
    phib/=mxabs;
    c.Coefficients_From_Interpolation(-(T).5,phia,0,phi0,1,phi1,(T)1.5,phib);
    return ITERATIVE_SOLVER<T>().Bisection_Secant_Root(c,0,1);
}
//#####################################################################
// Function Setup_Before_Compute
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_PHI<TV>::
Setup_Before_Compute(ARRAY<bool,TV_INT>& outside_fluid_input,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_input)
{
    outside_fluid=&outside_fluid_input;
    psi_N=&psi_N_input;

    for(typename HASHTABLE<TV_INT,typename BASE::CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next())
        (*outside_fluid)(it.Key())=false;

    const_cast<COLLISION_AWARE_INDEX_MAP<TV>&>(index_map).number_extra_cells=cut_cells.Size();
}
//#####################################################################
template class FLUID_TO_SOLID_INTERPOLATION_PHI<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FLUID_TO_SOLID_INTERPOLATION_PHI<VECTOR<double,2> >;
#endif
