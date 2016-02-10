//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM_SYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_VISCOSITY_UNIFORM<TV>::
LEVELSET_VISCOSITY_UNIFORM(BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input,const GRID<TV>& grid_input,T dt,T density,T viscosity)
    :index_map(grid_input,callback_input),system(index_map),scale(dt*viscosity/density),print_matrix(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_VISCOSITY_UNIFORM<TV>::
~LEVELSET_VISCOSITY_UNIFORM()
{
}
//#####################################################################
// Function Apply_Full_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Full_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,bool fully_explicit,bool fully_implicit,int axis)
{
    T local_scale=(fully_explicit || fully_implicit)?scale:scale/2;
    index_map.Compute(axis,periodic_boundary);
    std::stringstream ss;ss<<"Solve size "<<index_map.index_to_face.m<<std::endl;LOG::filecout(ss.str());
    system.Compute(axis,local_scale);
    Resize_Vectors(fully_explicit);
    index_map.Gather(u,b.v);

    static int solve_id=0;solve_id++;
    if(print_matrix){
        ss.clear();ss<<"viscosity solve id "<<solve_id<<std::endl;LOG::filecout(ss.str());
        if(fully_explicit){q.v.Resize(index_map.index_to_face.m);s.v.Resize(index_map.index_to_face.m);}
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-M-%i.txt",solve_id).c_str()).Write("M",system,q,s);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-b-%i.txt",solve_id).c_str()).Write("b",b);}

    if(!fully_implicit) Apply_Explicit_Viscosity(u,axis);
    if(!fully_implicit && !fully_explicit) exchange(x.v.x,b.v.x);
    if(!fully_explicit) Apply_Implicit_Viscosity(u,axis);
    if(print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("visc-x-%i.txt",solve_id).c_str()).Write("x",x);
    index_map.Scatter(x.v,u);
}
//#####################################################################
// Function Apply_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,bool fully_explicit,bool fully_implicit,bool coupled)
{
    if(coupled) Apply_Full_Viscosity(u,fully_explicit,fully_implicit,0);
    else for(int a=1;a<=d;a++) Apply_Full_Viscosity(u,fully_explicit,fully_implicit,a);
}
//#####################################################################
// Function Apply_Implicit_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Implicit_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,int axis)
{
    system.Add_Constant_Part(b);
    x.v=b.v;
    CONJUGATE_GRADIENT<T> cg;
    cg.print_diagnostics=true;
    bool result=cg.Solve(system,x,b,q,s,r,k,z,(T)1e-4,1,1000);
    PHYSBAM_ASSERT(result);
}
//#####################################################################
// Function Apply_Explicit_Viscosity
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Apply_Explicit_Viscosity(ARRAY<T,FACE_INDEX<d> >& u,int axis)
{
    system.poisson.P.Times(x.v,b.v);
    system.Add_Constant_Part(x);
    b.v=x.v+scale*b.v;
    x.v=b.v;
}
//#####################################################################
// Function Resize_Vectors
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM<TV>::
Resize_Vectors(bool minimal)
{
    x.v.Resize(index_map.index_to_face.m);
    b.v.Resize(index_map.index_to_face.m);
    if(!minimal){
        q.v.Resize(index_map.index_to_face.m);
        s.v.Resize(index_map.index_to_face.m);
        r.v.Resize(index_map.index_to_face.m);
        k.v.Resize(index_map.index_to_face.m);
        z.v.Resize(index_map.index_to_face.m);}
}
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<float,1> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<float,2> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<double,1> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<double,2> >;
template class LEVELSET_VISCOSITY_UNIFORM<VECTOR<double,3> >;
#endif
