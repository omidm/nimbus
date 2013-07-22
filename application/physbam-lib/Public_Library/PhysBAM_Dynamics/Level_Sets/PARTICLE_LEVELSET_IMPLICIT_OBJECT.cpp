//#####################################################################
// Copyright 2006-2007, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>::
PARTICLE_LEVELSET_IMPLICIT_OBJECT(T_PARTICLE_LEVELSET& particle_levelset_input)
    :LEVELSET_IMPLICIT_OBJECT<TV>(particle_levelset_input.levelset.grid,particle_levelset_input.levelset.phi),particle_levelset(particle_levelset_input)
{
    particle_influence.Resize(levelset.grid.Node_Indices(1));p_grid=levelset.grid.Get_Regular_Grid();
}
//#####################################################################
// Create
//#####################################################################
template<class TV> PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>* PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>::
Create()
{
    int ghost_cells=3;
    PARTICLE_LEVELSET_IMPLICIT_OBJECT* levelset_implicit_object=new PARTICLE_LEVELSET_IMPLICIT_OBJECT(*(new T_PARTICLE_LEVELSET(*(new GRID<TV>),*(new T_ARRAYS_SCALAR),ghost_cells)));
    levelset_implicit_object->need_destroy_data=true;return levelset_implicit_object;
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template PARTICLE_LEVELSET_IMPLICIT_OBJECT<VECTOR<T,d> >::PARTICLE_LEVELSET_IMPLICIT_OBJECT(LEVELSET_POLICY<GRID<VECTOR<T,d> > >::PARTICLE_LEVELSET&); \
    template PARTICLE_LEVELSET_IMPLICIT_OBJECT<VECTOR<T,d> >* PARTICLE_LEVELSET_IMPLICIT_OBJECT<VECTOR<T,d> >::Create();
INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif
