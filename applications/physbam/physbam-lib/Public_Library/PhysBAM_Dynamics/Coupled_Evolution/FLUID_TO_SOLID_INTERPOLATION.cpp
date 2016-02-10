//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION<TV>::
FLUID_TO_SOLID_INTERPOLATION(const COLLISION_AWARE_INDEX_MAP<TV>& map,const PARTICLES<TV>& particles_input)
    :FLUID_TO_SOLID_INTERPOLATION_BASE<TV>(map),max_dist(2),particles(particles_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION<TV>::
~FLUID_TO_SOLID_INTERPOLATION()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Compute_Weights(const TV& X,int axis,ARRAY<ENTRY>& array)
{
    T limit=max_dist*index_map.grid.dX.Max(),total_weight=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(index_map.grid,RANGE<TV_INT>(index_map.grid.Index(X)).Thickened(3),axis);it.Valid();it.Next()){
        T dist=(it.Location()-X).Magnitude();
        if(dist>limit) continue;
        int index=index_map.face_indices(it.Full_Index());
        if(!index) continue;
        ENTRY e={LEVELSET_UTILITIES<T>::Delta(dist,limit),index};
        total_weight+=e.w;
        array.Append(e);}
    PHYSBAM_ASSERT(total_weight>0);
    for(int i=1;i<=array.m;i++) array(i).w/=total_weight;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Compute(const int ghost_cells)
{
    entries.Resize(coupled_particles.m);
    for(int i=1;i<=entries.m;i++)
        for(int a=1;a<=TV::m;a++)
            Compute_Weights(particles.X(coupled_particles(i)),a,entries(i)(a));
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Times_Add(const VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const
{
    for(int j=1;j<=coupled_particles.m;j++){int p=coupled_particles(j);
        for(int a=1;a<=TV::m;a++){
            T& v=solid_velocity.V.array(p)(a);
            const ARRAY<ENTRY>& array=entries(p)(a);
            for(int i=1;i<=array.m;i++)
                v+=array(i).w*fluid_velocity(array(i).i);}}
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Transpose_Times_Add(const GENERALIZED_VELOCITY<TV>& solid_force,VECTOR_ND<T>& fluid_force) const
{
    for(int j=1;j<=coupled_particles.m;j++){int p=coupled_particles(j);
        for(int a=1;a<=TV::m;a++){
            T v=solid_force.V.array(p)(a);
            const ARRAY<ENTRY>& array=entries(p)(a);
            for(int i=1;i<=array.m;i++)
                fluid_force(array(i).i)+=array(i).w*v;}}
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Print_Each_Matrix(int n,int fluid_faces,GENERALIZED_VELOCITY<TV>& G) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("H-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("H",G.Raw_Size(),fluid_faces);
    ARRAY<int> reverse_map_deformable(G.V.array.Size());
    reverse_map_deformable.Subset(G.V.indices)=IDENTITY_ARRAY<>(G.V.Size());

    for(int j=1;j<=coupled_particles.m;j++){int p=coupled_particles(j);
        for(int a=1;a<=TV::m;a++){
            const ARRAY<ENTRY>& array=entries(p)(a);
            for(int i=1;i<=array.m;i++)
                oo.Add_Sparse_Entry((reverse_map_deformable(p)-1)*TV::m+a,array(i).i,array(i).w);}}

    oo.End_Sparse_Matrix();
}
//#####################################################################
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<float,1> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<float,2> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<double,1> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<double,2> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<double,3> >;
#endif
