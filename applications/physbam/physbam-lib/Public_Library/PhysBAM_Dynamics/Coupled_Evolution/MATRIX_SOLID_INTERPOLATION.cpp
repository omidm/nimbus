//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_INTERPOLATION
//##################################################################### 
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION<TV>::
MATRIX_SOLID_INTERPOLATION(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info)
    :MATRIX_SOLID_INTERPOLATION_BASE<TV>(info)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION<TV>::
~MATRIX_SOLID_INTERPOLATION()
{
}
//#####################################################################
// Function Number_Of_Constraints
//#####################################################################
template<class TV> COUPLING_CONSTRAINT_ID MATRIX_SOLID_INTERPOLATION<TV>::
Number_Of_Constraints() const
{
    return rows.m;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION<TV>::
Compute(const int ghost_cells)
{
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_SIMPLEX;typedef VECTOR<int,TV::dimension> TV_INT;
    // TODO: save ghost cells if no MPI?
    ARRAY<T_SIMPLEX> clipped_simplices;
    typename GRID<TV>::REGION region_type=ghost_cells?GRID<TV>::INTERIOR_REGION:GRID<TV>::WHOLE_REGION;
    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(iterator_info,ghost_cells,region_type);iterator.Valid();iterator.Next()){
        const ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >& simplices=iterator.Get_Simplices();
        TV accumulated_normal;
        T accumulated_flux=0;
        RANGE<TV> dual_cell(iterator.Dual_Cell().Thickened((T)2*iterator_info.iterator_rasterization_thickness));
        ROW& row=rows(rows.Append(ROW()));

        for(int s=1;s<=simplices.m;s++){
            if(!simplices(s).y) PHYSBAM_NOT_IMPLEMENTED();

            COLLISION_GEOMETRY<TV>& body=*iterator_info.coupling_bodies(simplices(s).x);
            const T_SIMPLEX& simplex=body.World_Space_Simplex(simplices(s).y);
            TV normal=simplex.Normal();
            T cosine=normal(iterator.axis),area=0;
            if(cosine<0){normal=-normal;cosine=-cosine;}

            simplex.Clip_To_Box(dual_cell,clipped_simplices);
            if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* deformable=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(&body)){
                TV simplex_weight;
                for(int c=1;c<=clipped_simplices.m;c++)
                    simplex_weight+=simplex.Sum_Barycentric_Coordinates(clipped_simplices(c))*clipped_simplices(c).Size();
                simplex_weight/=TV::dimension;
                area=simplex_weight.Sum();
                const TV_INT& element=deformable->object.mesh.elements(simplices(s).y);
                for(int i=1;i<=element.m;i++)
                    row.deformable_weights.Append(DEFORMABLE_WEIGHT(element(i),cosine*simplex_weight(i)));}
            else if(RIGID_COLLISION_GEOMETRY<TV>* rigid_body_wrapper=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(&body)){
                for(int i=1;i<=clipped_simplices.m;i++) area+=clipped_simplices(i).Size();
                row.rigid_weights.Append(RIGID_WEIGHT(rigid_body_wrapper->rigid_geometry.particle_index,cosine*area,iterator.Location()-dynamic_cast<RIGID_BODY<TV>&>(rigid_body_wrapper->rigid_geometry).X()));}
            else PHYSBAM_FATAL_ERROR("OMGWTFBBQ");

            accumulated_normal+=normal*area;
            accumulated_flux+=cosine*area;}

        row.normal=accumulated_normal.Normalized();
        row.axis=iterator.Axis();
        if(accumulated_flux){
            T factor=1/accumulated_flux;
            for(int i=1;i<=row.deformable_weights.m;i++) row.deformable_weights(i).weight*=factor;
            for(int i=1;i<=row.rigid_weights.m;i++) row.rigid_weights(i).weight*=factor;}
        else{
            T weight=(T)1/(row.deformable_weights.m+row.rigid_weights.m);
            for(int i=1;i<=row.deformable_weights.m;i++) row.deformable_weights(i).weight=weight;
            for(int i=1;i<=row.rigid_weights.m;i++) row.rigid_weights(i).weight=weight;}

        row.deformable_weights.Coalesce();
        row.rigid_weights.Coalesce();}
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION<TV>::
Times_Add(const GENERALIZED_VELOCITY<TV>& solids,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++){
        const ARRAY<DEFORMABLE_WEIGHT>& deformable_weights=rows(i).deformable_weights;
        const ARRAY<RIGID_WEIGHT>& rigid_weights=rows(i).rigid_weights;
        const int axis=rows(i).axis;

        T weighted_velocity_component=0;
        for(int j=1;j<=deformable_weights.m;j++)
            weighted_velocity_component+=solids.V.array(deformable_weights(j).index)(axis)*deformable_weights(j).weight;
        for(int j=1;j<=rigid_weights.m;j++){
            const TWIST<TV>& twist=solids.rigid_V.array(rigid_weights(j).index);
            weighted_velocity_component+=(twist.linear+TV::Cross_Product(twist.angular,rigid_weights(j).radius))(axis)*rigid_weights(j).weight;}
        constraints(i)+=weighted_velocity_component;}
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION<TV>::
Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,GENERALIZED_VELOCITY<TV>& solids) const
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++){
        const ARRAY<DEFORMABLE_WEIGHT>& deformable_weights=rows(i).deformable_weights;
        const ARRAY<RIGID_WEIGHT>& rigid_weights=rows(i).rigid_weights;
        const int axis=rows(i).axis;
        for(int j=1;j<=deformable_weights.m;j++)
            solids.V.array(deformable_weights(j).index)(axis)+=constraints(i)*deformable_weights(j).weight;
        for(int j=1;j<=rigid_weights.m;j++){
            TWIST<TV>& twist=solids.rigid_V.array(rigid_weights(j).index);
            T linear_contribution=constraints(i)*rigid_weights(j).weight;
            twist.linear(axis)+=linear_contribution;
            twist.angular+=TV::Cross_Product(rigid_weights(j).radius,TV::Axis_Vector(axis))*linear_contribution;}}
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION<TV>::
Print_Each_Matrix(int n,GENERALIZED_VELOCITY<TV>& G) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("J-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("J",Value(rows.m),G.Raw_Size());
    ARRAY<int> reverse_map_deformable(G.V.array.Size());
    reverse_map_deformable.Subset(G.V.indices)=IDENTITY_ARRAY<>(G.V.Size());
    ARRAY<int> reverse_map_rigid(G.rigid_V.array.Size());
    reverse_map_rigid.Subset(G.rigid_V.indices)=IDENTITY_ARRAY<>(G.rigid_V.Size());

    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++){
        const ARRAY<DEFORMABLE_WEIGHT>& deformable_weights=rows(i).deformable_weights;
        const ARRAY<RIGID_WEIGHT>& rigid_weights=rows(i).rigid_weights;
        const int axis=rows(i).axis;

        for(int j=1;j<=deformable_weights.m;j++)
            oo.Add_Sparse_Entry(Value(i),(reverse_map_deformable(deformable_weights(j).index)-1)*TV::dimension+axis,deformable_weights(j).weight);

        for(int j=1;j<=rigid_weights.m;j++) if(int index=reverse_map_rigid(rigid_weights(j).index)){ // Prune out static/kinematic
            int base=G.V.Size()*TV::dimension+(index-1)*TWIST<TV>::dimension;
            oo.Add_Sparse_Entry(Value(i),base+axis,rigid_weights(j).weight);
            VECTOR<T,TV::SPIN::m> cpm=MATRIX<T,TV::SPIN::m,TV::m>::Cross_Product_Matrix(rigid_weights(j).radius).Column(axis)*rigid_weights(j).weight;
            for(int k=1;k<=TV::SPIN::m;k++) oo.Add_Sparse_Entry(Value(i),base+TV::m+k,cpm(k));}}

    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Add_Diagonal
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION<TV>::
Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_MASS<TV>& solid_mass) const
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++){
        const ARRAY<DEFORMABLE_WEIGHT>& deformable_weights=rows(i).deformable_weights;
        const ARRAY<RIGID_WEIGHT>& rigid_weights=rows(i).rigid_weights;
        T diag=0;
        for(int j=1;j<=deformable_weights.m;j++){
            diag+=sqr(deformable_weights(j).weight)*solid_mass.one_over_mass.array(deformable_weights(j).index);}
        for(int j=1;j<=rigid_weights.m;j++){
            const RIGID_BODY_MASS<TV,true>& M=solid_mass.world_space_rigid_mass_inverse.array(rigid_weights(j).index);
            typename TV::SPIN r=TV::Cross_Product(rigid_weights(j).radius,TV::Axis_Vector(rows(i).axis));
            diag+=sqr(rigid_weights(j).weight)*(M.mass+TV::SPIN::Dot_Product(M.inertia_tensor*r,r));}
        diagonal(i)+=diag;}
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    ARRAY<int> reverse_map_deformable(this->V_size);
    reverse_map_deformable.Subset(*this->V_indices)=IDENTITY_ARRAY<>(this->V_indices->m);
    ARRAY<int> reverse_map_rigid(this->rigid_V_size);
    reverse_map_rigid.Subset(*this->rigid_V_indices)=IDENTITY_ARRAY<>(this->rigid_V_indices->m);

    for(COUPLING_CONSTRAINT_ID i(1);i<=rows.m;i++){
        const ARRAY<DEFORMABLE_WEIGHT>& deformable_weights=rows(i).deformable_weights;
        const ARRAY<RIGID_WEIGHT>& rigid_weights=rows(i).rigid_weights;
        const int axis=rows(i).axis;

        for(int j=1;j<=deformable_weights.m;j++)
            data.Append(TRIPLE<int,int,T>(Value(i),(reverse_map_deformable(deformable_weights(j).index)-1)*TV::dimension+axis,deformable_weights(j).weight));

        for(int j=1;j<=rigid_weights.m;j++) if(int index=reverse_map_rigid(rigid_weights(j).index)){ // Prune out static/kinematic
            int base=this->V_indices->m*TV::dimension+(index-1)*TWIST<TV>::dimension;
            data.Append(TRIPLE<int,int,T>(Value(i),base+axis,rigid_weights(j).weight));
            VECTOR<T,TV::SPIN::m> cpm=MATRIX<T,TV::SPIN::m,TV::m>::Cross_Product_Matrix(rigid_weights(j).radius).Column(axis)*rigid_weights(j).weight;
            for(int k=1;k<=TV::SPIN::m;k++) data.Append(TRIPLE<int,int,T>(Value(i),base+TV::m+k,cpm(k)));}}
}
//#####################################################################
template class MATRIX_SOLID_INTERPOLATION<VECTOR<float,1> >;
template class MATRIX_SOLID_INTERPOLATION<VECTOR<float,2> >;
template class MATRIX_SOLID_INTERPOLATION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_SOLID_INTERPOLATION<VECTOR<double,1> >;
template class MATRIX_SOLID_INTERPOLATION<VECTOR<double,2> >;
template class MATRIX_SOLID_INTERPOLATION<VECTOR<double,3> >;
#endif
