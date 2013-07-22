//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_PARTICLE_COUPLING
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/GRID_PARTICLE_COUPLING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_PARTICLE_COUPLING<TV>::
GRID_PARTICLE_COUPLING(GRID<TV>& grid_input,FLUID_PARTICLES<TV>& particles_input)
    :grid(grid_input),particles(particles_input),levelset(grid_input,phi),surface_tension_force(0),
    number_particle_per_cell(1<<TV::m),pic_ratio(0.5),viscosity(1),
    output_matrices(false),use_constant_interpolation(false),output_face_velocities(false),fix_volume_error(false),use_log_based_volume_correction(true),use_IC_preconditioner(false),feed_back_for_momentum_conservation(true),not_feed_back_to_surface_faces(true),use_simplified_mass(true),apply_pressure(true),recompute_mass(true),max_cg_iterations(1000)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GRID_PARTICLE_COUPLING<TV>::
~GRID_PARTICLE_COUPLING()
{}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Initialize_Grid()
{
    cell_weights.Resize(grid.Domain_Indices(1));
    phi.Resize(grid.Domain_Indices(1));
    pressure.Resize(grid.Domain_Indices(1));
    psi_D.Resize(grid.Domain_Indices(1));
    psi_N.Resize(grid,1);
    cell_indices.Resize(grid.Domain_Indices(1));
    face_indices.Resize(grid,1);
    face_velocities.Resize(grid,1);
    face_forces.Resize(grid,1);
    surface_cells.Resize(grid.Domain_Indices(1));
    surface_faces.Resize(grid,0);
    record_time=0;
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Solve(const T& dt)
{
    LOG::Time("Setup system");
    Compute_Grid_Gradient_Matrix();
    if(viscosity>0) Compute_Grid_Diffusion_Matrix(dt);
    if(surface_tension_force && surface_tension_force->apply_implicit_forces) surface_tension_force->Construct_Implicit_Matrix_First_Half(C);
    Compute_Full_Coupling_Matrix();
    LOG::Time("Solve system");
    //Project_Out_Null_Space_Of_Nullspace();
    Solve_Coupling_System(dt);
    record_time+=dt;
    if(pic_ratio>0 && record_time>=(T)1/24){
        Interpolate_From_Grid_To_Particles(dt);
        record_time-=(T)1/24;}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Update_Position_Based_State()
{
    Compute_Cell_Weights();
    Compute_Levelset();
    Setup_Boundary_Conditions();
    Build_Index_Mapping();
    if(use_constant_interpolation) Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Constant();
    else Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Linear();
    if(recompute_mass) Compute_Mass();
}
//#####################################################################
// Function Compute_Levelset_And_Mass
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Levelset()
{
    for(int i=1;i<=phi.array.m;i++) phi.array(i)=(T)(-cell_weights.array(i)+0.5);
    levelset.Fast_Marching_Method(0,grid.dX.Max()*3);
}
//#####################################################################
// Function Compute_Cell_Weights
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Cell_Weights()
{
    T average=1/number_particle_per_cell;
    T one_over_radius=grid.one_over_dX.Max();
    ARRAYS_COMPUTATIONS::Fill(cell_weights.array,(T)0);
    for(int i=1;i<=particles.Size();i++){
        const TV& p=particles.X(i);
        TV_INT cell0=grid.Cell(p,0);
        RANGE<TV_INT> range(cell0-1,cell0+1);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,range);iter.Valid();iter.Next())
            cell_weights(iter.Cell_Index())+=Linear_Kernel(p,iter.Location(),one_over_radius)*average;}
}
//#####################################################################
// Function Compute_Mass
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Mass()
{
    T full_mass=grid.Cell_Size()*density;
    face_mass.Resize(index_faces.m);
    particle_mass.Resize(particles.one_over_mass.m*TV::m);
    for(int i=1;i<=index_faces.m;i++){
        T_FACE_INDEX face=index_faces(i);
        TV_INT cell1=face.index,cell2=face.index;cell1(face.axis)--;
        T theta=1;
        if(psi_D(cell1)!=psi_D(cell2)) theta=LEVELSET_UTILITIES<T>::Theta(phi(cell1),phi(cell2));
        if(psi_D(cell1)) theta=1-theta;
        face_mass(i)=theta*full_mass;}
    HT.Times(face_mass,particle_mass);
}
//#####################################################################
// Function Setup_Boundary_Conditions
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Setup_Boundary_Conditions()
{
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,1);iter.Valid();iter.Next()){
            TV_INT cell=iter.Cell_Index();
            psi_D(cell)=(phi(cell)>0 || !grid.Inside_Domain(cell));}
    if(0) for(int i=1;i<=particles.Size();i++){
        const TV& p=particles.X(i);
        TV_INT cell=grid.Cell(p,0);
        psi_D(cell)=false;}
    for(int axis=1;axis<=TV::m;axis++) ARRAYS_COMPUTATIONS::Fill(psi_N.Component(axis).array,false);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iter(grid,0,GRID<TV>::BOUNDARY_REGION);iter.Valid();iter.Next()) psi_N(iter.Full_Index())=true;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,1);iter.Valid();iter.Next()){
        TV_INT cell=iter.Cell_Index();
        surface_cells(cell)=false;
        if(psi_D(cell)) continue;
        for(int axis=1;axis<=TV::m;axis++)for(int side=-1;side<=1;side+=2){
            TV_INT cell2=cell+TV_INT::Axis_Vector(axis)*side;
            TV_INT face=(side==-1)?cell:(cell+TV_INT::Axis_Vector(axis));
            if(psi_D(cell2) && !psi_N.Component(axis)(face)){surface_cells(cell)=true;break;}}}
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iter(grid);iter.Valid();iter.Next()){
        T_FACE_INDEX face=iter.Full_Index();
        surface_faces(face)=!psi_N(face) && psi_D(iter.First_Cell_Index())!=psi_D(iter.Second_Cell_Index());}
}
//#####################################################################
// Function Build_Index_Mapping
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Build_Index_Mapping()
{
    ARRAYS_COMPUTATIONS::Fill(cell_indices.array,0);
    for(int axis=1;axis<=TV::m;axis++) ARRAYS_COMPUTATIONS::Fill(face_indices.Component(axis).array,0);
    index_cells.Resize(0);index_faces.Resize(0);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid);iter.Valid();iter.Next()){
        TV_INT cell=iter.Cell_Index();
        if(psi_D(cell)) continue;
        cell_indices(cell)=index_cells.Append(cell);}
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iter(grid);iter.Valid();iter.Next()){
        T_FACE_INDEX face=iter.Full_Index();
        TV_INT cell1=iter.First_Cell_Index(),cell2=iter.Second_Cell_Index();
        if(psi_N(face) || (psi_D(cell1) && psi_D(cell2))) continue;
        face_indices(face)=index_faces.Append(face);}
}
//#####################################################################
// Function Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Linear
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Linear()
{
    T contribution=1/number_particle_per_cell,one_over_radius=grid.one_over_dX.Max();
    HT.Reset(index_faces.m);J.Reset(index_faces.m);
    ARRAY<T> sum(index_faces.m);ARRAYS_COMPUTATIONS::Fill(sum,(T)0);
    ARRAY<T> sum2(particles.Size()*TV::m);ARRAYS_COMPUTATIONS::Fill(sum2,(T)0);
    for(int axis=1;axis<=TV::m;axis++){
        TV offset=TV::Axis_Vector(axis)*grid.dX(axis)*0.5;
        GRID<TV> face_grid=grid.Get_Face_Grid(axis);
        for(int i=1;i<=particles.Size();i++){
            const TV& p=particles.X(i);
            TV_INT face0=face_grid.Closest_Node(p);
            RANGE<TV_INT> range(face0-1,face0+1);
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(face_grid,range);iter.Valid();iter.Next()){
                int j=face_indices.Component(axis)(iter.Cell_Index());
                T weight=contribution*Linear_Kernel(p,face_grid.Node(iter.Cell_Index()),one_over_radius);
                if(weight==0) continue;
                if(j!=0){
                    sum(j)+=weight;
                    sum2(Flatten_Particle_Index(i,axis))+=weight;
                    HT.Append_Entry_To_Current_Row(j,weight);
                    J.Append_Entry_To_Current_Row(j,weight);}
                if(psi_N.Component(axis)(iter.Cell_Index()))
                    sum2(Flatten_Particle_Index(i,axis))+=weight;}
            J.Finish_Row();
            HT.Finish_Row();}}
    HT.Sort_Entries();
    J.Sort_Entries();

    for(int i=1;i<=HT.A.m;i++) if(sum(HT.A(i).j)>0) HT.A(i).a/=sum(HT.A(i).j);
    for(int i=1;i<=J.m;i++) for(int j=J.offsets(i);j<J.offsets(i+1);j++) if(sum2(i)>0) J.A(j).a/=sum2(i);
    H.Reset(0);HT.Transpose(H);
}
//#####################################################################
// Function Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Constant
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Constant()
{
    HT.Reset(index_faces.m);J.Reset(index_faces.m);
    ARRAY<int> sum(index_faces.m);ARRAYS_COMPUTATIONS::Fill(sum,0);
    ARRAY<T> sum2(particles.Size()*TV::m);ARRAYS_COMPUTATIONS::Fill(sum2,(T)0);
    for(int axis=1;axis<=TV::m;axis++){
        TV offset=TV::Axis_Vector(axis)*grid.dX(axis)*0.5;
        GRID<TV> face_grid=grid.Get_Face_Grid(axis);
        for(int i=1;i<=particles.Size();i++){
            const TV& p=particles.X(i);
            TV_INT face=face_grid.Closest_Node(p);
            int j=face_indices.Component(axis)(face);
            if(j!=0){
                sum(j)++;
                HT.Append_Entry_To_Current_Row(j,1);
                J.Append_Entry_To_Current_Row(j,1);}
            if(psi_N.Component(axis)(face)) particles.V(i)(axis)=0;
            HT.Finish_Row();
            J.Finish_Row();}}
    HT.Sort_Entries();
    J.Sort_Entries();
    for(int i=1;i<=HT.A.m;i++) if(sum(HT.A(i).j)>0) HT.A(i).a/=(T)sum(HT.A(i).j);
    for(int i=1;i<=J.m;i++) for(int j=J.offsets(i);j<J.offsets(i+1);j++) if(sum2(i)>0) J.A(j).a/=sum2(i);
    H.Reset(0);HT.Transpose(H);
}
//#####################################################################
// Function Compute_Grid_To_Surface_Mesh_Interpolation_Matrix
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Grid_To_Surface_Mesh_Interpolation_Matrix()
{
    if(surface_tension_force){
        W.Reset(index_faces.m);
        const GEOMETRY_PARTICLES<TV>& vertices=surface_tension_force->particles;
        for(int axis=1;axis<=TV::m;axis++){
            for(int i=1;i<=vertices.X.m;i++){
                const TV& p=vertices.X(i);
                TV_INT cell=grid.Cell(p,0);
                if(psi_D(cell)){
                    RANGE<TV_INT> range(cell-1,cell+1);
                    T min_distance=FLT_MAX;
                    for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,range);iter.Valid();iter.Next()){
                        TV_INT cell1=iter.Cell_Index();
                        if(grid.Inside_Domain(cell1) && !psi_D(cell1)){
                            T distance=(iter.Location()-p).Magnitude();
                            if(distance<min_distance){cell=cell1;min_distance=distance;}}}}
                TV_INT face1=cell,face2=cell+TV_INT::Axis_Vector(axis);
                int j1=face_indices.Component(axis)(face1),j2=face_indices.Component(axis)(face2);
                T alpha=(p(axis)-grid.Axis_X_Face(face1,axis)(axis))*grid.one_over_dX(axis);
                if(alpha<=0) W.Append_Entry_To_Current_Row(j1,1);
                else if(alpha>=1) W.Append_Entry_To_Current_Row(j2,1);
                else{
                    W.Append_Entry_To_Current_Row(j1,1-alpha);
                    W.Append_Entry_To_Current_Row(j2,alpha);}
                W.Finish_Row();}}
        W.Sort_Entries();
        WT.Reset(0);W.Transpose(WT);
        WH=W*H;}
}
//#####################################################################
// Function Compute_Grid_Gradient_Matrix
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Grid_Gradient_Matrix()
{
    TV face_size=grid.Face_Sizes();
    G.Reset(index_cells.m);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iter(grid);iter.Valid();iter.Next()){
        T_FACE_INDEX face=iter.Full_Index();
        int j=face_indices(face);
        if(j==0) continue;
        int i1=cell_indices(iter.First_Cell_Index()),i2=cell_indices(iter.Second_Cell_Index());
        if(i1!=0) G.Append_Entry_To_Current_Row(i1,face_size(face.axis));
        if(i2!=0) G.Append_Entry_To_Current_Row(i2,-face_size(face.axis));
        G.Finish_Row();}
    G.Sort_Entries();
}
//#####################################################################
// Function Compute_Grid_Diffusion_Matrix
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Grid_Diffusion_Matrix(const T dt)
{
    T coefficient=dt*sqrt(viscosity*grid.One_Over_Cell_Size());
    TV face_size=grid.Face_Sizes();
    Gu.Reset(index_faces.m);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iter(grid);iter.Valid();iter.Next()){
        T_FACE_INDEX face=iter.Full_Index();
        int j=face_indices(face);
        if(j==0) continue;
        for(int other_axis=1;other_axis<=TV::dimension;other_axis++){
            int i1=face_indices(T_FACE_INDEX(face.axis,face.index+TV_INT::Axis_Vector(other_axis)));
            if(i1!=0){
                Gu.Append_Entry_To_Current_Row(j,-face_size(face.axis)*coefficient);
                Gu.Append_Entry_To_Current_Row(i1,face_size(face.axis)*coefficient);}}
        Gu.Finish_Row();}
    Gu.Sort_Entries();
}
//#####################################################################
// Function Compute_Full_Coupling_Matrix
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Compute_Full_Coupling_Matrix()
{
    Mi.Resize(TV::m*particles.Size());
    if(recompute_mass) for(int i=1;i<=Mi.n;i++) if(particle_mass(i)==0) Mi(i)=0;else Mi(i)=Inverse(particle_mass(i));
    else for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) Mi(Flatten_Particle_Index(i,axis))=particles.one_over_mass(i);
    if(!use_simplified_mass){
        K.Reset(0);
        if(apply_pressure){KT=HT*G;KT.Transpose(K);}
        if(surface_tension_force && surface_tension_force->apply_implicit_forces){
            if(K.m>0) K=K.Concatenate_Rows(C*WH);
            else K=C*WH;}
        if(viscosity>0){
            if(K.m>0) K=K.Concatenate_Rows(Gu*H);
            else K=Gu*H;}
        KT.Reset(0);K.Transpose(KT);
        A=K.Times_Diagonal_Times(Mi,KT).Create_NXN_Matrix();
        //A=KT.Transpose_Times_Diagonal_Times(Mi);
        if(apply_pressure) A=Add_Diagonal_Terms_To_Matrix((T)1,G.n+1,KT.n,A);
        else A=Add_Diagonal_Terms_To_Matrix((T)1,1,KT.n,A);}
    else{
        VECTOR_ND<T> fMi(index_faces.m),delta(index_faces.m);delta.Fill((T)1);
        SPARSE_MATRIX_FLAT_MXN<T> HMiHT=H.Times_Diagonal_Times(Mi,HT);
        HMiHT.Times(delta,fMi);
        KT=G;K.Reset(0);KT.Transpose(K);
        if(surface_tension_force && surface_tension_force->apply_implicit_forces){
            K=K.Concatenate_Rows(C*W);
            KT.Reset(0);K.Transpose(KT);
            VECTOR_ND<T> sMi(C.n);delta.Resize(sMi.n);delta.Fill((T)1);
            SPARSE_MATRIX_FLAT_MXN<T> WHMiHTWT=W.Times_Diagonal_Times(fMi,WT);//W*HMiHT*WT;
            WHMiHTWT.Times(delta,sMi);
            SPARSE_MATRIX_FLAT_MXN<T> block11,block12,block21,block22,D,CT,A_temp;
            D.Reset(0);G.Transpose(D);
            CT.Reset(0);C.Transpose(CT);
            block11=D.Times_Diagonal_Times(fMi,G);
            block21=C*W.Times_Diagonal_Times(fMi,G);
            block12.Reset(0);block21.Transpose(block12);
            block22=C.Times_Diagonal_Times(sMi,CT);
            A_temp=Build_Matrix_From_Blocks(block11,block12,block21,block22);
            A_temp=Add_Diagonal_Terms_To_Matrix((T)1,G.n+1,A_temp.n,A_temp);
            A=A_temp.Create_NXN_Matrix();}
        else A=(K.Times_Diagonal_Times(fMi,KT)).Create_NXN_Matrix();}
    if(output_matrices) Output_Matrices();
}
//#####################################################################
// Function Build_Matrix_From_Block
//#####################################################################
template<class TV> SPARSE_MATRIX_FLAT_MXN<typename TV::SCALAR> GRID_PARTICLE_COUPLING<TV>::
Build_Matrix_From_Blocks(SPARSE_MATRIX_FLAT_MXN<T>& block11,SPARSE_MATRIX_FLAT_MXN<T>& block12,SPARSE_MATRIX_FLAT_MXN<T>& block21,SPARSE_MATRIX_FLAT_MXN<T>& block22)
{
    block12.n+=block11.n;for(int i=1;i<=block12.A.m;i++) block12.A(i).j+=block11.n;
    block22.n+=block21.n;for(int i=1;i<=block22.A.m;i++) block22.A(i).j+=block21.n;
    block11.n=block12.n;block21.n=block22.n;
    return (block11+block12).Concatenate_Rows(block21+block22);
}
//#####################################################################
// Function Add_Diagonal_Terms_To_Matrix
//#####################################################################
template<class TV> SPARSE_MATRIX_FLAT_NXN<typename TV::SCALAR> GRID_PARTICLE_COUPLING<TV>::
Add_Diagonal_Terms_To_Matrix(const T diagonal_term,const int start_row,const int end_row,const SPARSE_MATRIX_FLAT_NXN<T>& matrix)
{
    if(!(matrix.n<=end_row && start_row<=end_row && 1<=start_row)) return matrix;
    ARRAY<int> row_lengths(matrix.n);
    for(int i=1;i<=matrix.n;i++){
        int begin=matrix.offsets(i),end=matrix.offsets(i+1);
        row_lengths(i)=end-begin;
        if(i>=start_row && i<=end_row){
            bool has_diagonal=false;
            for(int j=begin;j<end;j++) if(matrix.A(j).j==i){
                    has_diagonal=true;break;}
            if(!has_diagonal) row_lengths(i)++;}}
    SPARSE_MATRIX_FLAT_NXN<T> result;
    result.Set_Row_Lengths(row_lengths);result.n=matrix.n;

    int index=matrix.offsets(1);
    for(int i=1;i<=matrix.n;i++){int end=matrix.offsets(i+1);
        for(;index<end;index++) result(i,matrix.A(index).j)=matrix.A(index).a;
        if(i>=start_row && i<=end_row) result(i,i)+=diagonal_term;}

    return result;
}
//#####################################################################
// Function Add_Diagonal_Terms_To_Matrix
//#####################################################################
template<class TV> SPARSE_MATRIX_FLAT_MXN<typename TV::SCALAR> GRID_PARTICLE_COUPLING<TV>::
Add_Diagonal_Terms_To_Matrix(const T diagonal_term,const int start_row,const int end_row,const SPARSE_MATRIX_FLAT_MXN<T>& matrix)
{
    if(!(matrix.n<=end_row && start_row<=end_row && 1<=start_row)) return matrix;
    ARRAY<int> row_lengths(matrix.m);
    for(int i=1;i<=matrix.m;i++){
        int begin=matrix.offsets(i),end=matrix.offsets(i+1);
        row_lengths(i)=end-begin;
        if(i>=start_row && i<=end_row){
            bool has_diagonal=false;
            for(int j=begin;j<end;j++) if(matrix.A(j).j==i){
                    has_diagonal=true;break;}
            if(!has_diagonal) row_lengths(i)++;}}
    SPARSE_MATRIX_FLAT_MXN<T> result;
    result.Set_Row_Lengths(row_lengths);result.n=matrix.n;

    int index=matrix.offsets(1);
    for(int i=1;i<=matrix.m;i++){int end=matrix.offsets(i+1);
        for(;index<end;index++) result(i,matrix.A(index).j)=matrix.A(index).a;
        if(i>=start_row && i<=end_row) result(i,i)+=diagonal_term;}

    return result;
}
//#####################################################################
// Function Solve_Coupling_System
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Solve_Coupling_System(const T dt)
{
    x.Resize(KT.n);x.Fill(0);
    b.Resize(KT.n);b.Fill(0);
    q.Resize(KT.n);q.Fill(0);
    s.Resize(KT.n);s.Fill(0);
    r.Resize(KT.n);r.Fill(0);
    k.Resize(KT.n);k.Fill(0);
    z.Resize(KT.n);z.Fill(0);
    v.Resize(TV::m*particles.Size());
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) v(Flatten_Particle_Index(i,axis))=particles.V(i)(axis);
    if(!use_simplified_mass) K.Times(v,b);
    else{
        VECTOR_ND<T> fv(index_faces.m);
        H.Times(v,fv);
        K.Times(fv,b);}
    if(fix_volume_error) Add_Volume_Error_To_RHS(b,dt);
    PCG_SPARSE<T> solver;
    if(!use_IC_preconditioner) solver.Use_Conjugate_Gradient();
    solver.Set_Maximum_Iterations(max_cg_iterations);
    solver.Show_Results();
    //solver.Show_Residuals();
    solver.Solve(A,x,b,q,s,r,k,z,(T)1e-8);
    if(!use_simplified_mass) KT.Times(x,v);
    else{
        VECTOR_ND<T> fv(index_faces.m);
        KT.Times(x,fv);
        HT.Times(fv,v);}
    v*=Mi;
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) particles.V(i)(axis)-=v(Flatten_Particle_Index(i,axis));

    if(1 || output_face_velocities){
        VECTOR_ND<T> face_velocity_vector(index_faces.m);
        for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) v(Flatten_Particle_Index(i,axis))=particles.V(i)(axis);
        H.Times(v,face_velocity_vector);
        for(int axis=1;axis<=TV::m;axis++) ARRAYS_COMPUTATIONS::Fill(face_velocities.Component(axis).array,(T)0);
        for(int i=1;i<=index_faces.m;i++) face_velocities(index_faces(i))=face_velocity_vector(i);
        ARRAYS_COMPUTATIONS::Fill(pressure.array,(T)0);
        for(int i=1;i<=index_cells.m;i++) pressure(index_cells(i))=x(i);}
}
//#####################################################################
// Function Project_Out_Null_Space_Of_Nullspace
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Project_Out_Null_Space_Of_Nullspace()
{
    x.Resize(H.m);x.Fill(0);
    b.Resize(H.m);b.Fill(0);
    q.Resize(H.m);q.Fill(0);
    s.Resize(H.m);s.Fill(0);
    r.Resize(H.m);r.Fill(0);
    k.Resize(H.m);k.Fill(0);
    z.Resize(H.m);z.Fill(0);
    v.Resize(H.n);
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) v(Flatten_Particle_Index(i,axis))=particles.V(i)(axis);
    H.Times(v,b);
    SPARSE_MATRIX_FLAT_NXN<T> HHT=(H*HT).Create_NXN_Matrix();
    PCG_SPARSE<T> solver;
    solver.Set_Maximum_Iterations(max_cg_iterations);
    solver.Show_Results();
    //solver.Show_Residuals();
    solver.Solve(HHT,x,b,q,s,r,k,z,(T)1e-8);
    HT.Times(x,v);
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) particles.V(i)(axis)=v(Flatten_Particle_Index(i,axis));
}
//#####################################################################
// Function Add_Volume_Error_To_RHS
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Add_Volume_Error_To_RHS(VECTOR_ND<T>& b,const T dt)
{
    if(use_log_based_volume_correction){
        T dtau=1,normalization=grid.Cell_Size()/dtau;
        for(int i=1;i<=index_cells.m;i++){TV_INT cell=index_cells(i);
            T weight=cell_weights(cell);
            if(!surface_cells(cell) || weight>1) b(i)-=std::log(weight)*normalization;}}
    else{
        T normalization=grid.Cell_Size()/dt;
        for(int i=1;i<=index_cells.m;i++){TV_INT cell=index_cells(i);
            T weight=cell_weights(cell);
            if(!surface_cells(cell) || weight>1) b(i)-=(weight-1)*normalization;}}
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Add_External_Forces(const T dt)
{
    ARRAYS_COMPUTATIONS::Fill(particles.F,TV());
    VECTOR_ND<T> face_force_vector(index_faces.m),particle_force_vector(particles.X.m*TV::m);
    for(int i=1;i<=index_faces.m;i++) face_force_vector(i)=face_forces(index_faces(i))*face_mass(i);
    HT.Times(face_force_vector,particle_force_vector);
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.F.m;i++){
        int j=Flatten_Particle_Index(i,axis);
        particles.F(i)(axis)=particle_force_vector(j)/particle_mass(j);}
    particles.Euler_Step_Velocity(dt);
}
//#####################################################################
// Function Add_Explicit_Forces
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Add_Explicit_Forces(const T time,const T dt)
{
    if(!surface_tension_force) return;
    ARRAYS_COMPUTATIONS::Fill(particles.F,TV());
    ARRAY<TV> vertex_force(surface_tension_force->particles.X.m);
    ARRAYS_COMPUTATIONS::Fill(vertex_force,TV());
    surface_tension_force->Add_Velocity_Independent_Forces(vertex_force,time);
    int index=WH.offsets(1);
    for(int i=1;i<=WH.m;i++){
        int axis=(i-1)/vertex_force.m+1;int j=(i-1)%vertex_force.m+1;
        int end=WH.offsets(i+1);T y=vertex_force(j)(axis);
        for(;index<end;index++){
            int axis2=(WH.A(index).j-1)/particles.F.m+1,j2=(WH.A(index).j-1)%particles.F.m+1;
            particles.F(j2)(axis2)+=WH.A(index).a*y;}}
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.F.m;i++){
        int j=Flatten_Particle_Index(i,axis);
        particles.F(i)(axis)/=particle_mass(j);}
    particles.Euler_Step_Velocity(dt);
}
//#####################################################################
// Function Interpolate_From_Grid_To_Particles
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Interpolate_From_Grid_To_Particles(const T dt)
{
    VECTOR_ND<T> fv(index_faces.m),pic_v(TV::m*particles.Size());
    dv.Resize(TV::m*particles.Size());
    v.Resize(TV::m*particles.Size());
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) v(Flatten_Particle_Index(i,axis))=particles.V(i)(axis);
    H.Times(v,fv);
    J.Times(fv,pic_v);
    for(int i=1;i<=J.offsets.m;i++) if(J.offsets(i)==J.offsets(i+1)) pic_v(i)=v(i);
    for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++){
        int j=Flatten_Particle_Index(i,axis);
        T delta=pic_ratio*(pic_v(j)-particles.V(i)(axis));
        particles.V(i)(axis)+=delta;
        dv(j)=-delta;}
    if(feed_back_for_momentum_conservation){
        H.Times(dv,fv);
        if(not_feed_back_to_surface_faces)
            for(int i=1;i<=index_faces.m;i++){
                const T_FACE_INDEX& face=index_faces(i);
                if(surface_faces(face)) fv(i)=0;}
        HT.Times(fv,dv);
        for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=particles.Size();i++) particles.V(i)(axis)+=dv(Flatten_Particle_Index(i,axis));}
}
//#####################################################################
// Function Output_Matrices
//#####################################################################
template<class TV> void GRID_PARTICLE_COUPLING<TV>::
Output_Matrices()const
{
    OCTAVE_OUTPUT<T>("A.txt").Write("A",A);
    OCTAVE_OUTPUT<T>("G.txt").Write("G",G);
    OCTAVE_OUTPUT<T>("HT.txt").Write("HT",HT);
    OCTAVE_OUTPUT<T>("C.txt").Write("C",C);
    OCTAVE_OUTPUT<T>("W.txt").Write("W",W);
    OCTAVE_OUTPUT<T>("WH.txt").Write("WH",WH);
}
//#####################################################################
template class GRID_PARTICLE_COUPLING<VECTOR<float,2> >;
template class GRID_PARTICLE_COUPLING<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_PARTICLE_COUPLING<VECTOR<double,2> >;
template class GRID_PARTICLE_COUPLING<VECTOR<double,3> >;
#endif
