//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Parallel_Computation/REFINEMENT_THREADS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
PROJECTION_REFINEMENT_UNIFORM(const T_GRID& mac_grid,const int scale,const T alpha_in,const bool flame_input,const bool multiphase,const bool use_variable_beta,const bool use_poisson,THREAD_QUEUE* thread_queue_input)
    :PROJECTION_DYNAMICS_UNIFORM<T_GRID>(GRID<TV>(TV_INT(),RANGE<TV>(),true),flame_input,multiphase,use_variable_beta,use_poisson,thread_queue_input),thread_queue(thread_queue_input),coarse_mpi_grid(0),fine_mpi_grid(0),
    fast_local_projection(scale),beta_face(use_poisson?poisson->beta_face:*new T_FACE_ARRAYS_SCALAR()),alpha(alpha_in),coarse_scale(scale)
{
    for(int i=1;i<=TV::dimension;i++) for(int j=1;j<=2;j++){domain_boundary(i)(j)=true;solid_wall(i)(j)=false;}
    scale_face_inverse=(TV::dimension==2)?((T)1./scale):((T)1./(scale*scale));
    if(poisson) poisson->Set_Variable_beta(true);
    this->Initialize_Grid(mac_grid); //TODO: Fix this hack
    assert(!poisson || &beta_face==&poisson->beta_face);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
PROJECTION_REFINEMENT_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input,const int scale,const T alpha_in)
    :PROJECTION_DYNAMICS_UNIFORM<T_GRID>(mac_grid,levelset_input),thread_queue(0),coarse_mpi_grid(0),fine_mpi_grid(0),fast_local_projection(scale),
     beta_face(poisson->beta_face),alpha(alpha_in),coarse_scale(scale)
{
    for(int i=1;i<=TV::dimension;i++) for(int j=1;j<=2;j++){domain_boundary(i)(j)=true;solid_wall(i)(j)=false;}
    scale_face_inverse=(TV::dimension==2)?((T)1./scale):((T)1./(scale*scale));
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
~PROJECTION_REFINEMENT_UNIFORM()
{
    if(!poisson) delete &beta_face;
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Initialize_Grid(const T_GRID& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
    fine_grid=mac_grid;
    if(elliptic_solver->mpi_grid){
        fine_mpi_grid=elliptic_solver->mpi_grid;
        coarse_grid.Initialize(fine_mpi_grid->global_grid.counts/coarse_scale,fine_mpi_grid->global_grid.domain,fine_grid.Is_MAC_Grid());
        if(elliptic_solver->mpi_grid->threaded_grid) coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(elliptic_solver->mpi_grid->threaded_grid->buffers,elliptic_solver->mpi_grid->threaded_grid->tid,elliptic_solver->mpi_grid->threaded_grid->number_of_processes,coarse_grid,3);            
        else coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(coarse_grid,3);
        coarse_mpi_grid->Initialize(domain_boundary);}
    else coarse_grid.Initialize(fine_grid.counts/coarse_scale,fine_grid.domain,fine_grid.Is_MAC_Grid());
    local_grid.Initialize(TV_INT::All_Ones_Vector()*coarse_scale,RANGE<TV>::Centered_Box(),fine_grid.Is_MAC_Grid());
    fast_local_projection.Initialize_Grid(local_grid);
    fast_local_projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
    fast_local_projection.elliptic_solver->pcg.Use_Modified_Incomplete_Cholesky();
    local_face_velocities.Resize(local_grid);
    coarse_face_velocities.Resize(coarse_grid,3);
    coarse_face_velocities_save.Resize(coarse_grid,3);
    face_velocities_save.Resize(fine_grid,3);
    if(!poisson) beta_face.Resize(coarse_grid);
    fine_psi_N.Resize(fine_grid);
}
//#####################################################################
// Function Average_Velocities_From_Fine_To_Coarse
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Average_Velocities_From_Fine_To_Coarse(ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& fine_face_velocities)
{
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()) coarse_face_velocities(iterator.Full_Index())=0;
    FACE_INDEX<TV::dimension> fine_index;
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){int axis=local_iterator.Axis();
            fine_index=FACE_INDEX<TV::dimension>(axis,(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*coarse_scale+local_iterator.Face_Index());
            if(local_iterator.First_Boundary()) coarse_face_velocities(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))+=fine_face_velocities(fine_index);
            else coarse_face_velocities(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))+=fine_face_velocities(fine_index);}}    
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        T factor=0;
        if(iterator.First_Cell_Index()(axis)>0){factor++;}
        if(iterator.Second_Cell_Index()(axis)<=coarse_grid.Counts()(axis)){factor++;}
        coarse_face_velocities.Component(axis)(iterator.Face_Index())*=factor==1?scale_face_inverse:(T)0.5*scale_face_inverse;}
}
//#####################################################################
// Function Set_Beta_Face_For_Boundary_Conditions
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Set_Beta_Face_For_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities)
{ 
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()) beta_face.Component(iterator.Axis())(iterator.Face_Index())=0;
    FACE_INDEX<TV::dimension> fine_index;
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){int axis=local_iterator.Axis();
            fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*coarse_scale+local_iterator.Face_Index());
            if(!fine_psi_N(fine_index)){
                if(local_iterator.First_Boundary()) beta_face(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))+=T(1);
                else beta_face(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))+=T(1);}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){
        if(beta_face.Component(iterator.Axis())(iterator.Face_Index())==0) elliptic_solver->psi_N(iterator.Full_Index())=true;
        else{
            int factor=0;
            if(iterator.First_Cell_Index()(iterator.Axis())>0){factor++;}
            if(iterator.Second_Cell_Index()(iterator.Axis())<=coarse_grid.Counts()(iterator.Axis())){factor++;}
            beta_face.Component(iterator.Axis())(iterator.Face_Index())*=(factor==1)?scale_face_inverse:(T)0.5*scale_face_inverse;}}
}
//#####################################################################
// Function Set_Coarse_Boundary_Conditions
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Set_Coarse_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities)
{
    for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
        if(domain_boundary(axis)(axis_side)){ //Need to check mpi as smaller solves are never using mpi (for now)
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(coarse_grid,1,GRID<TV>::BOUNDARY_REGION,side);local_iterator.Valid();local_iterator.Next()){
                TV_INT cell=local_iterator.Face_Index()+interior_cell_offset;
                TV_INT boundary_face=axis_side==1?local_iterator.Face_Index()+TV_INT::Axis_Vector(axis):local_iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(solid_wall(axis)(axis_side) && coarse_face_velocities.Component(axis).Valid_Index(boundary_face)){
                    elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,boundary_face))=true;coarse_face_velocities(FACE_INDEX<TV::dimension>(axis,boundary_face))=0;}
                else{
                    elliptic_solver->psi_D(cell)=true;p(cell)=0;}}}}
    Set_Beta_Face_For_Boundary_Conditions(coarse_face_velocities);
}
//#####################################################################
// Function Map_Fine_To_Coarse
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Map_Fine_To_Coarse(T_FACE_ARRAYS_SCALAR& coarse_face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities)
{
    fine_psi_N=elliptic_solver->psi_N;
    BASE::Initialize_Grid(coarse_grid);
    elliptic_solver->mpi_grid=coarse_mpi_grid;
    elliptic_solver->psi_D.Fill(false);elliptic_solver->psi_N.Fill(false);
    Set_Coarse_Boundary_Conditions(coarse_face_velocities);
    Average_Velocities_From_Fine_To_Coarse(coarse_face_velocities,face_velocities);
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()) coarse_face_velocities_save(iterator.Full_Index())=coarse_face_velocities(iterator.Full_Index());
    for(typename GRID<TV>::FACE_ITERATOR iterator(fine_grid);iterator.Valid();iterator.Next()) face_velocities_save(iterator.Full_Index())=face_velocities(iterator.Full_Index());
}
//#####################################################################
// Function Map_Coarse_To_Fine
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities)
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));sum.Fill(0);
        FACE_INDEX<TV::dimension> fine_index,sum_index=FACE_INDEX<TV::dimension>(1,TV_INT::All_Ones_Vector());
        if(alpha!=1){
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){int axis=local_iterator.Axis();
                TV_INT offset;offset(axis)++;
                sum_index=FACE_INDEX<TV::dimension>(axis,(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
                sum(sum_index)=0;}
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){int axis=local_iterator.Axis();
                TV_INT offset;offset(axis)++;
                fine_index=FACE_INDEX<TV::dimension>(axis,cell*coarse_scale+local_iterator.Face_Index()-coarse_scale*TV_INT::All_Ones_Vector());
                sum_index=FACE_INDEX<TV::dimension>(axis,(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
                if(fine_psi_N(fine_index)) sum(sum_index)+=face_velocities_save(fine_index);}}
        //update based on (1-a)*((1/A)*(Vco-avg Vs))+a*(Vfo+(Vcn-Vco)*(1/A))
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){int axis=local_iterator.Axis();
            FACE_INDEX<TV::dimension> coarse_index,fine_index=FACE_INDEX<TV::dimension>(axis,cell*coarse_scale+local_iterator.Face_Index()-coarse_scale*TV_INT::All_Ones_Vector());
            if(local_iterator.First_Boundary()) coarse_index=FACE_INDEX<TV::dimension>(axis,coarse_grid.First_Face_Index_In_Cell(axis,cell));
            else coarse_index=FACE_INDEX<TV::dimension>(axis,coarse_grid.Second_Face_Index_In_Cell(axis,cell));
            T area=beta_face(coarse_index),one_over_area=fine_psi_N(fine_index)?0:(T)1./area;
            T beta=fine_psi_N(fine_index)?0:(T)1;
            T total_area=TV::dimension==2?(T)coarse_scale:(T)coarse_scale*(T)coarse_scale;
            if(alpha!=1){TV_INT offset;offset(axis)++;sum_index=FACE_INDEX<TV::dimension>(axis,local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
            if(!fine_psi_N(fine_index))
                fine_face_velocities(fine_index)=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocities(coarse_index)-sum(sum_index)/total_area)):0)+
                    (alpha?(alpha*(face_velocities_save(fine_index)+one_over_area*beta*(coarse_face_velocities(coarse_index)-coarse_face_velocities_save(coarse_index)))):0);}}
}
//#####################################################################
// Function Map_Fine_To_Local_Boundary_For_Cell
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Map_Fine_To_Local_Boundary_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT coarse_cell_index)
{
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*coarse_scale+local_iterator.Face_Index()-coarse_scale*TV_INT::All_Ones_Vector());
        local_face_velocities(local_iterator.Full_Index())=fine_face_velocities(fine_index);}
}
//#####################################################################
// Function Map_Fine_To_Local_Interior_For_Cell
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Map_Fine_To_Local_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT coarse_cell_index,bool zero_out)
{
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*coarse_scale+local_iterator.Face_Index()-coarse_scale*TV_INT::All_Ones_Vector());
        if(zero_out) local_face_velocities(local_iterator.Full_Index())=0;
        else local_face_velocities(local_iterator.Full_Index())=face_velocities_save(fine_index);}
}
//#####################################################################
// Function Map_Fine_To_Local_Boundaries_For_Cell
//#####################################################################
template<class T_GRID> bool PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Map_Fine_To_Local_Boundaries_For_Cell(GRID<TV>& local_mac_grid,ARRAY<bool,FACE_INDEX<TV::dimension> >& local_psi_N,TV_INT cell_index)
{
    bool has_solids=false;
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        local_psi_N(local_iterator.Full_Index())=fine_psi_N(FACE_INDEX<TV::dimension>(local_iterator.Axis(),cell_index*coarse_scale+local_iterator.Face_Index()-coarse_scale*TV_INT::All_Ones_Vector()));
        if(local_psi_N(local_iterator.Full_Index())) has_solids=true;}
    return has_solids;
}
//#####################################################################
// Function Map_Fine_To_Local_Interior_For_Cell
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Map_Local_To_Fine_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT cell_index)
{
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
        fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),cell_index*coarse_scale+local_iterator.Face_Index()-coarse_scale*TV_INT::All_Ones_Vector()))=local_face_velocities(local_iterator.Full_Index());}
}
//#####################################################################
// Function Local_Projection_PCG
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,T_GRID& local_grid,T_FACE_ARRAYS_SCALAR& local_face_velocities,FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& local_projection,const T dt,const T time,TV_INT cell_index)
{ 
    Map_Fine_To_Local_Boundary_For_Cell(local_grid,local_face_velocities,fine_face_velocities,cell_index);
    Map_Fine_To_Local_Interior_For_Cell(local_grid,local_face_velocities,fine_face_velocities,cell_index,false);
    bool contains_solids=Map_Fine_To_Local_Boundaries_For_Cell(local_grid,local_projection.elliptic_solver->psi_N,cell_index);
    for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_grid);local_iterator.Valid();local_iterator.Next()){
        local_projection.p(local_iterator.Cell_Index())=p(cell_index);}
    local_projection.elliptic_solver->Set_Neumann_Outer_Boundaries();
    local_projection.p*=dt;        
    if(contains_solids) local_projection.Make_Divergence_Free(local_face_velocities,dt,time);
    else local_projection.Make_Divergence_Free_Fast(local_face_velocities,dt,time);
    local_projection.p/=dt;
    Map_Local_To_Fine_Interior_For_Cell(local_grid,local_face_velocities,fine_face_velocities,cell_index);
}
//#####################################################################
// Function Threaded_Local_Projection_PCG
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Threaded_Local_Projection_PCG(RANGE<TV_INT>& domain,T_FACE_ARRAYS_SCALAR& fine_face_velocities,const T dt,const T time)
{
    ARRAY<T,FACE_INDEX<TV::dimension> > domain_face_velocities(local_grid);
    FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > domain_local_projection(fast_local_projection);
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid,domain);iterator.Valid();iterator.Next()) Local_Projection_PCG(fine_face_velocities,local_grid,domain_face_velocities,domain_local_projection,dt,time,iterator.Cell_Index());
}
//#####################################################################
// Function Local_Projection_PCG
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,const T dt,const T time)
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()) Local_Projection_PCG(fine_face_velocities,local_grid,local_face_velocities,fast_local_projection,dt,time,iterator.Cell_Index());
}
//#####################################################################
// Function Map_Coarse_To_Fine
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    Map_Coarse_To_Fine(coarse_face_velocities,face_velocities);
    if(thread_queue) DOMAIN_ITERATOR_THREADED_ALPHA<PROJECTION_REFINEMENT_UNIFORM<T_GRID>,TV>(coarse_grid.Domain_Indices(),thread_queue,1,0,1,8).template Run<T_FACE_ARRAYS_SCALAR&,T,T>(*this,&PROJECTION_REFINEMENT_UNIFORM<T_GRID>::Threaded_Local_Projection_PCG,face_velocities,dt,time);        
    else Local_Projection_PCG(face_velocities,dt,time);
    elliptic_solver->mpi_grid=fine_mpi_grid;
    BASE::Initialize_Grid(fine_grid);
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class T_GRID> void PROJECTION_REFINEMENT_UNIFORM<T_GRID>::
Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    //LOG::Time("Map Fine To Coarse");
    Map_Fine_To_Coarse(coarse_face_velocities,face_velocities);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after mapping",0,1);
    //LOG::Time("Solve");
    BASE::Make_Divergence_Free(coarse_face_velocities,dt,time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after solve",0,1);
    //LOG::Time("Map Coarse To Fine");
    Map_Coarse_To_Fine(coarse_face_velocities,face_velocities,dt,time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after complete",0,1);
}
//#####################################################################
template class PROJECTION_REFINEMENT_UNIFORM<GRID<VECTOR<float,1> > >;
template class PROJECTION_REFINEMENT_UNIFORM<GRID<VECTOR<float,2> > >;
template class PROJECTION_REFINEMENT_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_REFINEMENT_UNIFORM<GRID<VECTOR<double,1> > >;
template class PROJECTION_REFINEMENT_UNIFORM<GRID<VECTOR<double,2> > >;
template class PROJECTION_REFINEMENT_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
