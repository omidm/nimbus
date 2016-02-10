//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_PARTICLE_STATE.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
namespace PhysBAM{
//#####################################################################
static const int self_collision_subcycles=4;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
INCOMPRESSIBLE_FINITE_VOLUME(STRAIN_MEASURE<TV,d>& strain_measure)
    :DEFORMABLES_FORCES<TV>(dynamic_cast<PARTICLES<TV>&>(strain_measure.particles)),strain_measure(strain_measure),
    disable_projection(false),minimum_volume_recovery_time_scale(0),max_cg_iterations(20),mpi_solids(0),merge_at_boundary(false),use_neumann(false),
    use_self_moving_projection(true),use_rigid_clamp_projection(true),use_diagonal_preconditioner(false),repulsions(0)
{
    T_MESH& mesh=strain_measure.mesh;
    if(!strain_measure.mesh.boundary_mesh) strain_measure.mesh.Initialize_Boundary_Mesh();
    if(!mesh.node_on_boundary) mesh.Initialize_Node_On_Boundary();
    if(!mesh.boundary_nodes) mesh.Initialize_Boundary_Nodes();

    std::stringstream ss;ss<<"self_collision_subcycles "<<self_collision_subcycles<<std::endl;LOG::filecout(ss.str());
    T_BOUNDARY_MESH& boundary_mesh=*strain_measure.mesh.boundary_mesh;
    if(!boundary_mesh.incident_elements) boundary_mesh.Initialize_Incident_Elements();
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
    boundary_to_element.Resize(boundary_mesh.elements.m);
    for(int t=1;t<=mesh.elements.m;t++){VECTOR<int,d+1>& element=mesh.elements(t);
        if(VECTOR<bool,d+1>(mesh.node_on_boundary->Subset(element)).Number_True()>=element.m-1) // using Number_True directly on the subset hits a compiler bug in gcc 4.1.1
            for(int i=1;i<=element.m;i++){
                int b=boundary_mesh.Simplex(element.Remove_Index(i));
                if(b) boundary_to_element(b).Set(t,i);}}

    node_regions.Resize(particles.array_collection->Size());
    for(int p=1;p<=particles.array_collection->Size();p++){ARRAY<int>& incident=(*mesh.incident_elements)(p);
        for(int j=1;j<=incident.m;j++) node_regions(p).Append_Unique_Elements(strain_measure.mesh.elements(incident(j)));}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
~INCOMPRESSIBLE_FINITE_VOLUME()
{}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> INCOMPRESSIBLE_FINITE_VOLUME<TV,d>* INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Create(T_OBJECT& object,const bool verbose)
{
    STRAIN_MEASURE<TV,d>* strain_measure=new STRAIN_MEASURE<TV,d>(object);
    if(verbose) strain_measure->Print_Altitude_Statistics();
    return new INCOMPRESSIBLE_FINITE_VOLUME(*strain_measure);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    // TODO: Do not compute force_particles_of_fragment here and in base class.
    force_elements.Update(strain_measure.mesh.elements,particle_is_simulated);
    force_dynamic_particles.Update(strain_measure.mesh.elements.Flattened(),particle_is_simulated);
    force_boundary_elements.Update(strain_measure.mesh.boundary_mesh->elements,particle_is_simulated);
    force_dynamic_particles_list.Remove_All();
    for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()) force_dynamic_particles_list.Append(iterator.Data());
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    BASE::Update_Position_Based_State(time,is_position_update);
    if(MPI_WORLD::Initialized() && !mpi_solids) PHYSBAM_FATAL_ERROR();
    ARRAY<T> element_volumes(strain_measure.mesh.elements.m,false); // TODO: this is not efficient.

    Bs_per_node.Resize(strain_measure.mesh.elements.m,false,false);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        T_MATRIX Ds=strain_measure.Ds(particles.X,t);
        element_volumes(t)=(T)1/FACTORIAL<d>::value*Ds.Parallelepiped_Measure();
        Bs_per_node(t)=(T)1/FACTORIAL<d+1>::value*Ds.Cofactor_Matrix();}

    boundary_normals.Resize(strain_measure.mesh.boundary_mesh->elements.m,false,false);

    for(ELEMENT_ITERATOR iterator(force_boundary_elements);iterator.Valid();iterator.Next()){int b=iterator.Data();
        int interior,i;boundary_to_element(b).Get(interior,i);
        if(i>1) boundary_normals(b)=-Bs_per_node(interior).Column(i-1);
        else boundary_normals(b)=Bs_per_node(interior)*VECTOR<T,d>::All_Ones_Vector();}

    volumes_full.Resize(particles.array_collection->Size(),false,false);
    INDIRECT_ARRAY<ARRAY<T>,ARRAY_VIEW<int>&> volumes_subset=volumes_full.Subset(strain_measure.mesh.elements.Flattened());
    ARRAYS_COMPUTATIONS::Fill(volumes_subset,(T)0);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        volumes_full.Subset(strain_measure.mesh.elements(t))+=element_volumes(t);}

    total_volume=0;
    for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        volumes_full(p)*=(T)1/(d+1);
        total_volume+=volumes_full(p);}

    if(mpi_solids) total_volume=mpi_solids->Reduce_Add_Global(total_volume);

    if(!rest_volumes_full.m){
        total_rest_volume=total_volume;
        rest_volumes_full.Resize(particles.array_collection->Size());
        for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            rest_volumes_full(p)=volumes_full(p);}}

        // TODO(jontg): Make it work bleh.
    //if(!mpi_solids) LOG::cout<<"boundary volume: "<<ARRAYS_COMPUTATIONS::Sum(rest_volumes_full.Subset(*strain_measure.mesh.boundary_nodes))<<std::endl;

    Update_Preconditioner();
}
//#####################################################################
// Class POISSON_SYSTEM
//#####################################################################
namespace{
template<class TV,int d>
class POISSON_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
protected:
    const INCOMPRESSIBLE_FINITE_VOLUME<TV,d>& fvm;
    const ARRAY<int>& dynamic_particles;
public:
    typedef INDIRECT_ARRAY<ARRAY<T> > VECTOR_T;
    typedef KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<T> > > KRYLOV_VECTOR_T;

    POISSON_SYSTEM(const INCOMPRESSIBLE_FINITE_VOLUME<TV,d>& fvm)
        :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),
        fvm(fvm),dynamic_particles(fvm.force_dynamic_particles_list)
    {}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& bp,KRYLOV_VECTOR_BASE<T>& bresult) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& p=debug_cast<const KRYLOV_VECTOR_T&>(bp);KRYLOV_VECTOR_T& result=debug_cast<KRYLOV_VECTOR_T&>(bresult);
    const KRYLOV_VECTOR_T* use_p=&p;
    INDIRECT_ARRAY<const ARRAY<T> > diagonal_preconditioner(fvm.diagonal_preconditioner_full,dynamic_particles);
    if(fvm.use_diagonal_preconditioner){use_p=&result;result.v=p.v;result.v*=diagonal_preconditioner;}
    fvm.Gradient(use_p->v.array,fvm.gradient_full);
    fvm.gradient_full.Subset(dynamic_particles)*=fvm.particles.one_over_mass.Subset(dynamic_particles);
    fvm.Project_Vector_Field(fvm.gradient_full);
    fvm.Negative_Divergence(fvm.gradient_full,result.v.array);
    if(fvm.use_diagonal_preconditioner) result.v*=diagonal_preconditioner;}

    void Project(KRYLOV_VECTOR_BASE<T>& p) const PHYSBAM_OVERRIDE
    {}

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& p) const PHYSBAM_OVERRIDE
    {Project(p);}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx),&y=debug_cast<const KRYLOV_VECTOR_T&>(by);
    assert(x.v.Size()==dynamic_particles.Size());
    T inner_product=ARRAYS_COMPUTATIONS::Dot_Product(x.v,y.v);
    if(fvm.mpi_solids) inner_product=fvm.mpi_solids->Reduce_Add(inner_product);
    return inner_product;}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);
    assert(x.v.Size()==dynamic_particles.Size());
    T convergence_norm=ARRAYS_COMPUTATIONS::Maximum_Magnitude(x.v);
    if(fvm.mpi_solids) convergence_norm=fvm.mpi_solids->Reduce_Max(convergence_norm);
    return convergence_norm;}

    T Magnitude(const KRYLOV_VECTOR_BASE<T>& x) const
    {return sqrt((T)Inner_Product(x,x));}
};
}
//#####################################################################
// Function Make_Incompressible
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Make_Incompressible(const T dt,const bool correct_volume)
{
    LOG::SCOPE scope("PROJECTION","projection, dt = %g",dt);
    if(correct_volume) {std::stringstream ss;ss<<"Correcting Volume"<<std::endl;LOG::filecout(ss.str());}
    else {std::stringstream ss;ss<<"Correcting Divergence"<<std::endl;LOG::filecout(ss.str());}

    POISSON_SYSTEM<TV,d> system(*this);
    T max_error=0;
    
    for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        if(rest_volumes_full(p)){T error=(volumes_full(p)-rest_volumes_full(p))/rest_volumes_full(p);max_error=max(max_error,abs(error));}}
    if(mpi_solids) max_error=mpi_solids->Reduce_Max(max_error);
    std::stringstream ss0;ss0<<"max error = "<<max_error<<", total volume = "<<total_volume<<", error = "<<Robust_Divide(total_volume-total_rest_volume,total_rest_volume)<<std::endl;LOG::filecout(ss0.str());

    if(disable_projection) return;

    if(boundary_pressures.m){
        if(TV::m!=d) PHYSBAM_FATAL_ERROR();
        gradient_full.Resize(particles.array_collection->Size(),false,false);
        INDIRECT_ARRAY<ARRAY<TV>,ARRAY<int>&> gradient_subset=gradient_full.Subset(force_dynamic_particles_list);
        ARRAYS_COMPUTATIONS::Fill(gradient_subset,TV());
        T_BOUNDARY_MESH& boundary_mesh=*strain_measure.mesh.boundary_mesh;
        for(ELEMENT_ITERATOR iterator(force_boundary_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            VECTOR<int,d>& element=boundary_mesh.elements(t);
            INDIRECT_ARRAY<ARRAY<T>,VECTOR<int,d>&> boundary_pressures_subset=boundary_pressures.Subset(element);
            T p_sum=ARRAYS_COMPUTATIONS::Sum(boundary_pressures_subset);
            for(int i=1;i<=d;i++) gradient_full(element[i])-=(p_sum+boundary_pressures(element[i]))*boundary_normals(t);}
        for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            particles.V(p)+=particles.one_over_mass(p)*gradient_full(p);}}

    Negative_Divergence(particles.V,divergence_full);
    KRYLOV_VECTOR_T divergence(divergence_full,force_dynamic_particles_list);

    if(correct_volume){
        T one_over_dt=1/dt,maximum_volume_recovery_fraction=dt/max(minimum_volume_recovery_time_scale,(T)1e-10);
        for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            T volume_error=rest_volumes_full(p)-volumes_full(p);
            volume_error=sign(volume_error)*min(abs(volume_error),maximum_volume_recovery_fraction*rest_volumes_full(p));
            divergence_full(p)+=one_over_dt*volume_error;}}

    system.Project(divergence);
    std::stringstream ss1;ss1<<"divergence magnitude = "<<system.Magnitude(divergence)<<std::endl;LOG::filecout(ss1.str());

    pressure_full.Resize(particles.array_collection->Size(),false,false);cg_q_full.Resize(particles.array_collection->Size(),false,false);cg_s_full.Resize(particles.array_collection->Size(),false,false);cg_r_full.Resize(particles.array_collection->Size(),false,false);
    KRYLOV_VECTOR_T pressure(pressure_full,force_dynamic_particles_list),cg_q(cg_q_full,force_dynamic_particles_list),
        cg_s(cg_s_full,force_dynamic_particles_list),cg_r(cg_r_full,force_dynamic_particles_list);
    cg_t_full.Resize(particles.array_collection->Size(),false,false);
    KRYLOV_VECTOR_T cg_t(cg_t_full,force_dynamic_particles_list);
    ARRAYS_COMPUTATIONS::Fill(pressure.v,(T)0);

    {CONJUGATE_RESIDUAL<T> cr;
    INDIRECT_ARRAY<ARRAY<T> > diagonal_preconditioner(diagonal_preconditioner_full,force_dynamic_particles_list);
    if(use_diagonal_preconditioner) divergence.v*=diagonal_preconditioner;
    T tolerance=max((T).01*system.Convergence_Norm(divergence),(T)1e-10);
    bool converged=cr.Solve(system,pressure,divergence,cg_q,cg_s,cg_t,cg_r,tolerance,0,max_cg_iterations);
    if(use_diagonal_preconditioner) pressure.v*=diagonal_preconditioner;
    LOG::Stat("divergence magnitude",sqrt(cr.residual_magnitude_squared));
    LOG::Stat("nullspace measure",cr.nullspace_measure);
    if(!converged) {std::stringstream ss;ss<<"CONJUGATE_RESIDUAL FAILED - GIVING UP"<<std::endl;LOG::filecout(ss.str());}}

    Gradient(pressure_full,gradient_full);
    gradient_full.Subset(force_dynamic_particles_list)*=particles.one_over_mass.Subset(force_dynamic_particles_list);
    Project_Vector_Field(gradient_full);
    for(int i=1;i<=force_dynamic_particles_list.m;i++){int p=force_dynamic_particles_list(i);particles.V(p)-=gradient_full(p);}
}
//#####################################################################
// Function Test_System
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Test_System()
{
    if(disable_projection) return;
    RANDOM_NUMBERS<T> random;random.Set_Seed(1823);
    Update_Position_Based_State(0,true);

    pressure_full.Resize(particles.array_collection->Size(),false,false);divergence_full.Resize(particles.array_collection->Size(),false,false);
    const ARRAY<int> &fragment_dynamic_particles=force_dynamic_particles_list,
        &fragment_particles=force_dynamic_particles_list;
    INDIRECT_ARRAY<ARRAY<T> > volumes(volumes_full,fragment_dynamic_particles);
    KRYLOV_VECTOR_T pressure(pressure_full,fragment_dynamic_particles),divergence(divergence_full,fragment_dynamic_particles);
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > X(particles.X,fragment_dynamic_particles),V(particles.V,fragment_dynamic_particles);
    INDIRECT_ARRAY<ARRAY<TV> > gradient(gradient_full,fragment_dynamic_particles);
    ARRAY<T> old_volumes(volumes);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&> V_subset=particles.V.Subset(fragment_particles);ARRAYS_COMPUTATIONS::Fill(V_subset,TV());
    INDIRECT_ARRAY<ARRAY<T>,ARRAY<int>&> pressure_subset=pressure_full.Subset(fragment_particles);ARRAYS_COMPUTATIONS::Fill(pressure_subset,T());

    T dt=(T)1e-3;
    POISSON_SYSTEM<TV,d> system(*this);

    for(int iteration=1;iteration<=10;iteration++){
        for(int p=1;p<=fragment_dynamic_particles.m;p++){
            random.Set_Seed(fragment_dynamic_particles(p)*123+20+23*iteration);
            T scale=random.Get_Uniform_Number(-(T)10,(T)10);TV direction=random.template Get_Direction<TV>();
            V(p)=exp(scale)*direction;}
        for(int i=1;i<=fragment_dynamic_particles.m;i++){int p=fragment_dynamic_particles(i);
            random.Set_Seed(p*123+21+23*iteration);
            pressure_full(p)=exp(random.Get_Uniform_Number(-(T)10,(T)10));}
        system.Project(pressure);

        T convergence_norm=ARRAYS_COMPUTATIONS::Maximum_Magnitude(V);
        if(mpi_solids) convergence_norm=mpi_solids->Reduce_Max(convergence_norm);
        {std::stringstream ss;ss<<"|V|_inf = "<<convergence_norm<<std::endl;LOG::filecout(ss.str());}
        V/=convergence_norm;

        T magnitude_squared=ARRAYS_COMPUTATIONS::Dot_Product(V,V);
        if(mpi_solids) magnitude_squared=mpi_solids->Reduce_Add(magnitude_squared);
        T magnitude=sqrt(magnitude_squared);
        {std::stringstream ss;ss<<"|V|_before = "<<magnitude<<std::endl;LOG::filecout(ss.str());}
        pressure.v/=system.Convergence_Norm(pressure);
        {std::stringstream ss;ss<<"|V| = "<<magnitude<<", |p| = "<<system.Magnitude(pressure)<<std::endl;LOG::filecout(ss.str());}

        Project_Vector_Field(particles.V);
        Negative_Divergence(particles.V,divergence_full);system.Project(divergence);
        Gradient(pressure_full,gradient_full);
        Project_Vector_Field(gradient_full);

        T gradient_dot_V=ARRAYS_COMPUTATIONS::Dot_Product(gradient,V),divergence_dot_pressure=(T)system.Inner_Product(divergence,pressure);
        if(mpi_solids) gradient_dot_V=mpi_solids->Reduce_Add(gradient_dot_V);
        {std::stringstream ss;ss<<"gradient_dot_V = "<<gradient_dot_V<<", divergence_dot_pressure = "<<divergence_dot_pressure
                                <<", relative error: "<<abs(Robust_Divide(gradient_dot_V,divergence_dot_pressure)-1)<<std::endl;LOG::filecout(ss.str());}

        ARRAY<TV> X_old(X);
        X+=dt*V;
        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(particles.X);
        Update_Position_Based_State(0,true);

        ARRAY<T> dv_full(volumes_full.m);KRYLOV_VECTOR_T dv(dv_full,fragment_dynamic_particles);
        dv.v=volumes-old_volumes;
        divergence*=dt;

        {std::stringstream ss;ss<<"|dv| = "<<system.Magnitude(dv)<<std::endl;
            ss<<"dv_dot_divergence / |divergence|^2 = "<<system.Inner_Product(dv,divergence)/system.Inner_Product(divergence,divergence)<<std::endl;LOG::filecout(ss.str());}

        // restore old positions
        X=X_old;
        Update_Position_Based_State(0,true);}

    exit(1);
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Gradient(ARRAY_VIEW<const T> p,ARRAY<TV>& gradient) const
{
    ARRAY_VIEW<T> modifiable_p(p.Const_Cast());
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(modifiable_p);
    gradient.Resize(particles.array_collection->Size(),false,false);
    INDIRECT_ARRAY<ARRAY<TV>,ARRAY<int>&> gradient_subset=gradient.Subset(force_dynamic_particles_list);ARRAYS_COMPUTATIONS::Fill(gradient_subset,TV());
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        VECTOR<int,d+1>& element=strain_measure.mesh.elements(t);
        INDIRECT_ARRAY<ARRAY_VIEW<const T>,VECTOR<int,d+1>&> p_subset=p.Subset(element);
        strain_measure.Distribute_Force(gradient,element,Bs_per_node(t)*-ARRAYS_COMPUTATIONS::Sum(p_subset));}
}
//#####################################################################
// Function Negative_Divergence
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Negative_Divergence(ARRAY_VIEW<const TV> V,ARRAY<T>& divergence) const
{
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(V.Const_Cast());
    divergence.Resize(particles.array_collection->Size(),false,false);
    INDIRECT_ARRAY<ARRAY<T>,ARRAY<int>&> divergence_subset=divergence.Subset(force_dynamic_particles_list);
    ARRAYS_COMPUTATIONS::Fill(divergence_subset,T());
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        VECTOR<int,d+1>& element=strain_measure.mesh.elements(t);
        T divergence_per_node=T_MATRIX::Inner_Product(Bs_per_node(t),strain_measure.Ds(V,t));
        divergence.Subset(element)-=divergence_per_node;}
}
//#####################################################################
// Function Diagonal_Elements
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Diagonal_Elements(ARRAY<T>& D) const
{
    ARRAY<TV> forces(particles.array_collection->Size());
    ARRAY<ARRAY<int> >& incident_elements=*strain_measure.mesh.incident_elements;
    for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        D(p)=T();
        INDIRECT_ARRAY<ARRAY<TV>,ARRAY<int>&> forces_subset=forces.Subset(node_regions(p));
        ARRAYS_COMPUTATIONS::Fill(forces_subset,TV());
        for(int j=1;j<=incident_elements(p).m;j++){int t=incident_elements(p)(j);
            strain_measure.Distribute_Force(forces,strain_measure.mesh.elements(t),-Bs_per_node(t));}
        forces.Subset(node_regions(p))*=particles.one_over_mass.Subset(node_regions(p));
        for(int j=1;j<=incident_elements(p).m;j++){int t=incident_elements(p)(j);
            D(p)-=T_MATRIX::Inner_Product(Bs_per_node(t),strain_measure.Ds(forces,t));}}
}
//#####################################################################
// Function Diagonal_Elements
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Update_Preconditioner()
{
    if(!use_diagonal_preconditioner) return;
    diagonal_preconditioner_full.Resize(particles.array_collection->Size(),false,false);
    Diagonal_Elements(diagonal_preconditioner_full);
    for(ELEMENT_ITERATOR iterator(force_dynamic_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        diagonal_preconditioner_full(p)=1/sqrt(diagonal_preconditioner_full(p));}
}
//#####################################################################
// Function Project_Vector_Field
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Project_Vector_Field(ARRAY_VIEW<TV> field) const
{
    if(!use_neumann) return;
    if(self_collision_subcycles==1) PHYSBAM_FATAL_ERROR();
    Project_All_Isolated_Clamping_Constraints(field,projection_data);
    Project_All_Clamping_Constraints(field,projection_data);
    if(repulsions && use_self_moving_projection && d==3 && (projection_data.point_face_pairs.m || projection_data.edge_edge_pairs.m))
        for(int i=1;i<=self_collision_subcycles;i++){
            TRIANGLE_REPULSIONS<TV>::Project_All_Moving_Constraints(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,field);
            Project_All_Clamping_Constraints(field,projection_data);}
}
//#####################################################################
// Function Project_All_Clamping_Constraints
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Project_All_Clamping_Constraints(ARRAY_VIEW<TV> field,const PROJECTION_DATA& data) const
{
    for(int i=1;i<=data.neumann_boundary_nodes.m;i++){int p=data.neumann_boundary_nodes(i);
        field(p)-=TV::Dot_Product(field(p),data.neumann_boundary_normals(p))*data.neumann_boundary_normals(p);}
    for(int i=1;i<=data.fixed_nodes.m;i++) field(data.fixed_nodes(i))=TV(); // TODO: generalize constructs
}
//#####################################################################
// Function Project_All_Isolated_Clamping_Constraints
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Project_All_Isolated_Clamping_Constraints(ARRAY_VIEW<TV> field,const PROJECTION_DATA& data) const
{
    for(int i=1;i<=data.neumann_boundary_nodes_isolated.m;i++){int p=data.neumann_boundary_nodes_isolated(i);
        field(p)-=TV::Dot_Product(field(p),data.neumann_boundary_normals(p))*data.neumann_boundary_normals(p);}
}
//#####################################################################
// Function Set_Neumann_Boundary_Conditions
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Set_Neumann_Boundary_Conditions(const ARRAY<COLLISION_PARTICLE_STATE<TV> >* particle_states,TRIANGLE_REPULSIONS<TV>* repulsions_input)
{
    if(repulsions_input) repulsions=repulsions_input;

    neumann_boundary_count.Resize(particles.array_collection->Size());
    // TODO: Unnecessarily expensive
    ARRAYS_COMPUTATIONS::Fill(neumann_boundary_count,0);
    // TODO: Assuming we are created new simplifies this a little.
    projection_data.point_face_pairs.Remove_All();
    projection_data.edge_edge_pairs.Remove_All();
    projection_data.neumann_boundary_nodes.Remove_All();
    // TODO: Unnecessarily expensive
    projection_data.neumann_boundary_normals.Resize(particles.array_collection->Size());ARRAYS_COMPUTATIONS::Fill(projection_data.neumann_boundary_normals,TV());

    if(particle_states && use_rigid_clamp_projection)
        for(ELEMENT_ITERATOR inner_iterator(force_dynamic_particles);inner_iterator.Valid();inner_iterator.Next()){int p=inner_iterator.Data();
            const COLLISION_PARTICLE_STATE<TV>& collision=(*particle_states)(p);
            if(collision.enforce){
                projection_data.neumann_boundary_normals(p)=collision.normal;
                if(!neumann_boundary_count(p))
                    projection_data.neumann_boundary_nodes.Append(p);
                neumann_boundary_count(p)=1;}}

    if(repulsions && use_self_moving_projection){
        repulsions->Set_Collision_Pairs(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,projection_data.point_face_pairs,projection_data.edge_edge_pairs,(T)2); // TODO: Fix me.
        for(int i=1;i<=projection_data.point_face_pairs.m;i++) neumann_boundary_count.Subset(projection_data.point_face_pairs(i).nodes)+=1;
        for(int i=1;i<=projection_data.edge_edge_pairs.m;i++) neumann_boundary_count.Subset(projection_data.edge_edge_pairs(i).nodes)+=1;}

    int j=0;
    for(int i=1;i<=projection_data.neumann_boundary_nodes.m;i++){int p=projection_data.neumann_boundary_nodes(i);
        if(neumann_boundary_count(p)==1) projection_data.neumann_boundary_nodes_isolated.Append(p);
        else projection_data.neumann_boundary_nodes(++j)=p;}
    projection_data.neumann_boundary_nodes.m=j;

    std::stringstream ss;ss<<"moving constraints: "<<projection_data.point_face_pairs.m<<std::endl;LOG::filecout(ss.str());
}
//#####################################################################
// Function Max_Relative_Velocity_Error
//#####################################################################
template<class TV,int d> typename TV::SCALAR INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Max_Relative_Velocity_Error()
{
    Update_Position_Based_State(0,true);
    T max_error=0;
    
    for(int i=1;i<=force_dynamic_particles_list.m;i++){int p=force_dynamic_particles_list(i);
        if(rest_volumes_full(p)){T error=(volumes_full(p)-rest_volumes_full(p))/rest_volumes_full(p);max_error=max(max_error,abs(error));}}
    if(mpi_solids) max_error=mpi_solids->Reduce_Max_Global(max_error);
    return max_error;
}
//#####################################################################
// Function Save_Volumes
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Save_Volumes()
{
    Update_Position_Based_State(0,true);
    saved_volumes_full=volumes_full;
}
//#####################################################################
// Function Check_Improvement
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Check_Improvement()
{
    Update_Position_Based_State(0,true);
    T ave_improve=0,max_improve=-FLT_MAX,min_improve=FLT_MAX,old_total_volume_accumulated=0,new_total_volume_accumulated=0;
    int better=0,worse=0,same=0;
    const INDIRECT_ARRAY<const ARRAY<T> > volumes(volumes_full.Subset(force_dynamic_particles_list));
    const INDIRECT_ARRAY<const ARRAY<T> > saved_volumes(saved_volumes_full.Subset(force_dynamic_particles_list));
    const INDIRECT_ARRAY<const ARRAY<T> > rest_volumes(rest_volumes_full.Subset(force_dynamic_particles_list));
    T new_total_volume=ARRAYS_COMPUTATIONS::Sum(volumes);
    if(mpi_solids) new_total_volume=mpi_solids->Reduce_Add_Global(new_total_volume);
    new_total_volume_accumulated+=new_total_volume;
    T old_total_volume=ARRAYS_COMPUTATIONS::Sum(saved_volumes);
    if(mpi_solids) old_total_volume=mpi_solids->Reduce_Add_Global(old_total_volume);
    old_total_volume_accumulated+=old_total_volume;
    if(volumes.Size()!=rest_volumes.Size() || volumes.Size()!=saved_volumes.Size()){std::stringstream ss;ss<<"Volume arrays have different sizes!"<<std::endl;PHYSBAM_FATAL_ERROR();LOG::filecout(ss.str());}
    for(int i=1;i<=volumes.Size();i++) if(rest_volumes(i)){
        T old_v=abs((saved_volumes(i)-rest_volumes(i))/rest_volumes(i)),new_v=abs((volumes(i)-rest_volumes(i))/rest_volumes(i)),diff=new_v-old_v;
        ave_improve+=diff;
        max_improve=max(max_improve,diff);
        min_improve=min(min_improve,diff);
        if(abs(new_v-old_v)<1e-2) same++;
        else if(new_v<old_v) better++;
        else worse++;}

     if(mpi_solids) max_improve=mpi_solids->Reduce_Max_Global(max_improve);
     if(mpi_solids) min_improve=mpi_solids->Reduce_Min_Global(min_improve);
     if(mpi_solids) ave_improve=mpi_solids->Reduce_Add_Global(ave_improve);
     if(mpi_solids) better=mpi_solids->Reduce_Add_Global(better);
     if(mpi_solids) same=mpi_solids->Reduce_Add_Global(same);
     if(mpi_solids) worse=mpi_solids->Reduce_Add_Global(worse);
    ave_improve/=particles.array_collection->Size();
    std::stringstream ss;ss<<"INCOMP STAT total volume error: "<<abs(old_total_volume_accumulated-total_rest_volume)<<" => "<<abs(new_total_volume_accumulated-total_rest_volume)<<std::endl;
    ss<<"INCOMP STAT better: "<<better<<"    same: "<<same<<"    worse: "<<worse<<std::endl;
    ss<<"INCOMP STAT improvements -   min: "<<min_improve<<"    ave: "<<ave_improve<<"    max: "<<max_improve<<std::endl;LOG::filecout(ss.str());
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    strain_measure.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<float,2>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<float,3>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<double,2>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<double,3>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<double,3>,3>;
#endif
}
