//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLED_EVOLUTION
//#####################################################################
#ifndef __SOLID_FLUID_COUPLED_EVOLUTION__
#define __SOLID_FLUID_COUPLED_EVOLUTION__

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class T_GRID> class FLUIDS_PARAMETERS_UNIFORM;
template<class T_GRID> class INCOMPRESSIBLE_FLUID_CONTAINER;
template<class T_GRID> class POISSON_COLLIDABLE_UNIFORM;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV_input>
class SOLID_FLUID_COUPLED_EVOLUTION:public NEWMARK_EVOLUTION<TV_input>
{
    typedef TV_input TV;typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;typedef GRID<TV> T_GRID;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef ARRAY<PAIR<int,T> > FACE_WEIGHT_ELEMENTS;typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS::template REBIND<FACE_WEIGHT_ELEMENTS*>::TYPE T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> > >::TYPE T_ARRAYS_STRUCTURE_SIMPLEX_LIST;
    typedef typename MESH_POLICY<TV::dimension-1>::MESH T_BOUNDARY_MESH;
    typedef typename MESH_POLICY<TV::dimension>::MESH T_MESH;typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension>::OBJECT T_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_BOUNDARY_SIMPLEX;typedef VECTOR<int,TV::dimension> T_BOUNDARY_ELEMENT;typedef VECTOR<int,TV::dimension+1> T_ELEMENT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX T_SIMPLEX;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename TV::SPIN T_SPIN;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_INERTIA_TENSOR;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension-1>::OBJECT T_THIN_SHELL;typedef typename T_THIN_SHELL::MESH T_THIN_SHELL_MESH;
    typedef VECTOR<int,TV::dimension> T_THIN_SHELL_ELEMENT;typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_THIN_SHELL_SIMPLEX;
    typedef typename MATRIX_POLICY<TV>::DIAGONAL_MATRIX T_DIAGONAL_MATRIX;
    typedef FIELD_PROJECTOR<SPARSE_MATRIX_ENTRY<T>,int,&SPARSE_MATRIX_ENTRY<T>::j> T_PROJECTED_INDEX;
    typedef FIELD_PROJECTOR<SPARSE_MATRIX_ENTRY<T>,T,&SPARSE_MATRIX_ENTRY<T>::a> T_PROJECTED_VALUE;
protected:
    typedef NEWMARK_EVOLUTION<TV> BASE;
    using BASE::solid_body_collection;using BASE::solids_parameters;using BASE::F_full;using BASE::rigid_F_full;using BASE::R_full;using BASE::rigid_R_full;using BASE::S_full;
    using BASE::rigid_S_full;using BASE::B_full;using BASE::rigid_B_full;using BASE::repulsions;using BASE::rigid_deformable_collisions;using BASE::Initialize_World_Space_Masses;
    using BASE::world_space_rigid_mass_inverse;using BASE::world_space_rigid_mass;using BASE::solids_evolution_callbacks;
    using BASE::X_save;using BASE::rigid_X_save;using BASE::rigid_rotation_save;using BASE::V_save;using BASE::rigid_velocity_save;using BASE::rigid_angular_momentum_save;

    static const int rows_per_rigid_body=TV::dimension+T_SPIN::dimension;

public:
    int rigid_body_count;
    bool print_matrix_rhs_and_solution;
protected:
    ARRAY<int> kinematic_rigid_bodies;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies;
    T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS dual_cell_weights;
    T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS rigid_body_dual_cell_weights;
    T_FACE_ARRAYS_SCALAR dual_cell_fluid_volume;
    T_FACE_ARRAYS_BOOL dual_cell_contains_solid;
    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters;
    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters;
    ARRAY<T_DIAGONAL_MATRIX> nodal_fluid_mass;
    ARRAY<T_DIAGONAL_MATRIX> rigid_body_fluid_mass;
    ARRAY<TV> rigid_body_updated_center_of_mass;
    ARRAY<T_INERTIA_TENSOR> rigid_body_fluid_inertia;
    ARRAY<TV> ar_full,z_full,zaq_full;ARRAY<TWIST<TV> > rigid_ar_full,rigid_z_full,rigid_zaq_full; // extra vectors for conjugate residual
    T_FACE_ARRAYS_SCALAR solid_projected_face_velocities_star;

    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* Get_Poisson()
    {return (fluids_parameters.compressible ? fluids_parameters.euler->euler_projection.elliptic_solver : fluids_parameters.incompressible->projection.poisson_collidable);}

    T_FACE_ARRAYS_SCALAR& Get_Face_Velocities()
    {return fluids_parameters.compressible ? fluids_parameters.euler->euler_projection.face_velocities:incompressible_fluid_container.face_velocities;}

    const GRID<TV>& Get_Grid()
    {return (fluids_parameters.compressible ? fluids_parameters.euler->grid : fluids_parameters.incompressible->projection.p_grid);}

public:
    
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > J_deformable;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > J_rigid;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > J_rigid_kinematic;
    ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array;

    T_ARRAYS_SCALAR& Get_Pressure()
    {return (fluids_parameters.compressible ? fluids_parameters.euler->euler_projection.p : fluids_parameters.incompressible->projection.p);}
    
    bool Negative(const GRID<TV>& grid,const int axis,const TV_INT& face_index,const T_ARRAYS_SCALAR& phi)
    {// this check is potentially expensive! leave inefficient for the moment
    TV_INT cells[GRID<TV>::number_of_cells_per_node];
    for(int node=1;node<=GRID<TV>::number_of_nodes_per_face;node++){
        grid.Cells_Neighboring_Node(grid.Face_Node_Index(axis,face_index,node),cells);
        for(int c=0;c<GRID<TV>::number_of_cells_per_node;c++) if(phi(cells[c])<=0) return true;}
    return false;}

    bool Simulate_Fluids() const
    {return (solids_fluids_parameters.mpi_solid_fluid || fluids_parameters.simulate) && (fluids_parameters.smoke || fluids_parameters.fire || fluids_parameters.water || fluids_parameters.sph || fluids_parameters.compressible);}

    bool Simulate_Solids() const
    {return ((solids_fluids_parameters.mpi_solid_fluid || solid_body_collection.deformable_body_collection.simulate) && solid_body_collection.deformable_body_collection.particles.array_collection->Size()) || ((solids_fluids_parameters.mpi_solid_fluid || solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies) && solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());}

    SOLID_FLUID_COUPLED_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,
        FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_input,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container_input,SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters);
    virtual ~SOLID_FLUID_COUPLED_EVOLUTION();

//#####################################################################
    T Get_Density_At_Face(const int axis,const TV_INT& face_index);
    void Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update) PHYSBAM_OVERRIDE;
    void Transfer_Momentum_And_Set_Boundary_Conditions(const T time,GENERALIZED_VELOCITY<TV>* B=0);
    void Set_Dirichlet_Boundary_Conditions(const T time);
    void Compute_W(const T current_position_time);
    void Compute_Coupling_Terms_Deformable(const T_ARRAYS_INT& cell_index_to_matrix_index,const ARRAY<INTERVAL<int> >& interior_regions,const int colors);
    void Compute_Coupling_Terms_Rigid(const T_ARRAYS_INT& cell_index_to_matrix_index,const ARRAY<INTERVAL<int> >& interior_regions,const int colors);
    void Add_Nondynamic_Solids_To_Right_Hand_Side(ARRAY<VECTOR_ND<T> >& right_hand_side,const ARRAY<INTERVAL<int> >& interior_regions,const int colors);
    void Apply_Pressure(const T dt,const T time);
    void Average_Solid_Projected_Face_Velocities_For_Energy_Update(const T_FACE_ARRAYS_SCALAR& solid_projected_face_velocities_star,const T_FACE_ARRAYS_SCALAR& solid_projected_face_velocities_np1,T_FACE_ARRAYS_SCALAR& face_velocities);
    void Apply_Solid_Boundary_Conditions(const T time,const bool use_pseudo_velocities,T_FACE_ARRAYS_SCALAR& face_velocities);
    void Set_Neumann_and_Dirichlet_Boundary_Conditions(const T time);
//#####################################################################
};
}
#endif
