//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __PARTICLE_LEVELSET_IMPLICIT_OBJECT__
#define __PARTICLE_LEVELSET_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class PARTICLE_LEVELSET_IMPLICIT_OBJECT:public LEVELSET_IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename GRID<TV>::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename LEVELSET_POLICY<GRID<TV> >::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID<TV>::BLOCK T_BLOCK;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<bool> >::TYPE T_ARRAYS_ARRAY_BOOL;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::LINEAR_INTERPOLATION_SCALAR::template REBIND<TV>::TYPE T_LINEAR_INTERPOLATION_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
public:
    typedef LEVELSET_IMPLICIT_OBJECT<TV> BASE;
    using BASE::levelset;using BASE::need_destroy_data;

    T_PARTICLE_LEVELSET& particle_levelset;
    T_ARRAYS_BOOL particle_influence;
    GRID<TV> p_grid;
    T_ARRAYS_ARRAY_BOOL influencing_positive_particles;
    T_ARRAYS_ARRAY_BOOL influencing_negative_particles;
    T_LINEAR_INTERPOLATION_VECTOR linear_interpolation;

    PARTICLE_LEVELSET_IMPLICIT_OBJECT(T_PARTICLE_LEVELSET& particle_levelset_input);
    virtual ~PARTICLE_LEVELSET_IMPLICIT_OBJECT()
    {}

    void Precompute_Cell_Particle_Influence()
    {influencing_positive_particles.Resize(levelset.grid.Domain_Indices(3));influencing_negative_particles.Resize(levelset.grid.Domain_Indices(3));particle_influence.Fill(false);
    Precompute_Cell_Particle_Influence(particle_levelset.positive_particles,influencing_positive_particles,1);
    Precompute_Cell_Particle_Influence(particle_levelset.negative_particles,influencing_negative_particles,-1);}

    void Precompute_Cell_Particle_Influence(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_ARRAY_BOOL& influencing_particles,int sign)
    {T one_over_radius_multiplier=-(T)sign;
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){
        TV_INT block_index=iterator.Node_Index();
        ARRAY<bool>& influence=influencing_particles(block_index);
        if(particles(block_index)){
            PARTICLE_LEVELSET_PARTICLES<TV>& block_particles=*particles(block_index);
            influence.Resize(block_particles.array_collection->Size());
            for(int k=1;k<=block_particles.array_collection->Size();k++){
                if(levelset.Phi(block_particles.X(k))*one_over_radius_multiplier>block_particles.radius(k)){
                    influence(k)=true;
                    for(int cell=1;cell<=GRID<TV>::number_of_cells_per_node;cell++)
                        particle_influence(p_grid.Node_Cell_Index(block_index,cell))=true;}
                else influence(k)=false;}}
        else influence.Resize(0);}}

    void Particle_Phi_Value(const PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles,const ARRAY<bool>& influence,int sign,const TV& location,T& phi) const
    {for(int k=1;k<=cell_particles.array_collection->Size();k++){
        if(!influence(k)) continue;
        T radius_minus_sign_phi=cell_particles.radius(k)-sign*phi;
        if(radius_minus_sign_phi > 0){
            T distance_squared=(location-cell_particles.X(k)).Magnitude_Squared();
            if(distance_squared < sqr(radius_minus_sign_phi)){
                phi=sign*(cell_particles.radius(k)-sqrt(distance_squared));}}}}
    
    T operator()(const TV& location) const PHYSBAM_OVERRIDE
    {TV_INT cell_index=p_grid.Clamped_Index(location);
    T phi=levelset.Phi(location);
    if(particle_influence(cell_index)){
        static TV_INT nodes[GRID<TV>::number_of_nodes_per_cell];
        T phi_plus=phi,phi_minus=phi;
        levelset.grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell_index,nodes);
        for(int node_index=0;node_index<GRID<TV>::number_of_nodes_per_cell;node_index++){
            TV_INT node=nodes[node_index];
            if(!particle_influence(node)) continue;
            if(particle_levelset.positive_particles(node)) Particle_Phi_Value(*particle_levelset.positive_particles(node),influencing_positive_particles(node),1,location,phi_plus);
            if(particle_levelset.negative_particles(node)) Particle_Phi_Value(*particle_levelset.negative_particles(node),influencing_negative_particles(node),-1,location,phi_minus);}
        phi=minmag(phi,phi_plus,phi_minus);}
    return phi;}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("PARTICLE_LEVELSET_IMPLICIT_OBJECT<T,VECTOR<T,%d> >",TV::dimension);}
    
//###########################################################################
    static PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>* Create();
//###########################################################################
};

}
#endif
