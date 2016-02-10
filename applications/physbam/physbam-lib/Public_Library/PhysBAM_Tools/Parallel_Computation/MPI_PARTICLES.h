//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_PARTICLES
//#####################################################################
#ifndef __MPI_PARTICLES__
#define __MPI_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>

namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;
template<class TV> class RANGE;
template<class TV> class GRID;

template<class T_GRID>
class MPI_PARTICLES:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef MPI_GRID<T_GRID> T_MPI_GRID;
public:
    T_MPI_GRID& mpi_grid;
    ARRAY<RANGE<TV_INT> > send_regions;
    ARRAY<ARRAY<int> > send_particles;
    ARRAY<RANGE<VECTOR<int,1> > > recv_particle_range;
    ARRAY<int> neighbor_ranks;
    ARRAY<TV_INT> neighbor_directions;
    ARRAY<int> send_particle_counts,recv_particle_counts;
    
    MPI_PARTICLES(T_MPI_GRID& mpi_grid_input);
    ~MPI_PARTICLES();

    bool Inside_Local_Domain(const TV& position,const RANGE<TV>& local_domain)const
    {
        bool flag=true;
        for(int axis=1;axis<=TV::m;axis++) for(int axis_side=1;axis_side<=2;axis_side++){//half open if it is not global boundary
            if(axis_side==1 && position(axis)<local_domain.min_corner(axis)){flag=false;break;}
            else if(axis_side==2){
                if(local_domain.max_corner(axis)!=mpi_grid.global_grid.domain.max_corner(axis) && position(axis)>=local_domain.max_corner(axis)){flag=false;break;}
                if(local_domain.max_corner(axis)==mpi_grid.global_grid.domain.max_corner(axis) && position(axis)>local_domain.max_corner(axis)){flag=false;break;}}}
        return flag;
    }

    bool Inside_Local_Domain(const TV& position)const
    {
        bool flag=true;
        for(int axis=1;axis<=TV::m;axis++) for(int axis_side=1;axis_side<=2;axis_side++)
            if(mpi_grid.Neighbor(axis,axis_side)){ //half open if it has neighbor
                if(axis_side==1 && position(axis)<mpi_grid.local_grid.domain.min_corner(axis)){flag=false;break;}
                else if(axis_side==2 && position(axis)>=mpi_grid.local_grid.domain.max_corner(axis)){flag=false;break;}}
            else{
                if(axis_side==1 && position(axis)<mpi_grid.local_grid.domain.min_corner(axis)){flag=false;break;}
                else if(axis_side==2 && position(axis)>mpi_grid.local_grid.domain.max_corner(axis)){flag=false;break;}}
        return flag;
    }

    template<class T_ARRAY> void Update_Boundary_Particle_Mapping(T_ARRAY& target_array,ARRAY_VIEW<TV> positions,const int bandwidth,const ARRAY<int>* indices=0,const bool include_ghost_regions=false,const bool include_corners=true);
    template<class T2> void Exchange_Boundary_Particle_Values(ARRAY_VIEW<T2> values) const;
    template<class T2> void Backwards_Exchange_Boundary_Particle_Values(ARRAY_VIEW<T2> send_values,ARRAY<ARRAY_VIEW<T2>*> recv_values) const;
    void Remove_Particles_Outside_Local_Domain(ARRAY_COLLECTION& array_collection,ARRAY_VIEW<TV> positions,const ARRAY<int>* indices=0);
    void Synchronize_Flag(bool& flag) const;
    void Output_Send_And_Recv_Indices(bool verbose=false)const;
};
}
#endif
