//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_CONNECTIVITY
//#####################################################################
#ifndef __PARTICLE_CONNECTIVITY__
#define __PARTICLE_CONNECTIVITY__

#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class PARTICLES;

template<class TV>
class PARTICLE_CONNECTIVITY
{
    typedef unsigned char T_RANK;
public:
    const int particles_number;
private:
    const RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    SPARSE_UNION_FIND<> union_find;

    friend class SOLID_BODY_COLLECTION<TV>;
public:

    PARTICLE_CONNECTIVITY(const PARTICLES<TV>& particles,const RIGID_BODY_COLLECTION<TV>& rigid_body_collection);

    ~PARTICLE_CONNECTIVITY()
    {}

    int Size() const
    {return union_find.Size();}

    void Clear_Connectivity()
    {union_find.Clear_Connectivity();}

    int Find(const int i) const
    {return union_find.Find(i);}

    int Union(const int i,const int j)
    {if(Exclude_Particle(i) || Exclude_Particle(j)) return 0;return union_find.Union(i,j);}

    template<class T_ARRAY>
    int Union(const T_ARRAY& array)
    {int root=0;typename T_ARRAY::ELEMENT i(1);for(;i<=array.Size();i++) if(!Exclude_Particle(array(i))){root=Find(array(i));break;}if(!root) return 0;
    for(;i<=array.Size();i++) if(!Exclude_Particle(array(i))) union_find.Union(root,array(i));return union_find.Find(root);}

    void Merge(const PARTICLE_CONNECTIVITY<TV>& particle_connectivity)
    {union_find.Merge(particle_connectivity.union_find);}

    // Okay for map to yield invalid indices for isolated elements
    template<class T_ARRAY>
    void Mapped_Merge(const PARTICLE_CONNECTIVITY<TV>& particle_connectivity,const T_ARRAY& map)
    {union_find.Mapped_Merge(particle_connectivity.union_find);}

    void Forest_Edges(ARRAY<PAIR<int,int> >& pairs) const
    {union_find.Forest_Edges(pairs);}

    void Merge_Forest_Edges(const ARRAY<PAIR<int,int> >& pairs)
    {union_find.Merge_Forest_Edges(pairs);}

    void Union_All()
    {Union(IDENTITY_ARRAY<>(Size()));}

    void Union_All_Deformable_Particles()
    {Union(IDENTITY_ARRAY<>(particles_number));}

//#####################################################################
    bool Exclude_Particle(const int i) const;
    void Union_All_Rigid_Body_Particles();
//#####################################################################
};
}
#endif
