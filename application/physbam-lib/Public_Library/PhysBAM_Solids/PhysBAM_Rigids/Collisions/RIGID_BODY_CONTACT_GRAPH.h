//#####################################################################
// Copyright 2003-2008, Ronald Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_CONTACT_GRAPH
//##################################################################### 
#ifndef __RIGID_BODY_CONTACT_GRAPH__
#define __RIGID_BODY_CONTACT_GRAPH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_CONTACT_GRAPH
{
public:
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles;
    DIRECTED_GRAPH<int> directed_graph; // nodes point to objects which are "above" it

    RIGID_BODY_CONTACT_GRAPH(RIGID_BODY_PARTICLES<TV>& rigid_body_particles_input) 
        :rigid_body_particles(rigid_body_particles_input),directed_graph(rigid_body_particles.array_collection->Size())
    {}

    void Reset()
    {directed_graph.Reset();}

    void Initialize() // Call when number of rigid bodies has changed
    {directed_graph.Initialize(rigid_body_particles.array_collection->Size());}

    void Add_Edge(const int body_below,const int body_above) 
    {directed_graph.Add_Edge(body_below,body_above);}

    int Number_Of_Levels() const
    {return directed_graph.Number_Of_Levels();}

    void Print() const
    {{std::stringstream ss;ss<<"CONTACT_GRAPH:"<<std::endl;LOG::filecout(ss.str());}
    for(int i(1);i<=rigid_body_particles.array_collection->Size();i++){
        {std::stringstream ss;ss<<"CONTACT_GRAPH\""<<rigid_body_particles.Rigid_Body(i).name<<"\" (LEVEL = "<<directed_graph.Level_Of_Node(i)<<"): ";LOG::filecout(ss.str());}
        if(directed_graph.Parents(i).m>0){std::stringstream ss;ss<<" IS ABOVE=";RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Parents(i));LOG::filecout(ss.str());}
        if(directed_graph.Children(i).m>0){std::stringstream ss;ss<<" IS BELOW=";RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Children(i));LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<std::endl;LOG::filecout(ss.str());}}{std::stringstream ss;ss<<std::endl;LOG::filecout(ss.str());}
    for(int i=1;i<=directed_graph.Number_Of_Levels();i++){
        {std::stringstream ss;ss<<"CONTACT_GRAPH Level "<<i<<": ";LOG::filecout(ss.str());}RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Nodes_In_Level(i));{std::stringstream ss;ss<<std::endl;LOG::filecout(ss.str());}}{std::stringstream ss;ss<<std::endl;LOG::filecout(ss.str());}}

    void Print_Statistics(const ARRAY<ARRAY<int> >& contact_pairs_for_level,const bool verbose=false)
    {int max_level_size=0,max_level_index=0,max_pairs_in_level=0,max_pairs_in_level_index=0,number_levels=Number_Of_Levels();int total=0;
    for(int i=1;i<=number_levels;i++) {
        int level_size=directed_graph.Nodes_In_Level(i);total+=level_size;
        if(level_size > max_level_size){max_level_size=level_size;max_level_index=i;}
        if(contact_pairs_for_level(i).m > max_pairs_in_level){max_pairs_in_level=contact_pairs_for_level(i).m;max_pairs_in_level_index=i;}}
    assert(total==rigid_body_particles.array_collection->Size());
    {std::stringstream ss;ss<<"Contact graph statistics: levels="<<number_levels<<", max_level_size="<<max_level_size<<", max_pairs_in_level="<<max_pairs_in_level<<std::endl;LOG::filecout(ss.str());}
    if(verbose){
        if(max_level_index > 0){
            {std::stringstream ss;ss<<"Bodies in biggest level (level="<<max_level_index<<"): ";LOG::filecout(ss.str());}
            RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Nodes_In_Level(max_level_index));{std::stringstream ss;ss<<std::endl;LOG::filecout(ss.str());}}
        if(max_pairs_in_level_index > 0){
            {std::stringstream ss;ss<<"Pairs in biggest level (level="<<max_pairs_in_level_index<<"): ";LOG::filecout(ss.str());}
            RIGID_BODY<TV>::Print_Pairs(rigid_body_particles,contact_pairs_for_level(max_pairs_in_level_index));{std::stringstream ss;ss<<std::endl;LOG::filecout(ss.str());}}}}

//#####################################################################
};   
}
#endif

