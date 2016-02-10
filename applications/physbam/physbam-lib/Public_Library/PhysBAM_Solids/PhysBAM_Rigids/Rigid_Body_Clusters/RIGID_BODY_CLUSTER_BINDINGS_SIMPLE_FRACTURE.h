//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE__
#define __RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE__
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/UNDIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>

namespace PhysBAM{

template<class TV>
class RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE:public RIGID_BODY_CLUSTER_BINDINGS<TV>::CALLBACKS
{
    typedef typename TV::SCALAR T;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& bindings;
    ARRAY<int> parents_to_rebuild;
public:
    struct FRACTURE_DATA
    {
        ARRAY<VECTOR<RIGID_CLUSTER_CONSTITUENT_ID,2> > connections;
        ARRAY<T> restlengths;
    };

    T allowed_strain;
    HASHTABLE<VECTOR<int,2>,T> allowed_strains,decay,decay_rate;
    UNDIRECTED_GRAPH<int,int>* graph;
    T local_dt;

    RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_BODY_CLUSTER_BINDINGS<TV>& bindings)
        :rigid_body_collection(rigid_body_collection),bindings(bindings),allowed_strain((T)FLT_MAX)
    {}

    HASHTABLE<int,FRACTURE_DATA> fracture_data; // from cluster parent to FRACTURE_DATA

    void Pre_Advance_Unclustered(const T dt,const T time){}
    void Post_Advance_Unclustered(const T dt,const T time){local_dt=dt;}

    void Create_Cluster(const int parent)
    {Initialize_Strain(parent,fracture_data.Get_Or_Insert(parent));}

    void Destroy_Cluster(const int parent)
    {fracture_data.Delete(parent);}

    void Initialize_Strain(const int parent,FRACTURE_DATA& data)
    {typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*bindings.reverse_bindings.Get(parent);
    // compute restlengths  
    for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++) for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<i;j++){
        RIGID_BODY<TV> &child_1=rigid_body_collection.Rigid_Body(cluster.children(i)),&child_2=rigid_body_collection.Rigid_Body(cluster.children(j));
        data.connections.Append(VECTOR<RIGID_CLUSTER_CONSTITUENT_ID,2>(i,j));
        data.restlengths.Append((child_2.Frame().t-child_1.Frame().t).Magnitude());
        {std::stringstream ss;ss<<"restlength on (local coords "<<i<<","<<j<<")  "<<cluster.children(i)<<","<<cluster.children(j)<<" is "<<data.restlengths(data.restlengths.Size())<<std::endl;LOG::filecout(ss.str());}
    }}

    void Find_Weakest_Links(int root,T min_strain,HASHTABLE<int>& visited,ARRAY<int>& edges)
    {if(!visited.Contains(root)) visited.Insert(root);
    ARRAY<int> next_nodes;
    for(int i=1;i<=graph->Adjacent_Edges(root).m;i++){
        PAIR<int,int> edge=graph->Edges(graph->Adjacent_Edges(root)(i));
        int other_node=edge.x;if(other_node==root) other_node=edge.y;
        if(visited.Contains(other_node)) continue;
        T strain=allowed_strains.Get_Default(VECTOR<int,2>(root,other_node).Sorted(),allowed_strain);
        if(strain<=min_strain){
            if(strain<min_strain){min_strain=strain;edges.Remove_All();next_nodes.Remove_All();}
            edges.Append(graph->Adjacent_Edges(root)(i));
            next_nodes.Append(other_node);}}
    for(int i=1;i<=next_nodes.m;i++) Find_Weakest_Links(next_nodes(i),min_strain,visited,edges);}

    void Compute_New_Clusters_Based_On_Unclustered_Strain()
    {parents_to_rebuild.Remove_All();
    ARRAY<int,int> components;if(graph) graph->Connected_Components(components);
    for(typename HASHTABLE<int,FRACTURE_DATA>::ITERATOR iterator(fracture_data);iterator.Valid();iterator.Next()){
        ARRAY<int> remove_connections;HASHTABLE<int> visited;bool need_rebuild=false;
        FRACTURE_DATA& data=iterator.Data();
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*bindings.reverse_bindings.Get(iterator.Key());
        for(int i=data.connections.m;i>=1;i--){
            const VECTOR<RIGID_CLUSTER_CONSTITUENT_ID,2>& edge=data.connections(i);
            T rl=data.restlengths(i);
            RIGID_BODY<TV> &child_1=rigid_body_collection.Rigid_Body(cluster.children(edge[1])),&child_2=rigid_body_collection.Rigid_Body(cluster.children(edge[2]));
            VECTOR<int,2> hash_index=VECTOR<int,2>(cluster.children(edge[1]),cluster.children(edge[2])).Sorted();
            T strain=abs((child_2.Frame().t-child_1.Frame().t).Magnitude()/rl-(T)1);
            T& local_decay=decay.Get_Or_Insert(hash_index,0);local_decay+=local_dt*decay_rate.Get_Or_Insert(hash_index,0);
            strain+=local_decay;
            {std::stringstream ss;ss<<"strain between "<<cluster.children(edge[1])<<","<<cluster.children(edge[2])<<" is "<<strain<<std::endl;LOG::filecout(ss.str());}
            T allowed_strain_local=allowed_strains.Get_Default(VECTOR<int,2>(cluster.children(edge[1]),cluster.children(edge[2])).Sorted(),allowed_strain);
            if(strain>allowed_strain_local){need_rebuild=true;
                if(graph){
                    ARRAY<int> break_connections;int root=cluster.children(edge[1]);
                    if(!visited.Contains(root)){
                        Find_Weakest_Links(root,allowed_strain_local,visited,break_connections);
                        for(int i=1;i<=break_connections.m;i++){
                            VECTOR<int,2> edge=VECTOR<int,2>(graph->Edges(break_connections(i)).x,graph->Edges(break_connections(i)).y).Sorted();
                            graph->Remove_Edge(break_connections(i));
                            for(int j=data.connections.m;j>=1;j--){
                                VECTOR<int,2> tmp_edge(Value(data.connections(j)(1)),Value(data.connections(j)(2)));
                                if(edge==tmp_edge.Sorted()) remove_connections.Append(j);}}}}
                else data.connections.Remove_Index_Lazy(i);}}
        if(need_rebuild) parents_to_rebuild.Append(iterator.Key());
        Sort(remove_connections);
        for(int i=remove_connections.m;i>=1;i--) data.connections.Remove_Index_Lazy(remove_connections(i));}}

    bool Create_New_Clusters()
    {for(int i_dummy=1;i_dummy<=parents_to_rebuild.m;i_dummy++){int parent=parents_to_rebuild(i_dummy);
        FRACTURE_DATA& data=fracture_data.Get(parent);
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*bindings.reverse_bindings.Get(parent);
        UNDIRECTED_GRAPH<RIGID_CLUSTER_CONSTITUENT_ID,int> undirected(cluster.children.Size());
        for(int i=1;i<=data.connections.m;i++) undirected.Add_Edge(data.connections(i).x,data.connections(i).y,undirected.Last_Edge()+1);
        ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> components;
        int number_components=undirected.Connected_Components(components);
        if(number_components==1) continue;
        ARRAY<ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> > new_cluster_constituents(number_components);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++)
            new_cluster_constituents(components(i)).Append(cluster.children(i));
        bindings.Delete_Binding(parent);
        for(int i=1;i<=new_cluster_constituents.m;i++) if(new_cluster_constituents(i).Size()>RIGID_CLUSTER_CONSTITUENT_ID(1)) bindings.Add_Binding(new_cluster_constituents(i));}
    //for(typename HASHTABLE<int,typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER>::ITERATOR iterator(bindings.reverse_bindings);iterator.Valid();iterator.Next()){
    //    {std::stringstream ss;ss<<"cluster parent "<<iterator.Key()<<std::endl;LOG::filecout(ss.str());} 
    //    for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<=iterator.Data().children.Size();j++)
    //        {std::stringstream ss;ss<<"   constituent "<<iterator.Data().children(j)<<std::endl;LOG::filecout(ss.str());}}
    //for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
    //    {std::stringstream ss;ss<<"rigid "<<i<<" -> "<<(i<=bindings.binding_index.m?bindings.binding_index(i).x:0)<<std::endl;LOG::filecout(ss.str());}}
    return parents_to_rebuild.m>0;}

//#####################################################################
};
}
#endif
