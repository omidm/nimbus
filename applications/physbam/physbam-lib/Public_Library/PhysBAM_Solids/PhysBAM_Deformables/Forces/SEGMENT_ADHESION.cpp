//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_ADHESION
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/EDGE_EDGE_ADHESION_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/EDGE_EDGE_INITIAL_CULL_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_ADHESION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;

// TODO: use worst case CFL
// TODO: add implicit velocity independent for unconditional stability

//#####################################################################
// Function SEGMENT_ADHESION
//#####################################################################
template<class TV> SEGMENT_ADHESION<TV>::
SEGMENT_ADHESION(PARTICLES<TV>& particles,SEGMENT_MESH& mesh,ARRAY<HAIR_ID>& particle_to_spring_id,HASHTABLE<VECTOR<int,4> >& intersecting_edge_edge_pairs)
    :DEFORMABLES_FORCES<TV>(particles),mpi_solids(0),mesh(mesh),curve(mesh,dynamic_cast<GEOMETRY_PARTICLES<TV>&>(particles)),restlength((T).001),max_connections(5),
    particle_to_spring_id(particle_to_spring_id),internal_curve(internal_mesh,particles),external_curve(external_mesh,particles),springs(new T_SPRING_HASH()),
    segments_with_springs(mesh.elements.m),intersecting_edge_edge_pairs(intersecting_edge_edge_pairs)
{
    for(int i=1;i<=mesh.elements.m;i++) segments_with_springs(i).x.Preallocate(max_connections);
}
//#####################################################################
// Function SEGMENT_ADHESION
//#####################################################################
template<class TV> SEGMENT_ADHESION<TV>::
SEGMENT_ADHESION(PARTICLES<TV>& particles,SEGMENT_MESH& mesh,ARRAY<HAIR_ID>& particle_to_spring_id)
    :DEFORMABLES_FORCES<TV>(particles),mpi_solids(0),mesh(mesh),curve(mesh,particles),restlength((T).001),max_connections(5),
    particle_to_spring_id(particle_to_spring_id),internal_curve(internal_mesh,particles),external_curve(external_mesh,particles),springs(new T_SPRING_HASH()),
    segments_with_springs(mesh.elements.m),intersecting_edge_edge_pairs(default_intersecting_edge_edge_pairs)
{
    for(int i=1;i<=mesh.elements.m;i++) segments_with_springs(i).x.Preallocate(max_connections);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    this->cfl_initialized=false;
}
//#####################################################################
// Function Set_Parameters
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Set_Parameters(const T youngs_modulus_input,const T overdamping_fraction_input,const T on_distance_input,const T off_distance_input, const int max_connections_input)
{
    this->youngs_modulus=youngs_modulus_input;
    this->overdamping_fraction=overdamping_fraction_input;
    this->on_distance=on_distance_input;
    this->off_distance=off_distance_input;
    this->max_connections=max_connections_input;
}
//#####################################################################
// Function Set_Parameters
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Set_Restlength(const T restlength_input)
{
    this->restlength=restlength_input;
}
//#####################################################################
// Function Update_Hierarchy
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Update_Hierarchy()
{
    // Initialize/Update Hierarchy
    if(!internal_curve.hierarchy || !external_curve.hierarchy){
        internal_curve.Initialize_Hierarchy(false);external_curve.Initialize_Hierarchy(false);}
    internal_curve.hierarchy->Update_Leaf_Boxes(particles.X);
    for(int s=1;s<=internal_curve.mesh.elements.m;s++) internal_curve.hierarchy->box_hierarchy(s).Change_Size(on_distance);
    internal_curve.hierarchy->Update_Nonleaf_Boxes();
    external_curve.hierarchy->Update_Leaf_Boxes(particles.X);
    for(int s=1;s<=external_curve.mesh.elements.m;s++) external_curve.hierarchy->box_hierarchy(s).Change_Size(on_distance);
    external_curve.hierarchy->Update_Nonleaf_Boxes();
}
//#####################################################################
// Function Update_Springs
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Update_Springs(const bool search_hierarchy)
{
    LOG::SCOPE scope("update adhesion springs","Update Adhesion Springs");
    // TODO: synchronize X's
    if(mpi_solids){
        LOG::Time("Gather particles");
        mpi_solids->All_Gather_Particles(particles.X,particles.V);
    }

    //PHYSBAM_DEBUG_PRINT("states",youngs_modulus,overdamping_fraction,on_distance,off_distance);
    // prune broken springs
    LOG::Time("Prune Particles");
    ARRAY<VECTOR<int,2> > deletion_list;
    //int first=true;
    //T old_off=off_distance;
    //VECTOR<int,2> min_elem;
    T_SPRING_HASH* new_springs=search_hierarchy?new T_SPRING_HASH():0; // many deletes can hurt open addressed hash table, use "generational" garbage collection
    if(new_springs) for(int i=1;i<=internal_segment_indices.m;i++){int segment=internal_segment_indices(i);segments_with_springs(segment).x.Remove_All();}
    for(typename HASHTABLE<VECTOR<int,2>,SPRING_STATE>::ITERATOR i(*springs);i.Valid();i.Next()){
        const VECTOR<int,2> &segment1_nodes=curve.mesh.elements(i.Key()[1]),&segment2_nodes=curve.mesh.elements(i.Key()[2]);
        SPRING_STATE& state=i.Data();
        SEGMENT_3D<T> segment1(particles.X.Subset(segment1_nodes)),segment2(particles.X.Subset(segment2_nodes));
        TV X1=(1-state.weights[1])*segment1.x1+state.weights[1]*segment1.x2;
        TV X2=(1-state.weights[2])*segment2.x1+state.weights[2]*segment2.x2;
        state.normal=X1-X2;
        state.distance=state.normal.Normalize();
        //if (off_distance>=state.distance){
        //    LOG::cout<<"ADH: Deleting and Setting"<<min_elem<<std::endl;
        //    if(first) first=false; else deletion_list.Append(min_elem);
        //    min_elem=i.Key();off_distance=state.distance;}
        //if(state.distance>off_distance || (state.marked && (segments_with_springs(i.Key()[1])>max_connections || segments_with_springs(i.Key()[2])>max_connections))){
        if(state.distance>off_distance) {
            deletion_list.Append(i.Key());
            intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]));
            intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment2_nodes[1],segment2_nodes[2],segment1_nodes[1],segment1_nodes[2]));}
        //else if(new_springs){
        //    if(!state.external || mpi_solids->partition_id_from_particle_index(state.nodes[1])==mpi_solids->Partition()) segments_with_springs(i.Key()[1]).x.Append(PAIR<T,int>(state.distance,i.Key()[2]));
        //    else segments_with_springs(i.Key()[2]).x.Append(PAIR<T,int>(state.distance,i.Key()[1]));
        else if(new_springs){
            segments_with_springs(i.Key()[1]).x.Append(PAIR<T,int>(state.distance,i.Key()[2]));
            new_springs->Insert(i.Key(),i.Data());}}
    //off_distance=old_off;
    if(!new_springs){ //not making new array, need to delete (as we haven't copied), also need to cull things that are no longer in range from external pairs (which will happen the same
                      //on ownership proc through above
        for(int i=1;i<=deletion_list.m;i++)springs->Delete(deletion_list(i));
        
        // note in the case of new_springs then we will delete the entire boundary list and make it again so we don't do this
        for(int i=external_springs.m;i>=1;i--){
            const VECTOR<int,2> &segment1_nodes=curve.mesh.elements(external_spring_segments(i)[1]),&segment2_nodes=curve.mesh.elements(external_spring_segments(i)[2]);
            SPRING_STATE& state=external_springs(i);
            SEGMENT_3D<T> segment1(particles.X.Subset(segment1_nodes)),segment2(particles.X.Subset(segment2_nodes)); 
            TV X1=(1-state.weights[1])*segment1.x1+state.weights[1]*segment1.x2;
            TV X2=(1-state.weights[2])*segment2.x1+state.weights[2]*segment2.x2;
            state.normal=X1-X2;
            state.distance=state.normal.Normalize();
            if(state.distance>off_distance){
                external_spring_segments.Remove_Index_Lazy(i);external_springs.Remove_Index_Lazy(i);
                intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]));
                intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment2_nodes[1],segment2_nodes[2],segment1_nodes[1],segment1_nodes[2]));}}}
    else{// swap springs
        for(int i=external_springs.m;i>=1;i--){
            const VECTOR<int,2> &segment1_nodes=curve.mesh.elements(external_spring_segments(i)[1]),&segment2_nodes=curve.mesh.elements(external_spring_segments(i)[2]);
            intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]));
            intersecting_edge_edge_pairs.Delete_If_Present(VECTOR<int,4>(segment2_nodes[1],segment2_nodes[2],segment1_nodes[1],segment1_nodes[2]));}
        delete springs;springs=new_springs;}

    //if(springs->Size()>0) return;

    if(search_hierarchy){ 
        LOG::Time("Update Hierarchy");
        Update_Hierarchy();
        
        LOG::Time("Search Hierarchy");
        // Find new pairs
        std::stringstream ss;ss<<"segments "<<internal_curve.mesh.elements.m<<std::endl;LOG::filecout(ss.str());
        for(int i=1;i<=internal_segment_indices.m;i++){int segment=internal_segment_indices(i);segments_with_springs(segment).y=false;} // we no longer have heaps
        EDGE_EDGE_ADHESION_VISITOR<TV> internal_visitor(*this,internal_segment_indices,internal_segment_indices,true);
        internal_curve.hierarchy->Intersection_List(*internal_curve.hierarchy,internal_visitor,ZERO());
        if(mpi_solids){
            EDGE_EDGE_ADHESION_VISITOR<TV> external_visitor(*this,internal_segment_indices,external_segment_indices,false);
            internal_curve.hierarchy->Intersection_List(*external_curve.hierarchy,external_visitor,ZERO());}

        if(mpi_solids){
            LOG::Time("Send and Receive Pairs");
            // Prepare to send pairs to other processors
            typedef PAIR<VECTOR<int,2>,VECTOR<T,2> > T_PAIR;
            ARRAY<ARRAY<T_PAIR>,PARTITION_ID> output_pairs(mpi_solids->Number_Of_Partitions()),input_pairs;
    
    
            // make a matrix of the data pairs
            ARRAY<ARRAY<int>,PARTITION_ID> my_boundary(mpi_solids->Number_Of_Partitions()); // for each other partition, the nodes of this processor that they need
            ARRAY<ARRAY<int>,PARTITION_ID> my_ghost(mpi_solids->Number_Of_Partitions()); // for each other partition, the nodes that this processor needs
    
            for(typename HASHTABLE<VECTOR<int,2>,SPRING_STATE>::ITERATOR i(*springs);i.Valid();i.Next()){
                if(i.Data().external){
                    const VECTOR<int,2>& segments=i.Key();
                    const VECTOR<int,4>& nodes=i.Data().nodes;
                    PARTITION_ID partition1=mpi_solids->partition_id_from_particle_index(nodes[1]),partition2=mpi_solids->partition_id_from_particle_index(nodes[3]);
                    if(partition1==mpi_solids->Partition()){
                        my_boundary(partition2).Append_Elements(nodes.template Slice<1,2>());
                        my_ghost(partition2).Append_Elements(nodes.template Slice<3,4>());
                        output_pairs(partition2).Append(T_PAIR(segments,i.Data().weights));}
                    else{
                        my_boundary(partition1).Append_Elements(nodes.template Slice<3,4>());
                        my_ghost(partition1).Append_Elements(nodes.template Slice<1,2>());
                        output_pairs(partition1).Append(T_PAIR(segments,i.Data().weights));}}}
            // exchange boundary pairs
            mpi_solids->All_Scatter_Adhesion_Pairs(output_pairs,input_pairs);
            // bin all the new springs and compute their parameters
            LOG::Time("Compute External Pair State");
            external_springs.Remove_All();external_spring_segments.Remove_All();
            const ARRAY_VIEW<T>& one_over_effective_mass=particles.one_over_effective_mass;
            for(PARTITION_ID partition(1);partition<=input_pairs.Size();partition++){
                ARRAY<T_PAIR>& pairs=input_pairs(partition);
                for(int j=1;j<=pairs.m;j++){
                    const VECTOR<int,2>& segment_indices=pairs(j).x;const VECTOR<T,2>& segment_weights=pairs(j).y;
                    SPRING_STATE state;
                    const VECTOR<int,2> &segment1_nodes=curve.mesh.elements(segment_indices[1]),&segment2_nodes=curve.mesh.elements(segment_indices[2]);
                    SEGMENT_3D<T> segment1(particles.X.Subset(segment1_nodes)),segment2(particles.X.Subset(segment2_nodes)); 
                    state.weights=segment_weights;
                    state.external=true;
                    state.nodes=segment1_nodes.Append_Elements(segment2_nodes);
                    TV X1=(1-state.weights[1])*segment1.x1+state.weights[1]*segment1.x2;
                    TV X2=(1-state.weights[2])*segment2.x1+state.weights[2]*segment2.x2;
                    state.normal=X1-X2;
                    state.distance=state.normal.Normalize();
                    T harmonic_mass=Pseudo_Inverse((T).25*(one_over_effective_mass(state.nodes[1])+one_over_effective_mass(state.nodes[2])
                            +one_over_effective_mass(state.nodes[3])+one_over_effective_mass(state.nodes[4])));
                    state.damping=overdamping_fraction*2*sqrt(youngs_modulus*restlength*harmonic_mass);
                    external_spring_segments.Append(segment_indices);
                    external_springs.Append(state);

                    // add them to collision ignore
                    intersecting_edge_edge_pairs.Set(VECTOR<int,4>(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]));
                    intersecting_edge_edge_pairs.Set(VECTOR<int,4>(segment2_nodes[1],segment2_nodes[2],segment1_nodes[1],segment1_nodes[2]));

                    // add the dependencies
                    PARTITION_ID partition1=mpi_solids->partition_id_from_particle_index(state.nodes[1]),partition2=mpi_solids->partition_id_from_particle_index(state.nodes[3]);
                    if(partition1==mpi_solids->Partition()){
                        my_boundary(partition2).Append_Elements(state.nodes.template Slice<1,2>());
                        my_ghost(partition2).Append_Elements(state.nodes.template Slice<3,4>());}
                    else{
                        my_boundary(partition1).Append_Elements(state.nodes.template Slice<3,4>());
                        my_ghost(partition1).Append_Elements(state.nodes.template Slice<1,2>());}
                }
            }
#if 0
            LOG::Time("update ghost/boundary particles");
            for(FRAGMENT_ID fragment(1);fragment<=mpi_solids->force_boundaries_of_fragment.Size();fragment++)
                mpi_solids->force_boundaries_of_fragment.Remove_All(); // TODO: NOT GENERAL ONLY FOR HAIR
            PHYSBAM_ASSERT(mpi_solids->fragments_of_partition(mpi_solids->Partition()).Size()==1);
            FRAGMENT_ID my_fragment=mpi_solids->fragments_of_partition(mpi_solids->Partition())(1);
            mpi_solids->force_boundaries_of_fragment.Resize(mpi_solids->partition_from_fragment.Size());
            ARRAY<PAIR<FRAGMENT_ID,ARRAY<int> > >& local_boundary_list=mpi_solids->force_boundaries_of_fragment(my_fragment);
            for(PARTITION_ID other_partition(1);other_partition<=my_boundary.Size();other_partition++) if(other_partition!=mpi_solids->Partition()){
                PHYSBAM_ASSERT(mpi_solids->fragments_of_partition(other_partition).Size()==1);
                FRAGMENT_ID other_fragment=mpi_solids->fragments_of_partition(other_partition)(1);
                local_boundary_list.Append(PAIR<FRAGMENT_ID,ARRAY<int> >(other_fragment,my_boundary(other_partition)));
                mpi_solids->force_boundaries_of_fragment(other_fragment).Append(PAIR<FRAGMENT_ID,ARRAY<int> >(my_fragment,my_ghost(other_partition)));
            }
            // TODO: This looks broken...
            ARRAY<PAIR<SUPER_FRAGMENT_ID,SUPER_FRAGMENT_ID> > swap_pairs;ARRAY<SUPER_FRAGMENT_ID> rebuild;rebuild.Append(SUPER_FRAGMENT_ID(1));
            LOG::Stop_Time();
#endif            
        }
    }

    // make cache of springs
    internal_springs.Remove_All();
    for(typename HASHTABLE<VECTOR<int,2>,SPRING_STATE>::ITERATOR i(*springs);i.Valid();i.Next()){
//        const VECTOR<int,2> &segment1_nodes=curve.mesh.elements(i.Key()[1]),&segment2_nodes=curve.mesh.elements(i.Key()[2]);
        SPRING_STATE& state=i.Data();
        internal_springs.Append(state);}    
}
//#####################################################################
// Function Update_Springs
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Update_Collisions_List()
{
    return;
}
//#####################################################################
// Function Update_Partitions
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Update_Partitions(bool restart,MPI_SOLIDS<TV>* mpi_solids,const std::string output_directory)
{
    this->mpi_solids=mpi_solids;
    if(mpi_solids) for(int i=1;i<=mesh.elements.m;i++){
        const VECTOR<int,2>& nodes=mesh.elements(i);
        VECTOR<PARTITION_ID,2> partitions(mpi_solids->partition_id_from_particle_index.Subset(nodes));
        PHYSBAM_ASSERT(partitions[1]==partitions[2]);
        if(partitions[1]==mpi_solids->Partition()){internal_segment_indices.Append(i);internal_mesh.elements.Append(nodes);}
        else{external_segment_indices.Append(i);external_mesh.elements.Append(nodes);}}
    else for(int i=1;i<=mesh.elements.m;i++){
        const VECTOR<int,2>& nodes=mesh.elements(i);
        internal_segment_indices.Append(i);internal_mesh.elements.Append(nodes);}

    // Cull hairs that are close at beginning
    if(restart){
        FILE_UTILITIES::Read_From_File<float>(output_directory+"/adhesion_existing",existing_pairs); // TODO: use real output_directory
        restart=false;}
    else{
        Update_Hierarchy();
        EDGE_EDGE_INITIAL_CULL_VISITOR<TV> internal_visitor(*this,internal_segment_indices,internal_segment_indices);
        EDGE_EDGE_INITIAL_CULL_VISITOR<TV> external_visitor(*this,internal_segment_indices,external_segment_indices);
        internal_curve.hierarchy->Intersection_List(*internal_curve.hierarchy,internal_visitor,ZERO());
        if(mpi_solids) internal_curve.hierarchy->Intersection_List(*external_curve.hierarchy,external_visitor,ZERO());
        FILE_UTILITIES::Write_To_File<float>(output_directory+"/adhesion_existing",existing_pairs);} // TODO: use real output_directory
}
//#####################################################################
// Function Update_Springs
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    Update_Springs(false);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    //for(HASHTABLE_ITERATOR<VECTOR<int,2>,const SPRING_STATE> i(*springs);i.Valid();i.Next()){
    //    const SPRING_STATE& state=i.Data();
    //    TV force=youngs_modulus*(state.distance/restlength-(T)1)*state.normal;
    //    F(state.nodes[1])-=((T)1-state.weights[1])*force;F(state.nodes[2])-=state.weights[1]*force;
    //    F(state.nodes[3])+=((T)1-state.weights[2])*force;F(state.nodes[4])+=state.weights[2]*force;}
    for(int i=1;i<=internal_springs.m;i++){
        const SPRING_STATE& state=internal_springs(i);
        TV force=youngs_modulus*(state.distance/restlength-(T)1)*state.normal;
        F(state.nodes[1])-=((T)1-state.weights[1])*force;F(state.nodes[2])-=state.weights[1]*force;
        F(state.nodes[3])+=((T)1-state.weights[2])*force;F(state.nodes[4])+=state.weights[2]*force;}
    for(int i=1;i<=external_springs.m;i++){
        const SPRING_STATE& state=external_springs(i);
        TV force=youngs_modulus*(state.distance/restlength-(T)1)*state.normal;
        F(state.nodes[1])-=((T)1-state.weights[1])*force;F(state.nodes[2])-=state.weights[1]*force;
        F(state.nodes[3])+=((T)1-state.weights[2])*force;F(state.nodes[4])+=state.weights[2]*force;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    //for(HASHTABLE_ITERATOR<VECTOR<int,2>,const SPRING_STATE> i(*springs);i.Valid();i.Next()){
    //    const SPRING_STATE& state=i.Data();
    //    TV force=state.damping/restlength*TV::Dot_Product(((T)1-state.weights[1])*V(state.nodes[1])+state.weights[1]*V(state.nodes[2])
    //                                   -((T)1-state.weights[2])*V(state.nodes[3])-state.weights[2]*V(state.nodes[4]),state.normal)*state.normal;
    //    F(state.nodes[1])-=((T)1-state.weights[1])*force;F(state.nodes[2])-=state.weights[1]*force;
    //    F(state.nodes[3])+=((T)1-state.weights[2])*force;F(state.nodes[4])+=state.weights[2]*force;
    //}
    for(int i=1;i<=internal_springs.m;i++){
        const SPRING_STATE& state=internal_springs(i);
        TV force=state.damping/restlength*TV::Dot_Product(((T)1-state.weights[1])*V(state.nodes[1])+state.weights[1]*V(state.nodes[2])
                                       -((T)1-state.weights[2])*V(state.nodes[3])-state.weights[2]*V(state.nodes[4]),state.normal)*state.normal;
        F(state.nodes[1])-=((T)1-state.weights[1])*force;F(state.nodes[2])-=state.weights[1]*force;
        F(state.nodes[3])+=((T)1-state.weights[2])*force;F(state.nodes[4])+=state.weights[2]*force;}
    for(int i=1;i<=external_springs.m;i++){
        const SPRING_STATE& state=external_springs(i);
        TV force=state.damping/restlength*TV::Dot_Product(((T)1-state.weights[1])*V(state.nodes[1])+state.weights[1]*V(state.nodes[2])
                                       -((T)1-state.weights[2])*V(state.nodes[3])-state.weights[2]*V(state.nodes[4]),state.normal)*state.normal;
        F(state.nodes[1])-=((T)1-state.weights[1])*force;F(state.nodes[2])-=state.weights[1]*force;
        F(state.nodes[3])+=((T)1-state.weights[2])*force;F(state.nodes[4])+=state.weights[2]*force;}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SEGMENT_ADHESION<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{   
    // TODO: Implement
}
//#####################################################################
// Function Write_To_File
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Write_State(STREAM_TYPE stream_type,const std::string& filename)
{   
    std::stringstream ss;ss<<"WRITING: Number of springs: "<<springs->Size()<<std::endl;LOG::filecout(ss.str());
    FILE_UTILITIES::Write_To_File(stream_type,filename,*springs);
}
//#####################################################################
// Function Write_To_File
//#####################################################################
template<class TV> void SEGMENT_ADHESION<TV>::
Read_State(STREAM_TYPE stream_type,const std::string& filename)
{
    FILE_UTILITIES::Read_From_File(stream_type,filename,*springs);
    std::stringstream ss;ss<<"READING: Number of springs: "<<springs->Size()<<std::endl;LOG::filecout(ss.str());
}
//#####################################################################
template class SEGMENT_ADHESION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEGMENT_ADHESION<VECTOR<double,3> >;
#endif
