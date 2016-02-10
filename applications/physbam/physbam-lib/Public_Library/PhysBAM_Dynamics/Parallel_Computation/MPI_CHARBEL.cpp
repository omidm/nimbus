//#####################################################################
// Copyright 2008, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_CHARBEL
//#####################################################################
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_CHARBEL.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#endif
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{

#ifdef USE_MPI

//#####################################################################
// Constructor
//#####################################################################
template<class T> MPI_CHARBEL<T>::
MPI_CHARBEL()
    :local_tet_volume(TETRAHEDRALIZED_VOLUME<T>::Create()),current_tag(0)
{
    LOG::SCOPE scope("MPI INITIALIZE","Initializing MPI Charbel");
    number_of_processes=MPI::COMM_WORLD.Get_size();
    std::stringstream ss;
    ss<<"number of processes = "<<number_of_processes<<std::endl;

    comm=new MPI::Intracomm;
    *comm=MPI::COMM_WORLD.Dup();
    rank=comm->Get_rank();
    ss<<"global rank = "<<rank<<std::endl;
    LOG::filecout(ss.str());

    group=new MPI::Group(comm->Get_group());
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> MPI_CHARBEL<T>::
~MPI_CHARBEL()
{
    if(comm){comm->Free();delete comm;}
    if(group){group->Free();delete group;}
}
//#####################################################################
// Function Reduce_Max
//#####################################################################
template<class T> void MPI_CHARBEL<T>::
Reduce_Max(int& local_value) const
{
    int global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MAX,*comm);
    local_value=global_value;
}
//#####################################################################
// Function Setup_AEROF_PhysBAM_Mapping
//#####################################################################
template<class T> void MPI_CHARBEL<T>::
Setup_AEROF_PhysBAM_Mapping(TETRAHEDRALIZED_VOLUME<T>& tet_volume,ARRAY<ARRAY<int> >& tets_to_send,ARRAY<int>& local_to_global_map,
                const int& global_particle_count,const RANGE<TV>& domain)
{
    global_to_local_aerof_map.Resize(global_particle_count);ARRAYS_COMPUTATIONS::Fill(global_to_local_aerof_map,0);
    for(int i=1;i<=local_to_global_map.m;i++) global_to_local_aerof_map(local_to_global_map(i))=i;

    int tag=Get_Unique_Tag();
    ARRAY<ARRAY<char> > send_buffers(number_of_processes);ARRAY<MPI::Request> requests;
    for(int i=1;i<=number_of_processes;i++){
        INDIRECT_ARRAY<ARRAY<VECTOR<int,4> > > proc_tets=tet_volume.mesh.elements.Subset(tets_to_send(i));
        ARRAY<int> local_indices,global_indices;
        for(int t=1;t<=proc_tets.Size();t++){
            global_indices.Append_Elements(local_to_global_map.Subset(proc_tets(t)));
            local_indices.Append_Elements(proc_tets(t));}
        INDIRECT_ARRAY<ARRAY_VIEW<TV> > positions=tet_volume.particles.X.Subset(local_indices);
        int buffer_size=MPI_UTILITIES::Pack_Size(global_indices,positions,*comm)+1;
        send_buffers(i).Resize(buffer_size);int position=0;
        MPI_UTILITIES::Pack(global_indices,positions,send_buffers(i),position,*comm);
        requests.Append(comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,i-1,tag));}
    ARRAY<ARRAY<char> > recv_buffers(number_of_processes);
    for(int i=1;i<=number_of_processes;i++){
        MPI::Status status;
        comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(comm->Irecv(&(recv_buffers(i)(1)),recv_buffers(i).m,MPI_PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);

    interior_particles_to_recv.Resize(number_of_processes);ghost_particles_to_recv.Resize(number_of_processes);
    global_to_local_physbam_map.Resize(global_particle_count);ARRAYS_COMPUTATIONS::Fill(global_to_local_physbam_map,0);
    local_tet_volume->particles.array_collection->Delete_All_Elements();
    for(int i=1;i<=number_of_processes;i++){
        ARRAY<int> indices;int position=0;
        MPI_UTILITIES::Unpack(indices,recv_buffers(i),position,*comm);
        ARRAY<TV> positions(indices.m);INDIRECT_ARRAY<ARRAY<TV>,IDENTITY_ARRAY<int> > recv_positions(positions,IDENTITY_ARRAY<int>(indices.m));
        MPI_UTILITIES::Unpack(recv_positions,recv_buffers(i),position,*comm);
        for(int new_particle=1;new_particle<=indices.m;new_particle+=4){
            for(int p=0;p<4;p++){
                if(!global_to_local_physbam_map(indices(new_particle+p))){
                    global_to_local_physbam_map(indices(new_particle+p))=local_tet_volume->particles.array_collection->Add_Element();
                    local_tet_volume->particles.X(global_to_local_physbam_map(indices(new_particle+p)))=positions(new_particle+p);}
                if(domain.Inside(positions(new_particle+p),(T)0)) interior_particles_to_recv(i).Append_Unique(indices(new_particle+p));
                else ghost_particles_to_recv(i).Append_Unique(indices(new_particle+p));}
            local_tet_volume->mesh.elements.Append(VECTOR<int,4>(global_to_local_physbam_map(indices(new_particle)),global_to_local_physbam_map(indices(new_particle+1)),
                    global_to_local_physbam_map(indices(new_particle+2)),global_to_local_physbam_map(indices(new_particle+3))));}}
    local_tet_volume->Update_Number_Nodes();

    tag=Get_Unique_Tag();
    for(int i=1;i<=number_of_processes;i++){
        int buffer_size=MPI_UTILITIES::Pack_Size(interior_particles_to_recv(i),ghost_particles_to_recv(i),*comm)+1;
        send_buffers(i).Resize(buffer_size);int position=0;
        MPI_UTILITIES::Pack(interior_particles_to_recv(i),ghost_particles_to_recv(i),send_buffers(i),position,*comm);
        requests.Append(comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,i-1,tag));}
    for(int i=1;i<=number_of_processes;i++){
        MPI::Status status;
        comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(comm->Irecv(&(recv_buffers(i)(1)),recv_buffers(i).m,MPI_PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);
    interior_particles_to_send.Resize(number_of_processes);ghost_particles_to_send.Resize(number_of_processes);
    for(int i=1;i<=number_of_processes;i++){
        int position=0;
        MPI_UTILITIES::Unpack(interior_particles_to_send(i),ghost_particles_to_send(i),recv_buffers(i),position,*comm);}
    
    physbam_particles.array_collection->Add_Elements(local_tet_volume->particles.array_collection->Size());
    for(int i=1;i<=local_tet_volume->particles.array_collection->Size();i++) physbam_particles.X(i)=local_tet_volume->particles.X(i);
}
//#####################################################################
// Function Exchange_Compressible_Data
//#####################################################################
template<class T> void MPI_CHARBEL<T>::
Exchange_Compressible_Data(COMPRESSIBLE_FLUID_PARTICLES<TV>& particles_aerof)
{ 
    int interior_tag=Get_Unique_Tag();
    int ghost_tag=Get_Unique_Tag();
    ARRAY<ARRAY<char> > interior_send_buffers(number_of_processes),ghost_send_buffers(number_of_processes);ARRAY<MPI::Request> requests;
    for(int i=1;i<=number_of_processes;i++){
        if(i-1!=rank){
            INDIRECT_ARRAY<ARRAY<int> > proc_interior_particles_to_send=global_to_local_aerof_map.Subset(interior_particles_to_send(i)),
                proc_ghost_particles_to_send=global_to_local_aerof_map.Subset(ghost_particles_to_send(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> >& > interior_rho=particles_aerof.rho.Subset(proc_interior_particles_to_send),
                ghost_rho=particles_aerof.rho.Subset(proc_ghost_particles_to_send),interior_E=particles_aerof.E.Subset(proc_interior_particles_to_send),
                ghost_E=particles_aerof.E.Subset(proc_ghost_particles_to_send),
                interior_phi=particles_aerof.phi.Subset(proc_interior_particles_to_send),ghost_phi=particles_aerof.phi.Subset(proc_ghost_particles_to_send);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> >& > interior_u=particles_aerof.V.Subset(proc_interior_particles_to_send),
                ghost_u=particles_aerof.V.Subset(proc_ghost_particles_to_send),
                interior_grad_phi=particles_aerof.grad_phi.Subset(proc_interior_particles_to_send),ghost_grad_phi=particles_aerof.grad_phi.Subset(proc_ghost_particles_to_send);
            int buffer_size=MPI_UTILITIES::Pack_Size(interior_rho,interior_u,interior_E,interior_phi,interior_grad_phi,*comm)+1;
            interior_send_buffers(i).Resize(buffer_size);int position=0;
            MPI_UTILITIES::Pack(interior_rho,interior_u,interior_E,interior_phi,interior_grad_phi,interior_send_buffers(i),position,*comm);
            requests.Append(comm->Isend(&(interior_send_buffers(i)(1)),position,MPI::PACKED,i-1,interior_tag));
            buffer_size=MPI_UTILITIES::Pack_Size(ghost_rho,ghost_u,ghost_E,ghost_phi,ghost_grad_phi,*comm)+1;
            ghost_send_buffers(i).Resize(buffer_size);position=0;
            MPI_UTILITIES::Pack(ghost_rho,ghost_u,ghost_E,ghost_phi,ghost_grad_phi,ghost_send_buffers(i),position,*comm);
            requests.Append(comm->Isend(&(ghost_send_buffers(i)(1)),position,MPI::PACKED,i-1,ghost_tag));}}
    ARRAY<ARRAY<char> > interior_recv_buffers(number_of_processes),ghost_recv_buffers(number_of_processes);
    for(int i=1;i<=number_of_processes;i++){
        if(i-1!=rank){
            MPI::Status status;
            comm->Probe(i-1,interior_tag,status);
            interior_recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
            requests.Append(comm->Irecv(&(interior_recv_buffers(i)(1)),interior_recv_buffers(i).m,MPI_PACKED,i-1,interior_tag));
            comm->Probe(i-1,ghost_tag,status);
            ghost_recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
            requests.Append(comm->Irecv(&(ghost_recv_buffers(i)(1)),ghost_recv_buffers(i).m,MPI_PACKED,i-1,ghost_tag));}}
    MPI_UTILITIES::Wait_All(requests);
    for(int i=1;i<=number_of_processes;i++){
        if(i-1!=rank){
            int position=0;
            INDIRECT_ARRAY<ARRAY<int> > local_indices_of_interior_particles_to_recv=global_to_local_physbam_map.Subset(interior_particles_to_recv(i)),
                local_indices_of_ghost_particles_to_recv=global_to_local_physbam_map.Subset(ghost_particles_to_recv(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> > > interior_rho(physbam_particles.rho,local_indices_of_interior_particles_to_recv),
                ghost_rho(physbam_particles.rho,local_indices_of_ghost_particles_to_recv),
                interior_E(physbam_particles.E,local_indices_of_interior_particles_to_recv),
                ghost_E(physbam_particles.E,local_indices_of_ghost_particles_to_recv),
                interior_phi(physbam_particles.phi,local_indices_of_interior_particles_to_recv),
                ghost_phi(physbam_particles.phi,local_indices_of_ghost_particles_to_recv);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> > > interior_u(physbam_particles.V,local_indices_of_interior_particles_to_recv),
                ghost_u(physbam_particles.V,local_indices_of_ghost_particles_to_recv),
                interior_grad_phi(physbam_particles.grad_phi,local_indices_of_interior_particles_to_recv),
                ghost_grad_phi(physbam_particles.grad_phi,local_indices_of_ghost_particles_to_recv);
            MPI_UTILITIES::Unpack(interior_rho,interior_u,interior_E,interior_phi,interior_grad_phi,interior_recv_buffers(i),position,*comm);
            position=0;MPI_UTILITIES::Unpack(ghost_rho,ghost_u,ghost_E,ghost_phi,ghost_grad_phi,ghost_recv_buffers(i),position,*comm);}
        else{
            INDIRECT_ARRAY<ARRAY<int> > proc_interior_particles_to_send=global_to_local_aerof_map.Subset(interior_particles_to_send(i)),
                proc_ghost_particles_to_send=global_to_local_aerof_map.Subset(ghost_particles_to_send(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> >& > interior_rho_s=particles_aerof.rho.Subset(proc_interior_particles_to_send),
                ghost_rho_s=particles_aerof.rho.Subset(proc_ghost_particles_to_send),interior_E_s=particles_aerof.E.Subset(proc_interior_particles_to_send),
                ghost_E_s=particles_aerof.E.Subset(proc_ghost_particles_to_send),interior_phi_s=particles_aerof.phi.Subset(proc_interior_particles_to_send),
                ghost_phi_s=particles_aerof.phi.Subset(proc_ghost_particles_to_send);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> >& > interior_u_s=particles_aerof.V.Subset(proc_interior_particles_to_send),
                ghost_u_s=particles_aerof.V.Subset(proc_ghost_particles_to_send),interior_grad_phi_s=particles_aerof.grad_phi.Subset(proc_interior_particles_to_send),
                ghost_grad_phi_s=particles_aerof.grad_phi.Subset(proc_ghost_particles_to_send);
            INDIRECT_ARRAY<ARRAY<int> > local_indices_of_interior_particles_to_recv=global_to_local_physbam_map.Subset(interior_particles_to_recv(i)),
                local_indices_of_ghost_particles_to_recv=global_to_local_physbam_map.Subset(ghost_particles_to_recv(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> > > interior_rho_r(physbam_particles.rho,local_indices_of_interior_particles_to_recv),
                ghost_rho_r(physbam_particles.rho,local_indices_of_ghost_particles_to_recv),
                interior_E_r(physbam_particles.E,local_indices_of_interior_particles_to_recv),
                ghost_E_r(physbam_particles.E,local_indices_of_ghost_particles_to_recv),
                interior_phi_r(physbam_particles.phi,local_indices_of_interior_particles_to_recv),
                ghost_phi_r(physbam_particles.phi,local_indices_of_ghost_particles_to_recv);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> > > interior_u_r(physbam_particles.V,local_indices_of_interior_particles_to_recv),
                ghost_u_r(physbam_particles.V,local_indices_of_ghost_particles_to_recv),interior_grad_phi_r(physbam_particles.grad_phi,local_indices_of_interior_particles_to_recv),
                ghost_grad_phi_r(physbam_particles.grad_phi,local_indices_of_ghost_particles_to_recv);
            interior_rho_r=interior_rho_s;ghost_rho_r=ghost_rho_s;interior_E_r=interior_E_s;ghost_E_r=ghost_E_s;interior_phi_r=interior_phi_s;ghost_phi_r=ghost_phi_s;
            interior_u_r=interior_u_s;ghost_u_r=ghost_u_s;interior_grad_phi_r=interior_grad_phi_s;ghost_grad_phi_r=ghost_grad_phi_s;}}
}
//#####################################################################
// Function Exchange_Back_Compressible_Data
//#####################################################################
template<class T> void MPI_CHARBEL<T>::
Exchange_Back_Compressible_Data(COMPRESSIBLE_FLUID_PARTICLES<TV>& particles_aerof)
{ 
    int tag=Get_Unique_Tag();
    ARRAY<ARRAY<char> > send_buffers(number_of_processes);ARRAY<MPI::Request> requests;
    for(int i=1;i<=number_of_processes;i++){
        if(i-1!=rank){
            INDIRECT_ARRAY<ARRAY<int> > local_indices_of_interior_particles_to_recv=global_to_local_physbam_map.Subset(interior_particles_to_recv(i)),
                local_indices_of_ghost_particles_to_recv=global_to_local_physbam_map.Subset(ghost_particles_to_recv(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> > > interior_rho(physbam_particles.rho,local_indices_of_interior_particles_to_recv),
                interior_E(physbam_particles.E,local_indices_of_interior_particles_to_recv),interior_phi(physbam_particles.phi,local_indices_of_interior_particles_to_recv);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> > > interior_u(physbam_particles.V,local_indices_of_interior_particles_to_recv),interior_grad_phi(physbam_particles.grad_phi,local_indices_of_interior_particles_to_recv);
            int buffer_size=MPI_UTILITIES::Pack_Size(interior_rho,interior_u,interior_E,interior_phi,interior_grad_phi,*comm)+1;
            send_buffers(i).Resize(buffer_size);int position=0;
            MPI_UTILITIES::Pack(interior_rho,interior_u,interior_E,interior_phi,interior_grad_phi,send_buffers(i),position,*comm);
            requests.Append(comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,i-1,tag));}}
    ARRAY<ARRAY<char> > recv_buffers(number_of_processes);
    for(int i=1;i<=number_of_processes;i++){
        if(i-1!=rank){
            MPI::Status status;
            comm->Probe(i-1,tag,status);
            recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
            requests.Append(comm->Irecv(&(recv_buffers(i)(1)),recv_buffers(i).m,MPI_PACKED,i-1,tag));}}
    MPI_UTILITIES::Wait_All(requests);
    for(int i=1;i<=number_of_processes;i++){
        if(i-1!=rank){
            int position=0;
            INDIRECT_ARRAY<ARRAY<int> > proc_interior_particles_to_send=global_to_local_aerof_map.Subset(interior_particles_to_send(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> >& > interior_rho=particles_aerof.rho.Subset(proc_interior_particles_to_send),
                interior_E=particles_aerof.E.Subset(proc_interior_particles_to_send),interior_phi=particles_aerof.phi.Subset(proc_interior_particles_to_send);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> >& > interior_u=particles_aerof.V.Subset(proc_interior_particles_to_send),
                interior_grad_phi=particles_aerof.grad_phi.Subset(proc_interior_particles_to_send);
            MPI_UTILITIES::Unpack(interior_rho,interior_u,interior_E,interior_phi,interior_grad_phi,recv_buffers(i),position,*comm);}
        else{
            INDIRECT_ARRAY<ARRAY<int> > local_indices_of_interior_particles_to_recv=global_to_local_physbam_map.Subset(interior_particles_to_recv(i)),
                local_indices_of_ghost_particles_to_recv=global_to_local_physbam_map.Subset(ghost_particles_to_recv(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> > > interior_rho_r(physbam_particles.rho,local_indices_of_interior_particles_to_recv),
                interior_E_r(physbam_particles.E,local_indices_of_interior_particles_to_recv),interior_phi_r(physbam_particles.phi,local_indices_of_interior_particles_to_recv);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> > > interior_u_r(physbam_particles.V,local_indices_of_interior_particles_to_recv),
                interior_grad_phi_r(physbam_particles.grad_phi,local_indices_of_interior_particles_to_recv);
            INDIRECT_ARRAY<ARRAY<int> > proc_interior_particles_to_send=global_to_local_aerof_map.Subset(interior_particles_to_send(i));
            INDIRECT_ARRAY<ARRAY_VIEW<T>,INDIRECT_ARRAY<ARRAY<int> >& > interior_rho_s=particles_aerof.rho.Subset(proc_interior_particles_to_send),
                interior_E_s=particles_aerof.E.Subset(proc_interior_particles_to_send),interior_phi_s=particles_aerof.phi.Subset(proc_interior_particles_to_send);
            INDIRECT_ARRAY<ARRAY_VIEW<TV>,INDIRECT_ARRAY<ARRAY<int> >& > interior_u_s=particles_aerof.V.Subset(proc_interior_particles_to_send),
                interior_grad_phi_s=particles_aerof.grad_phi.Subset(proc_interior_particles_to_send);
            interior_rho_s=interior_rho_r;interior_E_s=interior_E_r;interior_phi_r=interior_phi_s;interior_u_s=interior_u_r;interior_grad_phi_r=interior_grad_phi_s;}}
}
//#####################################################################
#else
//#####################################################################
template<class T> MPI_CHARBEL<T>::MPI_CHARBEL(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T> MPI_CHARBEL<T>::~MPI_CHARBEL(){}
template<class T> void MPI_CHARBEL<T>::Reduce_Max(int& local_value) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T> void MPI_CHARBEL<T>::Setup_AEROF_PhysBAM_Mapping(TETRAHEDRALIZED_VOLUME<T>& tet_volume,ARRAY<ARRAY<int> >& tets_to_send,
                                   ARRAY<int>& local_to_global_map,const int& global_particle_count,
                                   const RANGE<TV>& domain){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T> void MPI_CHARBEL<T>::Exchange_Compressible_Data(COMPRESSIBLE_FLUID_PARTICLES<TV>& particles_aerof){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T> void MPI_CHARBEL<T>::Exchange_Back_Compressible_Data(COMPRESSIBLE_FLUID_PARTICLES<TV>& particles_aerof){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
#endif
template class MPI_CHARBEL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MPI_CHARBEL<double>;
#endif
}
