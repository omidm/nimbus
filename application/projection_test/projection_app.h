#include "shared/nimbus.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

class ProjectionApp : public Application {
    public:
    	ProjectionApp();
        virtual void Load();
};

class Main : public Job {
    public:
        Main(Application *app);
        virtual void Execute(std::string params, const DataArray& da);
        virtual Job* Clone();
};

namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;

template<class T_GRID>
class PCG_SPARSE_MPI:public NONCOPYABLE
{
public:
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    PCG_SPARSE<T>& pcg;
    SPARSE_MATRIX_PARTITION& partition;
    ARRAY<MPI::Datatype> boundary_datatypes,ghost_datatypes;
    ARRAY<ARRAY<int> > columns_to_send;
    ARRAY<ARRAY<int> > columns_to_receive;

    PCG_SPARSE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input)
        :pcg(pcg_input),comm(comm_input),thread_grid(0),mpi_threaded_grid(0),partition(partition_input)
    {}

    ~PCG_SPARSE_MPI()
    {MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);}

    template<class TYPE> TYPE Global_Sum(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::SUM,comm);return output;}

    template<class TYPE> TYPE Global_Max(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::MAX,comm);return output;}

    virtual void Fill_Ghost_Cells(VECTOR_ND<T>& v)
    {ARRAY<MPI::Request> requests;requests.Preallocate(2*partition.number_of_sides);
    for(int s=1;s<=partition.number_of_sides;s++)if(boundary_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append(comm.Isend(v.x-1,1,boundary_datatypes(s),partition.neighbor_ranks(s),s));
    for(int s=1;s<=partition.number_of_sides;s++)if(ghost_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append(comm.Irecv(v.x-1,1,ghost_datatypes(s),partition.neighbor_ranks(s),((s-1)^1)+1));
    MPI_UTILITIES::Wait_All(requests);}
    
};
}
