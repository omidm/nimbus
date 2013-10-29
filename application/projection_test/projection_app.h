#ifndef NIMBUS_APPLICATION_PROJECTION_TEST_PROJECTION_APP_H_
#define NIMBUS_APPLICATION_PROJECTION_TEST_PROJECTION_APP_H_
#include "shared/nimbus.h"
#include "shared/parser.h"
#include "shared/nimbus_types.h"
#include "physbam_include.h"
#include "worker/application.h"
#include "worker/job.h"
#include "worker/data.h"
#include "protocol_buffer/vector_msg.pb.h"
#define LEN 4
#define DESIRED_ITERATIONS 100
#define GLOBAL_TOLERANCE 1e-3
#define NUM_OF_FORLOOP_INPUTS 9

using nimbus::Job;
using nimbus::Data;
using nimbus::Application;
using namespace PhysBAM;
typedef float T;

class ProjectionApp : public Application {
    public:
    	ProjectionApp();
        virtual void Load();
};

class Main : public Job {
    public:
        Main(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Init : public Job {
    public:
    	Init(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Condition : public Job {
    public:
    	Project_Forloop_Condition(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part1 : public Job {
    public:
    	Project_Forloop_Part1(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part2 : public Job {
    public:
    	Project_Forloop_Part2(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part3 : public Job {
    public:
    	Project_Forloop_Part3(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Project_Forloop_Part4 : public Job {
    public:
    	Project_Forloop_Part4(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Global_Sum : public Job {
    public:
		Global_Sum(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Global_Max : public Job {
    public:
		Global_Max(Application *app);
        virtual void Execute(Parameter params, const DataArray& da);
        virtual Job* Clone();
};

class Vec : public Data {
  public:
    explicit Vec(int size);
    virtual ~Vec();

    virtual void Create();
    virtual void Destroy();
    virtual Data * Clone();
    virtual void Copy(Data* from);
    virtual bool Serialize(SerializedData* ser_data);
    virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

    int size();
    T* arr();

  private:
    int size_;
    T *arr_;
};

#endif

/*
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


Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance,const bool recompute_preconditioner)
{
    Initialize_Datatypes();
    int local_n=A.n,interior_n=partition.interior_indices.Size()+1;
    int global_n=Global_Sum(interior_n);
    T global_tolerance=Global_Max(tolerance);
    int desired_iterations=global_n;
    if(pcg.maximum_iterations) desired_iterations=min(desired_iterations,pcg.maximum_iterations);

    VECTOR_ND<T> temp(local_n,false),p(local_n,false),z_interior(interior_n,false);

    // build interior views of x,b,p,z,temp
    VECTOR_ND<T> x_interior,b_interior,p_interior,temp_interior;
    x_interior.Set_Subvector_View(x,partition.interior_indices);
    b_interior.Set_Subvector_View(b,partition.interior_indices);
    p_interior.Set_Subvector_View(p,partition.interior_indices);
    temp_interior.Set_Subvector_View(temp,partition.interior_indices);

    // find initial residual, r=b-Ax - reusing b for the residual
    Fill_Ghost_Cells(x);
    A.Times(x,temp);b_interior-=temp_interior;    
    if(Global_Max(b_interior.Max_Abs())<=global_tolerance){return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    if(pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)){
        delete A.C;A.C=A.Create_Submatrix(partition.interior_indices);
        A.C->In_Place_Incomplete_Cholesky_Factorization(pcg.modified_incomplete_cholesky,pcg.modified_incomplete_cholesky_coefficient,
            pcg.preconditioner_zero_tolerance,pcg.preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        if(pcg.incomplete_cholesky){
            // solve Mz=r
            A.C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
            A.C->Solve_Backward_Substitution(temp_interior,z_interior,false,true);} // diagonal is inverted to save on divides
        else z_interior=b_interior; // set z=r when there is no preconditioner

        // update search direction
        rho_old=rho;rho=Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior,b_interior));
        T beta=0;if(iteration==1) p_interior=z_interior;else{beta=(T)(rho/rho_old);for(int i=1;i<=interior_n;i++) p_interior(i)=z_interior(i)+beta*p_interior(i);} // when iteration=1, beta=0

        // update solution and residual
        Fill_Ghost_Cells(p);
        A.Times(p,temp);
        T alpha=(T)(rho/Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior)));
        for(int i=1;i<=interior_n;i++){x_interior(i)+=alpha*p_interior(i);b_interior(i)-=alpha*temp_interior(i);}

        T residual=Global_Max(b_interior.Max_Abs());

        if(residual<=global_tolerance){break;}
        if(iteration==desired_iterations){break;}
    }
  
    Fill_Ghost_Cells(x);
}
*/
