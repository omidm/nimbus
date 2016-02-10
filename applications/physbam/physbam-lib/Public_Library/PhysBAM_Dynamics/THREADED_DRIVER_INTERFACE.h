//#####################################################################
// Copyright 2011, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __THREADED_DRIVER_INTERFACE__
#define __THREADED_DRIVER_INTERFACE__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parallel_Computation/INT_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class T_DRIVER,class TV>
class THREADED_DRIVER_INTERFACE
{
    typedef typename TV::SCALAR T;
public:
    ARRAY<THREAD_PACKAGE> buffers;
    ARRAY<T_DRIVER*> drivers;
    THREAD_QUEUE thread_queue;

    THREADED_DRIVER_INTERFACE(int number_of_threads)
        :drivers(number_of_threads),thread_queue(number_of_threads)
    {
        //For now the caller needs to setup the driver list as we cannot know the example type
    }

    virtual ~THREADED_DRIVER_INTERFACE()
    {
        for(int i=2;i<=thread_queue.Number_Of_Threads();i++) delete drivers(i);
    }

    void Execute_Main_Program()
    {
        INT_ITERATOR_THREADED_ALPHA<THREADED_DRIVER_INTERFACE<T_DRIVER,TV> > thread_iterator(&thread_queue);
        thread_iterator.Run(*this,&THREADED_DRIVER_INTERFACE<T_DRIVER,TV>::Execute_Main_Program_Threaded);
    }
    void Execute_Main_Program_Threaded(int tid)
    {
        drivers(tid)->example.fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(buffers,tid,thread_queue.Number_Of_Threads(),*drivers(tid)->example.fluids_parameters.grid,3);    
        drivers(tid)->Execute_Main_Program();
    }

    void Initialize()
    {
        INT_ITERATOR_THREADED_ALPHA<THREADED_DRIVER_INTERFACE<T_DRIVER,TV> > thread_iterator(&thread_queue);
        thread_iterator.Run(*this,&THREADED_DRIVER_INTERFACE<T_DRIVER,TV>::Initialize_Threaded);
    }
    void Initialize_Threaded(int tid)
    {
        drivers(tid)->example.fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(buffers,tid,thread_queue.Number_Of_Threads(),*drivers(tid)->example.fluids_parameters.grid,3);    
        drivers(tid)->Initialize();
    }

    void Advance_To_Target_Time(const T target_time)
    {
        INT_ITERATOR_THREADED_ALPHA<THREADED_DRIVER_INTERFACE<T_DRIVER,TV> > thread_iterator(&thread_queue);
        thread_iterator.template Run<T>(*this,&THREADED_DRIVER_INTERFACE<T_DRIVER,TV>::Advance_To_Target_Time_Threaded,target_time);
    }
    void Advance_To_Target_Time_Threaded(const T target_time,int tid)
    {
        drivers(tid)->Advance_To_Target_Time(target_time);
    }

    void Simulate_To_Frame(const int frame)
    {
        INT_ITERATOR_THREADED_ALPHA<THREADED_DRIVER_INTERFACE<T_DRIVER,TV> > thread_iterator(&thread_queue);
        thread_iterator.template Run<int>(*this,&THREADED_DRIVER_INTERFACE<T_DRIVER,TV>::Simulate_To_Frame_Threaded,frame);
    }
    void Simulate_To_Frame_Threaded(const int frame,int tid)
    {
        drivers(tid)->Simulate_To_Frame(frame);
    }
 
    void Write_Output_Files(const int frame)
    {
        INT_ITERATOR_THREADED_ALPHA<THREADED_DRIVER_INTERFACE<T_DRIVER,TV> > thread_iterator(&thread_queue);
        thread_iterator.template Run<int>(*this,&THREADED_DRIVER_INTERFACE<T_DRIVER,TV>::Write_Output_Files_Threaded,frame);
    }
    void Write_Output_Files_Threaded(const int frame,int tid)
    {
        drivers(tid)->Write_Output_Files(frame);
    }

//#####################################################################
};
}
#endif
