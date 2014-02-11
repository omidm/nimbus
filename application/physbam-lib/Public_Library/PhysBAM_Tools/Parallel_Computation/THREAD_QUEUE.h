#ifndef __THREAD_QUEUE__
#define __THREAD_QUEUE__
#include <cstdio>
#include <cstdlib>
namespace PhysBAM {

class THREAD_QUEUE {
 public:
  struct TASK {
    virtual ~TASK() {}
    virtual void Run(const int tid) = 0;
  };
  THREAD_QUEUE(const int thread_count, const bool set_affinity=false) {
    if (thread_count != -1) {
      printf("ERROR: class THREAD_QUEUE should never be initialized!\n");
      exit(-1);
    }
  }
  virtual ~THREAD_QUEUE() {}

  virtual void Queue(TASK* task) {}
  virtual void Wait() {}
  virtual int Number_Of_Threads() {
    return -1;
  }
};

}  // namespace PhysBAM
#endif

/*
#ifdef USE_PTHREADS
#include <pthread.h>
#endif
class THREAD_QUEUE
{
#ifdef USE_PTHREADS
    ARRAY<pthread_t> threads;
    pthread_attr_t attr;
    pthread_mutex_t queue_lock;
    pthread_cond_t done_condition,todo_condition;
    int active_threads,inactive_threads;

    ARRAY<PAIR<THREAD_QUEUE*,int> > threadStates;
#endif

public:
    struct TASK
    {
        virtual ~TASK(){};
        virtual void Run(const int tid)=0;
    };

#ifdef USE_PTHREADS
    struct EXITER:public TASK
    {
        void Run(const int threadid)
        {pthread_exit(0);}
    };

private:
    QUEUE<TASK*> queue;
public:
#endif

    THREAD_QUEUE(const int thread_count,const bool set_affinity=false);
    ~THREAD_QUEUE();

    void Queue(TASK* task);
    void Wait();
    int Number_Of_Threads();

#ifdef USE_PTHREADS
    static void* Thread_Routine(void* data)
    {
        PAIR<THREAD_QUEUE*,int>* threadState=(PAIR<THREAD_QUEUE*,int>*)data;
        THREAD_QUEUE& queue=*threadState->x;
        while(1){
            pthread_mutex_lock(&queue.queue_lock);
            while(queue.queue.Empty()){
                queue.active_threads--;
                if(queue.active_threads==0) pthread_cond_signal(&queue.done_condition);
                queue.inactive_threads++;
                pthread_cond_wait(&queue.todo_condition,&queue.queue_lock);
                queue.active_threads++;queue.inactive_threads--;}
            TASK* work=queue.queue.Dequeue();
            pthread_mutex_unlock(&queue.queue_lock);
            work->Run(threadState->y);
            delete work;
        }
        return 0;
    }
#endif
//#####################################################################
};
}
*/
