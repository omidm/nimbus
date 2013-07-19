#ifndef _JOB
#define _JOB

#include <vector>
#include <set>
#include <map>
#include "data.h"

#define DataSetP std::set<Data *>
#define JobSetP std::set<Job *>

typedef std::vector<Data*> dataArray;
typedef void (*FuncToRun)(const dataArray&);


class Job
{
  public:
    Job(int, FuncToRun);
    FuncToRun handler;
    int id;

    DataSetP inputData;
    DataSetP outputData;

    JobSetP waitFor;
    JobSetP runBefore;

    void Execute(); 
    void Sleep();
    void Kill();

};



typedef std::map<int, Job*> JobMap;





#endif
