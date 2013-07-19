#ifndef _DATA
#define _DATA

#include <set>
#include <map>
#include "cluster.h"

#define Hosts std::set<Computer*>
#define Neighbors std::set<Data*>

//typedef std::vector<void *> argArray;
//typedef void (*FuncToCreate)(const argArray&);


class Data
{
  public:
    Data();

    Hosts hosts;
    //Neighbors
    int id; 


    void Create();
    void Destroy(Computer);
    void Copy(Computer, Computer);
    void Migrate(Computer, Computer);

    /* 
     * If you would like to help the scheduler,
     * also provide followinf functions
     */
    bool advanceData;
    void (*funcSplit)(Data *, Data *);
    void (*funcMerge)(Data *, Data *);


};



typedef std::map<int, Data*> DataMap;




#endif
