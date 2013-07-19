#include <iostream>
#include "app.h"

#define ML 4
#define GL 1

using namespace std;

Vec::Vec(int size)
{
  this->size = size;
};

void Vec::Create()
{
  arr = new int[size];
};

App::App()
{

};

void App::loadApp()
{

  cout << "Start Creating Data and Job Objects" << endl;
  
  Vec * mainLeft = new Vec (ML);
  Vec * mainRight = new Vec (ML);
  Vec * ghostLeft = new Vec (GL);
  Vec * ghostRight = new Vec (GL);

  Job * initialize = new Job (0, &init);
  Job * printing = new Job (1, &print);
  Job * applyLeft = new Job (2, &stenLeft);
  Job * applyRight = new Job (3, &stenRight);
  Job * updateLeft = new Job (4, &updateghostLeft);
  Job * updateRight = new Job (5, &updateghostRight);

  cout << "Finished Creating Data and Job Objects" << endl;

  /*
   * Now based on this do the followings
   *
   * 1) Populate the job dependency graph.
   *    (i.e. JobSets and DataSets for each Job)
   *
   * 2) Populate the dataMap and jobMap of the application.
   *    (i.e. key value pairs to map the data and job from scheduler commands.)
   *
   */

  initialize->runBefore.insert(applyLeft);
  initialize->runBefore.insert(applyRight);

  applyLeft->waitFor.insert(initialize);
  applyLeft->runBefore.insert(updateLeft);
  applyLeft->runBefore.insert(updateRight);

  // Continue like this



  dataMap[0] = mainLeft;
  dataMap[1] = mainRight;
  dataMap[2] = ghostLeft;
  dataMap[3] = ghostRight;

  jobMap[0] = initialize;
  jobMap[1] = printing;
  jobMap[2] = applyLeft;
  jobMap[3] = applyRight;
  jobMap[4] = updateLeft;
  jobMap[5] = updateRight;


};



void init (const dataArray& da)
{
  Vec *d = (Vec*)(da[0]);
  for(int i = 0; i < d->size ; i++)
    d->arr[i] = 0;
};


void print (const dataArray& da)
{
  Vec *d = (Vec*)(da[0]);
  for(int i = 0; i < d->size; i++)
    cout << d->arr[i] << ", ";
};


void stenLeft (const dataArray& da)
{
  int sten [] = {-1, +2, -1, 0};

  Vec *d1 = (Vec*)(da[0]);
  Vec *d2 = (Vec*)(da[1]);

  int len = d1->size;
  int *main = d1->arr;
  int *ghost = d2->arr;
  int temp[len];

  temp[0] = sten[0] * sten[3] + sten[1] * main[0] + sten[2] * main[1]; 
  temp[len-1] = sten[0] * main[len-2] + sten[1] * main[len-1] + sten[2] * ghost[0]; 
  for (int i = 1; i < (len - 1); i++)
    temp[i] = sten[0] * main[i-1] + sten[1] * main[i] + sten[2] * main[i+1];
  
  for(int i = 0; i < len; i++)
    main[i] = temp[i];

};

void stenRight (const dataArray& da)
{

};

void updateghostLeft (const dataArray& da)
{
  Vec *d1 = (Vec*)(da[0]);
  Vec *d2 = (Vec*)(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[0];
}; 

void updateghostRight (const dataArray& da)
{

};





