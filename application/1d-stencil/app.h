#ifndef _1DSTENCIL
#define _1DSTENCIL

#include "application.h"
#include "job.h"
#include "data.h"

void init (const dataArray&);
void print (const dataArray&);
void stenLeft (const dataArray&);
void stenRight (const dataArray&);
void updateghostLeft (const dataArray&); 
void updateghostRight (const dataArray&); 

void applySten (int *, int *, int);


class Vec : public Data 
{
  public:
    Vec(int);
    int size;
    int *arr;
    void Create();
};

class App : public Application
{
  public:
    App();
    virtual void loadApp();

};











#endif
