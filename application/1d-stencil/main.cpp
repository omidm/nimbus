#include <iostream>
#include "scheduler.h"
#include "data.h"

using namespace std;

#define ML 4
#define GL 1
#define ITNUM 5

void createArray (int **, int);
void initArray (int *, int, int);
void printArray (int *, int);
void applyStenLeft (int *, int *, int *, int);
void applyStenRight (int *, int *, int *, int);
void updateGhostLeft (int *, int *, int);
void updateGhostRight (int *, int *, int);

void applySten (int *, int *, int);

int main (int argc, char *argv[]) {

  std::cout << "Nimbus is up!" << std::endl;
  Scheduler s;
  s.Run();
  
  int sten [] = {-1, +2, -1, 0};

  int * mainLeft;
  int * mainRight;
  int * ghostLeft;
  int * ghostRight;

  createArray(&mainLeft, ML);
  createArray(&mainRight, ML);
  createArray(&ghostLeft, GL);
  createArray(&ghostRight, GL);
  
  initArray(mainLeft, ML, 5);
  initArray(mainRight, ML, 3);
  
  initArray(ghostLeft, GL, 3);
  initArray(ghostRight, GL, 5);


  cout << "*** Before Stencil ***" << endl;
  printArray(mainLeft, ML);
  printArray(mainRight, ML);
  cout << endl;

  for(int i = 0; i < ITNUM; i ++)
  {
    updateGhostLeft(mainRight, ghostLeft, ML);
    updateGhostRight(mainLeft, ghostRight, ML);

    applyStenLeft(mainLeft, ghostLeft, sten, ML);
    applyStenRight(mainRight, ghostRight, sten, ML);

    cout << "***After Stencil num " << i+1 << endl;
    printArray(mainLeft, ML);
    printArray(mainRight, ML);
    cout << endl;

  }

  cout << "******************************************************" << endl;
  cout << "******************************************************" << endl;


  int total [] = {5, 5, 5, 5, 3, 3, 3, 3};
  
  cout << "*** Before Stencil ***" << endl;
  printArray(total, 2 * ML);
  cout << endl;

  for(int i = 0; i < ITNUM; i ++)
  {
    applySten(total, sten, 2 * ML);
    cout << "***After Stencil num " << i+1 << endl;
    printArray(total, 2 * ML);
    cout << endl;
  }




}


void createArray (int ** arrayP, int size)
{
  *arrayP = new int[size]; 
}

void initArray (int * arrayP, int len, int val)
{
  for(int i = 0; i < len; i++)
    arrayP[i] = val;
};

void printArray(int * arrayP, int len)
{
  for(int i = 0; i < len; i++)
    cout << arrayP[i] << ", ";
  //cout << "\b\b " << endl;
};


void applyStenLeft (int * main, int * ghost, int * sten, int len)
{
  int temp[len];

  temp[0] = sten[0] * sten[3] + sten[1] * main[0] + sten[2] * main[1]; 
  temp[len-1] = sten[0] * main[len-2] + sten[1] * main[len-1] + sten[2] * ghost[0]; 
  for (int i = 1; i < (len - 1); i++)
    temp[i] = sten[0] * main[i-1] + sten[1] * main[i] + sten[2] * main[i+1];
  
  for(int i = 0; i < len; i++)
    main[i] = temp[i];

};


void applyStenRight (int * main, int * ghost, int * sten, int len)
{
  int temp[len];

  temp[0] = sten[0] * ghost[0] + sten[1] * main[0] + sten[2] * main[1]; 
  temp[len-1] = sten[0] * main[len-2] + sten[1] * main[len-1] + sten[2] * sten[3]; 
  for (int i = 1; i < (len - 1); i++)
    temp[i] = sten[0] * main[i-1] + sten[1] * main[i] + sten[2] * main[i+1];
  
  for(int i = 0; i < len; i++)
    main[i] = temp[i];

};


void updateGhostLeft (int * main, int * ghost, int len)
{
  ghost[0] = main[0];
};


void updateGhostRight (int * main, int * ghost, int len)
{
  ghost[0] = main[len-1];
};



void applySten (int * main, int * sten, int len)
{
  int temp[len];

  temp[0] = sten[0] * sten[3] + sten[1] * main[0] + sten[2] * main[1]; 
  temp[len-1] = sten[0] * main[len-2] + sten[1] * main[len-1] + sten[2] * sten[3]; 
  for (int i = 1; i < (len - 1); i++)
    temp[i] = sten[0] * main[i-1] + sten[1] * main[i] + sten[2] * main[i+1];
  
  for(int i = 0; i < len; i++)
    main[i] = temp[i];

};











