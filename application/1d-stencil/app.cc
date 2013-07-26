#include <iostream>
#include "app.h"

#define ML 4
#define GL 1

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

void App::load()
{

  std::cout << "Start Creating Data and Job Tables" << std::endl;
  
  registerJob("main", new Main(this, COMP));
  registerJob("init", new Init(this, COMP));
  registerJob("forloop", new ForLoop(this, COMP));
  registerJob("print", new Print(this, COMP));
  registerJob("applyLeft", new ApplyLeft(this, COMP));
  registerJob("applyRight", new ApplyRight(this, COMP));
  registerJob("updateLeft", new UpdateLeft(this, SYNC));
  registerJob("updateRight", new UpdateRight(this, SYNC));

  registerData("mainLeft", new Vec (ML));
  registerData("mainRight", new Vec (ML));
  registerData("ghostLeft", new Vec (GL));
  registerData("ghostRight", new Vec (GL));

  std::cout << "Finished Creating Data and Job Tables" << std::endl;

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

  /*
  initialize->after.insert(applyLeft);
  initialize->after.insert(applyRight);

  applyLeft->before.insert(initialize);
  applyLeft->after.insert(updateLeft);
  applyLeft->after.insert(updateRight);
  */
  // Continue like this
};

Main::Main(Application* app, JobType type)
  : Job(app, type) {
};

void Main::Execute(std::string params, const dataArray& da) {
  std::vector<int> jobIDs;
  std::vector<int> dataIDs;
  //getNewJobID(5, &jobIDs);
  //getNewDataID(4, &dataIDs);

  //defineData("mainLeft", dataIDs[0]);
  //defineData("mainRight", dataIDs[1]);
  //defineData("ghostLeft", dataIDs[2]);
  //defineData("ghostRight", dataIDs[3]);

  //spawnJob





};

Init::Init(Application* app, JobType type)
  : Job(app, type) {
};

void Init::Execute(std::string params, const dataArray& da)
{
  Vec *d = (Vec*)(da[0]);
  for(int i = 0; i < d->size ; i++)
    d->arr[i] = 0;
};


ForLoop::ForLoop(Application* app, JobType type)
  : Job(app, type) {
};

void ForLoop::Execute(std::string params, const dataArray& da) {
};

Print::Print(Application* app, JobType type)
  : Job(app, type) {
};

void Print::Execute(std::string params, const dataArray& da)
{
  Vec *d = (Vec*)(da[0]);
  for(int i = 0; i < d->size; i++)
    std::cout << d->arr[i] << ", ";
};


ApplyLeft::ApplyLeft(Application* app, JobType type)
  : Job(app, type) {
};

void ApplyLeft::Execute(std::string params, const dataArray& da)
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

ApplyRight::ApplyRight(Application* app, JobType type)
  : Job(app, type) {
};

void ApplyRight::Execute(std::string params, const dataArray& da)
{

};

UpdateLeft::UpdateLeft(Application* app, JobType type)
  : Job(app, type) {
};

void UpdateLeft::Execute(std::string params, const dataArray& da)
{
  Vec *d1 = (Vec*)(da[0]);
  Vec *d2 = (Vec*)(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[0];
}; 

UpdateRight::UpdateRight(Application* app, JobType type)
  : Job(app, type) {
};

void UpdateRight::Execute(std::string params, const dataArray& da)
{
};





