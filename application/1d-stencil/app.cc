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
  registerJob("forLoop", new ForLoop(this, COMP));
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

Job * Main::Clone() {
  std::cout << "Cloning main job!\n";
  return new Main(app, type);
};

void Main::Execute(std::string params, const dataArray& da) {
  std::cout << "Executing the main job\n";
  std::vector<int> j;
  std::vector<int> d;
  IDSet before, after, read, write;
  std::string par = "none";
  
  app->getNewJobID(5, &j);
  app->getNewDataID(4, &d);

  app->defineData("mainLeft", d[0]);
  app->defineData("mainRight", d[1]);
  app->defineData("ghostLeft", d[2]);
  app->defineData("ghostRight", d[3]);

  before.clear();
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[0]);
  app->spawnJob("init", j[0], before, after, read, write, par);

  before.clear();
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[1]);
  app->spawnJob("init", j[1], before, after, read, write, par);

  before.clear();
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[2]);
  write.clear(); write.insert(d[2]);
  app->spawnJob("init", j[2], before, after, read, write, par);

  before.clear();
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[3]);
  write.clear(); write.insert(d[3]);
  app->spawnJob("init", j[3], before, after, read, write, par);

  before.clear(); before.insert(j[0]); before.insert(j[1]); before.insert(j[2]); before.insert(j[3]);
  after.clear();
  read.clear();
  write.clear();
  // TODO: Load the "par" with the ids of four defined data instances.
  // TODO: Load the for loop couter and condition in "par"
  app->spawnJob("forLoop", j[4], before, after, read, write, par);
};

ForLoop::ForLoop(Application* app, JobType type)
  : Job(app, type) {
};

void ForLoop::Execute(std::string params, const dataArray& da) {
  std::vector<int> j;
  std::vector<int> d;
  IDSet before, after, read, write;
  std::string par;
  int counter = 0;
  int condition = 0;
  
  app->getNewJobID(7, &j);
  // TODO: Load "d" with the ids of data instances from the "params".
  // TODO: Load the for loop couter and condition from "params"; update counter.

  before.clear();
  after.clear(); after.insert(j[2]); after.insert(j[3]);
  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[2]);
  app->spawnJob("updateLeft", j[0], before, after, read, write, par);

  before.clear();
  after.clear(); after.insert(j[2]); after.insert(j[3]);
  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[3]);
  app->spawnJob("updateRight", j[1], before, after, read, write, par);

  before.clear(); before.insert(j[0]); before.insert(j[1]);
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[0]); read.insert(d[3]);
  write.clear(); write.insert(d[0]);
  app->spawnJob("applyLeft", j[2], before, after, read, write, par);

  before.clear(); before.insert(j[0]); before.insert(j[1]);
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[1]); read.insert(d[2]);
  write.clear(); write.insert(d[1]);
  app->spawnJob("applyLeft", j[3], before, after, read, write, par);

  if (counter > condition) {
  before.clear(); before.insert(j[2]); before.insert(j[3]);
  after.clear();
  read.clear();
  write.clear();
  // TODO: Load the "par" with the ids of four defined data instances.  
  // TODO: Load the for loop couter and condition in "par"
  app->spawnJob("forLoop", j[4], before, after, read, write, par);
  }
  else {
  before.clear(); before.insert(j[2]);
  after.clear();
  read.clear(); read.insert(d[0]);
  write.clear();
  app->spawnJob("applyLeft", j[5], before, after, read, write, par);

  before.clear(); before.insert(j[3]);
  after.clear();
  read.clear(); read.insert(d[1]);
  write.clear();
  app->spawnJob("applyLeft", j[6], before, after, read, write, par);
  }
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
  Vec *d1 = (Vec*)(da[0]);
  Vec *d2 = (Vec*)(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[0];
};





