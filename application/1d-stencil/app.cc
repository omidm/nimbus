#include <iostream>
#include "app.h"

#define ML 4
#define GL 1

Vec::Vec(int size)
{
  this->size = size;
};

Data * Vec::Clone() {
  std::cout << "Cloning Vec data!\n";
  return new Vec(size);
};

void Vec::Create()
{
  arr = new int[size];
};

App::App()
{

};

void App::Load()
{

  std::cout << "Start Creating Data and Job Tables" << std::endl;
  
  RegisterJob("main", new Main(this, JOB_COMP));
  RegisterJob("init", new Init(this, JOB_COMP));
  RegisterJob("forLoop", new ForLoop(this, JOB_COMP));
  RegisterJob("print", new Print(this, JOB_COMP));
  RegisterJob("applyLeft", new ApplyLeft(this, JOB_COMP));
  RegisterJob("applyRight", new ApplyRight(this, JOB_COMP));
  RegisterJob("updateLeft", new UpdateLeft(this, JOB_SYNC));
  RegisterJob("updateRight", new UpdateRight(this, JOB_SYNC));

  RegisterData("mainLeft", new Vec (ML));
  RegisterData("mainRight", new Vec (ML));
  RegisterData("ghostLeft", new Vec (GL));
  RegisterData("ghostRight", new Vec (GL));

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
  return new Main(application(), type_);
};

void Main::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the main job\n";
  std::vector<int> j;
  std::vector<int> d;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_t> neighbor_partitions;
  partition_t partition_id = 0;
  std::string par;
  
  application()->GetNewJobID(5, &j);
  application()->GetNewDataID(4, &d);

  application()->DefineData("mainLeft", d[0], partition_id, neighbor_partitions, par);
  application()->DefineData("mainRight", d[1], partition_id, neighbor_partitions, par);
  application()->DefineData("ghostLeft", d[2], partition_id, neighbor_partitions, par);
  application()->DefineData("ghostRight", d[3], partition_id, neighbor_partitions, par);

  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[0]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application()->SpawnJob("init", j[0], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[1]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application()->SpawnJob("init", j[1], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[2]);
  write.clear(); write.insert(d[2]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application()->SpawnJob("init", j[2], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[3]);
  write.clear(); write.insert(d[3]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application()->SpawnJob("init", j[3], read, write, before, after, JOB_COMP, par);

  read.clear();
  write.clear();
  before.clear(); before.insert(j[0]); before.insert(j[1]); before.insert(j[2]); before.insert(j[3]);
  after.clear();
  // TODO: Load the "par" with the ids of four defined data instances.
  // TODO: Load the for loop couter and condition in "par"
  application()->SpawnJob("forLoop", j[4], read, write, before, after, JOB_COMP, par);
};

ForLoop::ForLoop(Application* app, JobType type)
  : Job(app, type) {
};

Job * ForLoop::Clone() {
  std::cout << "Cloning forLoop job!\n";
  return new ForLoop(application(), type_);
};

void ForLoop::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the forLoop job\n";
  std::vector<int> j;
  std::vector<int> d;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  std::string par;
  int counter = 0;
  int condition = 0;
  
  application()->GetNewJobID(7, &j);
  // TODO: Load "d" with the ids of data instances from the "params".
  d.push_back(1);
  d.push_back(2);
  d.push_back(3);
  d.push_back(4);

  // TODO: Load the for loop couter and condition from "params"; update counter.


  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[2]);
  before.clear();
  after.clear(); after.insert(j[2]); after.insert(j[3]);
  application()->SpawnJob("updateLeft", j[0], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[3]);
  before.clear();
  after.clear(); after.insert(j[2]); after.insert(j[3]);
  application()->SpawnJob("updateRight", j[1], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[0]); read.insert(d[3]);
  write.clear(); write.insert(d[0]);
  before.clear(); before.insert(j[0]); before.insert(j[1]);
  after.clear(); after.insert(j[4]);
  application()->SpawnJob("applyLeft", j[2], read, write, before, after, JOB_COMP, par);

  before.clear(); before.insert(j[0]); before.insert(j[1]);
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[1]); read.insert(d[2]);
  write.clear(); write.insert(d[1]);
  application()->SpawnJob("applyRight", j[3], read, write, before, after, JOB_COMP, par);

  if (counter > condition) {
  read.clear();
  write.clear();
  before.clear(); before.insert(j[2]); before.insert(j[3]);
  after.clear();
  // TODO: Load the "par" with the ids of four defined data instances.  
  // TODO: Load the for loop couter and condition in "par"
  application()->SpawnJob("forLoop", j[4], read, write, before, after, JOB_COMP, par);
  }
  else {
  before.clear(); before.insert(j[2]);
  after.clear();
  read.clear(); read.insert(d[0]);
  write.clear();
  application()->SpawnJob("print", j[5], read, write, before, after, JOB_COMP, par);

  before.clear(); before.insert(j[3]);
  after.clear();
  read.clear(); read.insert(d[1]);
  write.clear();
  application()->SpawnJob("print", j[6], read, write, before, after, JOB_COMP, par);
  }
};

Init::Init(Application* app, JobType type)
  : Job(app, type) {
};

Job * Init::Clone() {
  std::cout << "Cloning init job!\n";
  return new Init(application(), type_);
};

void Init::Execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the init job\n";
  Vec *d = (Vec*)(da[0]);
  for(int i = 0; i < d->size ; i++)
    d->arr[i] = 0;
};


Print::Print(Application* app, JobType type)
  : Job(app, type) {
};

Job * Print::Clone() {
  std::cout << "Cloning print job!\n";
  return new Print(application(), type_);
};

void Print::Execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the print job\n";
  Vec *d = (Vec*)(da[0]);
  for(int i = 0; i < d->size; i++)
    std::cout << d->arr[i] << ", ";
  std::cout << std::endl;
};


ApplyLeft::ApplyLeft(Application* app, JobType type)
  : Job(app, type) {
};

Job * ApplyLeft::Clone() {
  std::cout << "Cloning applyLeft job!\n";
  return new ApplyLeft(application(), type_);
};

void ApplyLeft::Execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the applyLeft job\n";
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

Job * ApplyRight::Clone() {
  std::cout << "Cloning applyRight job!\n";
  return new ApplyRight(application(), type_);
};

void ApplyRight::Execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the applyRight job\n";


};

UpdateLeft::UpdateLeft(Application* app, JobType type)
  : Job(app, type) {
};

Job * UpdateLeft::Clone() {
  std::cout << "Cloning updateLeft job!\n";
  return new UpdateLeft(application(), type_);
};

void UpdateLeft::Execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the updateLeft job\n";
  Vec *d1 = (Vec*)(da[0]);
  Vec *d2 = (Vec*)(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[0];
}; 

UpdateRight::UpdateRight(Application* app, JobType type)
  : Job(app, type) {
};

Job * UpdateRight::Clone() {
  std::cout << "Cloning updateRight job!\n";
  return new UpdateRight(application(), type_);
};

void UpdateRight::Execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the updateRight job\n";
  Vec *d1 = (Vec*)(da[0]);
  Vec *d2 = (Vec*)(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[0];
};





