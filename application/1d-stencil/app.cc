#include <iostream>
#include "app.h"

#define ML 4
#define GL 1

Vec::Vec(int size)
{
  this->size = size;
};

Data * Vec::clone() {
  std::cout << "Cloning Vec data!\n";
  return new Vec(size);
};

void Vec::create()
{
  arr = new int[size];
};

App::App()
{

};

void App::load()
{

  std::cout << "Start Creating Data and Job Tables" << std::endl;
  
  registerJob("main", new Main(this, JOB_COMP));
  registerJob("init", new Init(this, JOB_COMP));
  registerJob("forLoop", new ForLoop(this, JOB_COMP));
  registerJob("print", new Print(this, JOB_COMP));
  registerJob("applyLeft", new ApplyLeft(this, JOB_COMP));
  registerJob("applyRight", new ApplyRight(this, JOB_COMP));
  registerJob("updateLeft", new UpdateLeft(this, JOB_SYNC));
  registerJob("updateRight", new UpdateRight(this, JOB_SYNC));

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

Job * Main::clone() {
  std::cout << "Cloning main job!\n";
  return new Main(application_, type_);
};

void Main::execute(std::string params, const DataArray& da) {
  std::cout << "Executing the main job\n";
  std::vector<int> j;
  std::vector<int> d;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_t> neighbor_partitions;
  partition_t partition_id = 0;
  std::string par;
  
  application_->getNewJobID(5, &j);
  application_->getNewDataID(4, &d);

  application_->DefineData("mainLeft", d[0], partition_id, neighbor_partitions, par);
  application_->DefineData("mainRight", d[1], partition_id, neighbor_partitions, par);
  application_->DefineData("ghostLeft", d[2], partition_id, neighbor_partitions, par);
  application_->DefineData("ghostRight", d[3], partition_id, neighbor_partitions, par);

  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[0]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application_->SpawnJob("init", j[0], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[1]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application_->SpawnJob("init", j[1], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[2]);
  write.clear(); write.insert(d[2]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application_->SpawnJob("init", j[2], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[3]);
  write.clear(); write.insert(d[3]);
  before.clear();
  after.clear(); after.insert(j[4]);
  application_->SpawnJob("init", j[3], read, write, before, after, JOB_COMP, par);

  read.clear();
  write.clear();
  before.clear(); before.insert(j[0]); before.insert(j[1]); before.insert(j[2]); before.insert(j[3]);
  after.clear();
  // TODO: Load the "par" with the ids of four defined data instances.
  // TODO: Load the for loop couter and condition in "par"
  application_->SpawnJob("forLoop", j[4], read, write, before, after, JOB_COMP, par);
};

ForLoop::ForLoop(Application* app, JobType type)
  : Job(app, type) {
};

Job * ForLoop::clone() {
  std::cout << "Cloning forLoop job!\n";
  return new ForLoop(application_, type_);
};

void ForLoop::execute(std::string params, const DataArray& da) {
  std::cout << "Executing the forLoop job\n";
  std::vector<int> j;
  std::vector<int> d;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  std::string par;
  int counter = 0;
  int condition = 0;
  
  application_->getNewJobID(7, &j);
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
  application_->SpawnJob("updateLeft", j[0], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[3]);
  before.clear();
  after.clear(); after.insert(j[2]); after.insert(j[3]);
  application_->SpawnJob("updateRight", j[1], read, write, before, after, JOB_COMP, par);

  read.clear(); read.insert(d[0]); read.insert(d[3]);
  write.clear(); write.insert(d[0]);
  before.clear(); before.insert(j[0]); before.insert(j[1]);
  after.clear(); after.insert(j[4]);
  application_->SpawnJob("applyLeft", j[2], read, write, before, after, JOB_COMP, par);

  before.clear(); before.insert(j[0]); before.insert(j[1]);
  after.clear(); after.insert(j[4]);
  read.clear(); read.insert(d[1]); read.insert(d[2]);
  write.clear(); write.insert(d[1]);
  application_->SpawnJob("applyRight", j[3], read, write, before, after, JOB_COMP, par);

  if (counter > condition) {
  read.clear();
  write.clear();
  before.clear(); before.insert(j[2]); before.insert(j[3]);
  after.clear();
  // TODO: Load the "par" with the ids of four defined data instances.  
  // TODO: Load the for loop couter and condition in "par"
  application_->SpawnJob("forLoop", j[4], read, write, before, after, JOB_COMP, par);
  }
  else {
  before.clear(); before.insert(j[2]);
  after.clear();
  read.clear(); read.insert(d[0]);
  write.clear();
  application_->SpawnJob("print", j[5], read, write, before, after, JOB_COMP, par);

  before.clear(); before.insert(j[3]);
  after.clear();
  read.clear(); read.insert(d[1]);
  write.clear();
  application_->SpawnJob("print", j[6], read, write, before, after, JOB_COMP, par);
  }
};

Init::Init(Application* app, JobType type)
  : Job(app, type) {
};

Job * Init::clone() {
  std::cout << "Cloning init job!\n";
  return new Init(application_, type_);
};

void Init::execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the init job\n";
  Vec *d = (Vec*)(da[0]);
  for(int i = 0; i < d->size ; i++)
    d->arr[i] = 0;
};


Print::Print(Application* app, JobType type)
  : Job(app, type) {
};

Job * Print::clone() {
  std::cout << "Cloning print job!\n";
  return new Print(application_, type_);
};

void Print::execute(std::string params, const DataArray& da)
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

Job * ApplyLeft::clone() {
  std::cout << "Cloning applyLeft job!\n";
  return new ApplyLeft(application_, type_);
};

void ApplyLeft::execute(std::string params, const DataArray& da)
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

Job * ApplyRight::clone() {
  std::cout << "Cloning applyRight job!\n";
  return new ApplyRight(application_, type_);
};

void ApplyRight::execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the applyRight job\n";


};

UpdateLeft::UpdateLeft(Application* app, JobType type)
  : Job(app, type) {
};

Job * UpdateLeft::clone() {
  std::cout << "Cloning updateLeft job!\n";
  return new UpdateLeft(application_, type_);
};

void UpdateLeft::execute(std::string params, const DataArray& da)
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

Job * UpdateRight::clone() {
  std::cout << "Cloning updateRight job!\n";
  return new UpdateRight(application_, type_);
};

void UpdateRight::execute(std::string params, const DataArray& da)
{
  std::cout << "Executing the updateRight job\n";
  Vec *d1 = (Vec*)(da[0]);
  Vec *d2 = (Vec*)(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[0];
};





