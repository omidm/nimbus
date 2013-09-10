/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * An Example application that is meant to run over single worker.
 * It is simply applying a stencil over a one dimensional array.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include "./app.h"

Vec::Vec(int size) {
  this->size = size;
};

Data * Vec::Clone() {
  std::cout << "Cloning Vec data!\n";
  return new Vec(size);
};

void Vec::Create() {
  arr = new int[size];
};

App::App() {
};

void App::Load() {
  std::cout << "Start Creating Data and Job Tables" << std::endl;

  RegisterJob("main", new Main(this));
  RegisterJob("init", new Init());
  RegisterJob("forLoop", new ForLoop(this));
  RegisterJob("print", new Print());
  RegisterJob("applyLeft", new ApplyLeft());
  RegisterJob("applyRight", new ApplyRight());
  RegisterJob("updateLeft", new UpdateLeft());
  RegisterJob("updateRight", new UpdateRight());

  RegisterData("main", new Vec(ML));
  RegisterData("ghost", new Vec(GL));

  std::cout << "Finished Creating Data and Job Tables" << std::endl;
};

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  std::cout << "Cloning main job!\n";
  return new Main(application());
};

void Main::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the main job\n";
  std::vector<job_id_t> j;
  std::vector<data_id_t> d;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_t> neighbor_partitions;
  partition_t partition_id = 0;
  std::string par;

  GetNewJobID(&j, 7);
  GetNewDataID(&d, 4);

  DefineData("main", d[0], partition_id, neighbor_partitions, par);
  DefineData("main", d[1], partition_id, neighbor_partitions, par);
  DefineData("ghost", d[2], partition_id, neighbor_partitions, par);
  DefineData("ghost", d[3], partition_id, neighbor_partitions, par);

  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[0]);
  before.clear();
  after.clear(); after.insert(j[4]);
  SpawnComputeJob("init", j[0], read, write, before, after, par);

  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[1]);
  before.clear();
  after.clear(); after.insert(j[4]);
  SpawnComputeJob("init", j[1], read, write, before, after, par);

  read.clear(); read.insert(d[2]);
  write.clear(); write.insert(d[2]);
  before.clear();
  after.clear(); after.insert(j[4]);
  SpawnComputeJob("init", j[2], read, write, before, after, par);

  read.clear(); read.insert(d[3]);
  write.clear(); write.insert(d[3]);
  before.clear();
  after.clear(); after.insert(j[4]);
  SpawnComputeJob("init", j[3], read, write, before, after, par);

  read.clear();
  write.clear();
  before.clear(); before.insert(j[0]); before.insert(j[1]); before.insert(j[2]); before.insert(j[3]); // NOLINT
  after.clear();
  IDSet<data_id_t> temp_set;
  temp_set.insert(d[0]); temp_set.insert(d[1]); temp_set.insert(d[2]); temp_set.insert(d[3]);
  ID<data_id_t> counter(LOOP_COUNTER);
  par = temp_set.toString();
  par += ("-" + counter.toString());
  SpawnComputeJob("forLoop", j[4], read, write, before, after, par);
};

ForLoop::ForLoop(Application* app) {
  set_application(app);
};

Job * ForLoop::Clone() {
  std::cout << "Cloning forLoop job!\n";
  return new ForLoop(application());
};

void ForLoop::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the forLoop job\n";
  std::vector<job_id_t> j;
  std::vector<data_id_t> d;
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  std::string par;




  char_separator<char> separator("-");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();

  std::set<data_id_t> temp_set;
  nimbus::ParseIDSet(*iter, temp_set);
  std::set<data_id_t>::iterator it;
  for (it = temp_set.begin(); it != temp_set.end(); it++) {
    d.push_back(*it);
  }

  iter++;
  uint32_t counter;
  nimbus::ParseID(*iter, counter);

//  d.push_back(16777217);
//  d.push_back(16777218);
//  d.push_back(16777219);
//  d.push_back(16777220);


  if (counter > LOOP_CONDITION) {
    GetNewJobID(&j, 5);

    read.clear(); read.insert(d[1]);
    write.clear(); write.insert(d[3]);
    before.clear();
    after.clear(); after.insert(j[2]); after.insert(j[3]);
    SpawnComputeJob("updateLeft", j[0], read, write, before, after, par);

    read.clear(); read.insert(d[0]);
    write.clear(); write.insert(d[2]);
    before.clear();
    after.clear(); after.insert(j[2]); after.insert(j[3]);
    SpawnComputeJob("updateRight", j[1], read, write, before, after, par);

    read.clear(); read.insert(d[0]); read.insert(d[3]);
    write.clear(); write.insert(d[0]);
    before.clear(); before.insert(j[0]); before.insert(j[1]);
    after.clear(); after.insert(j[4]);
    SpawnComputeJob("applyLeft", j[2], read, write, before, after, par);

    read.clear(); read.insert(d[1]); read.insert(d[2]);
    write.clear(); write.insert(d[1]);
    before.clear(); before.insert(j[0]); before.insert(j[1]);
    after.clear(); after.insert(j[4]);
    SpawnComputeJob("applyRight", j[3], read, write, before, after, par);

    read.clear();
    write.clear();
    before.clear(); before.insert(j[2]); before.insert(j[3]);
    after.clear();
    IDSet<data_id_t> temp_set;
    temp_set.insert(d[0]); temp_set.insert(d[1]); temp_set.insert(d[2]); temp_set.insert(d[3]);
    ID<data_id_t> temp_counter(counter - 1);
    par = temp_set.toString();
    par += ("-" + temp_counter.toString());
    SpawnComputeJob("forLoop", j[4], read, write, before, after, par);
  } else {
    GetNewJobID(&j, 2);

    read.clear(); read.insert(d[0]);
    write.clear();
    before.clear();
    after.clear(); after.insert(j[1]);
    SpawnComputeJob("print", j[0], read, write, before, after, par);

    read.clear(); read.insert(d[1]);
    write.clear();
    before.clear(); before.insert(j[0]);
    after.clear();
    SpawnComputeJob("print", j[1], read, write, before, after, par);
  }
};

Init::Init() {
};

Job * Init::Clone() {
  std::cout << "Cloning init job!\n";
  return new Init();
};

void Init::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the init job\n";
  Vec *d = reinterpret_cast<Vec*>(da[0]);
  for (int i = 0; i < d->size ; i++)
    d->arr[i] = i;
};


Print::Print() {
};

Job * Print::Clone() {
  std::cout << "Cloning print job!\n";
  return new Print();
};

void Print::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the print job\n";
  Vec *d = reinterpret_cast<Vec*>(da[0]);
  for (int i = 0; i < d->size; i++)
    std::cout << d->arr[i] << ", ";
  std::cout << std::endl;
};


ApplyLeft::ApplyLeft() {
};

Job * ApplyLeft::Clone() {
  std::cout << "Cloning applyLeft job!\n";
  return new ApplyLeft();
};

void ApplyLeft::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the applyLeft job\n";
  int sten[] = {-1, +2, -1, 0};

  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);

  int len = d1->size;
  int *main = d1->arr;
  int *ghost = d2->arr;
  int temp[len]; // NOLINT

  temp[0] = sten[0] * sten[3] + sten[1] * main[0] + sten[2] * main[1];
  temp[len-1] = sten[0] * main[len-2] + sten[1] * main[len-1] + sten[2] * ghost[0];
  for (int i = 1; i < (len - 1); i++)
    temp[i] = sten[0] * main[i-1] + sten[1] * main[i] + sten[2] * main[i+1];

  for (int i = 0; i < len; i++)
    main[i] = temp[i];
};

ApplyRight::ApplyRight() {
};

Job * ApplyRight::Clone() {
  std::cout << "Cloning applyRight job!\n";
  return new ApplyRight();
};

void ApplyRight::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the applyRight job\n";
  int sten[] = {-1, +2, -1, 0};

  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);

  int len = d1->size;
  int *main = d1->arr;
  int *ghost = d2->arr;
  int temp[len]; // NOLINT

  temp[0] = sten[0] * ghost[0] + sten[1] * main[0] + sten[2] * main[1];
  temp[len-1] = sten[0] * main[len-2] + sten[1] * main[len-1] + sten[2] * sten[3];
  for (int i = 1; i < (len - 1); i++)
    temp[i] = sten[0] * main[i-1] + sten[1] * main[i] + sten[2] * main[i+1];

  for (int i = 0; i < len; i++)
    main[i] = temp[i];
};

UpdateLeft::UpdateLeft() {
};

Job * UpdateLeft::Clone() {
  std::cout << "Cloning updateLeft job!\n";
  return new UpdateLeft();
};

void UpdateLeft::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the updateLeft job\n";
  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[0];
};

UpdateRight::UpdateRight() {
};

Job * UpdateRight::Clone() {
  std::cout << "Cloning updateRight job!\n";
  return new UpdateRight();
};

void UpdateRight::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the updateRight job\n";
  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);

  int *main = d1->arr;
  int *ghost = d2->arr;
  ghost[0] = main[d1->size - 1];
};





