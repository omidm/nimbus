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
 * An Example application that is meant to run over multiple workers.
 * It is simply applying a stencil over a one dimensional array.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include "./app.h"

using vector_msg::VectorMsg;

Vec::Vec(int size) {
  size_ = size;
};

Vec::~Vec() {
};

Data * Vec::Clone() {
  std::cout << "Cloning Vec data!\n";
  return new Vec(size_);
};

void Vec::Create() {
  arr_ = new int[size_];
};

void Vec::Destroy() {
  delete arr_;
};

void Vec::Copy(Data* from) {
  Vec *d = reinterpret_cast<Vec*>(from);
  for (int i = 0; i < size_; i++)
    arr_[i] = d->arr()[i];
}

bool Vec::Serialize(SerializedData* ser_data) {
  VectorMsg vec_msg;
  for (int i = 0; i < size_; i++)
    vec_msg.add_elem(arr_[i]);
  std::string str;
  vec_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());
  return true;
}

bool Vec::DeSerialize(const SerializedData& ser_data, Data** result) {
  VectorMsg vec_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  vec_msg.ParseFromString(str);
  Vec* vec = new Vec(size_);
  vec->Create();
  for (int i = 0; (i < size_) && (i < vec_msg.elem_size()); i++)
     vec->arr()[i] = vec_msg.elem(i);

  *result = vec;
  return true;
}

int Vec::size() {
  return size_;
}

int* Vec::arr() {
  return arr_;
}

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

  RegisterData("middle", new Vec(ML - 1));
  RegisterData("side", new Vec(GL));

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
  partition_t p_1 = 1;
  partition_t p_2 = 2;
  std::string par;

  GetNewJobID(&j, 7);
  GetNewDataID(&d, 8);

  DefineData("side", d[0], p_1, neighbor_partitions, par);
  DefineData("middle", d[1], p_1, neighbor_partitions, par);
  DefineData("side", d[2], p_1, neighbor_partitions, par);
  DefineData("side", d[3], p_2, neighbor_partitions, par);
  DefineData("middle", d[4], p_2, neighbor_partitions, par);
  DefineData("side", d[5], p_2, neighbor_partitions, par);
  DefineData("side", d[6], p_2, neighbor_partitions, par);
  DefineData("side", d[7], p_1, neighbor_partitions, par);

  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[0]);
  before.clear();
  after.clear(); after.insert(j[6]);
  par = ID<data_id_t>(0).toString();
  SpawnComputeJob("init", j[0], read, write, before, after, par);

  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[1]);
  before.clear();
  after.clear(); after.insert(j[6]);
  par = ID<data_id_t>(0).toString();
  SpawnComputeJob("init", j[1], read, write, before, after, par);

  read.clear(); read.insert(d[2]);
  write.clear(); write.insert(d[2]);
  before.clear();
  after.clear(); after.insert(j[6]);
  par = ID<data_id_t>(ML - 1).toString();
  SpawnComputeJob("init", j[2], read, write, before, after, par);

  read.clear(); read.insert(d[3]);
  write.clear(); write.insert(d[3]);
  before.clear();
  after.clear(); after.insert(j[6]);
  par = ID<data_id_t>(0).toString();
  SpawnComputeJob("init", j[3], read, write, before, after, par);

  read.clear(); read.insert(d[4]);
  write.clear(); write.insert(d[4]);
  before.clear();
  after.clear(); after.insert(j[6]);
  par = ID<data_id_t>(1).toString();
  SpawnComputeJob("init", j[4], read, write, before, after, par);

  read.clear(); read.insert(d[5]);
  write.clear(); write.insert(d[5]);
  before.clear();
  after.clear(); after.insert(j[6]);
  par = ID<data_id_t>(0).toString();
  SpawnComputeJob("init", j[5], read, write, before, after, par);

  read.clear();
  write.clear();
  before.clear(); before.insert(j[0]); before.insert(j[1]); before.insert(j[2]);
  before.insert(j[3]); before.insert(j[4]); before.insert(j[5]);
  after.clear();
  IDSet<data_id_t> temp_set;
  temp_set.insert(d[0]); temp_set.insert(d[1]); temp_set.insert(d[2]); temp_set.insert(d[3]);
  temp_set.insert(d[4]); temp_set.insert(d[5]); temp_set.insert(d[6]); temp_set.insert(d[7]);
  par = temp_set.toString();
  par += ("-" + ID<data_id_t>(LOOP_COUNTER).toString());
  SpawnComputeJob("forLoop", j[6], read, write, before, after, par);
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

  IDSet<data_id_t>::IDSetContainer temp_set;
  nimbus::ParseIDSet(*iter, temp_set);
  IDSet<data_id_t>::IDSetContainer::iterator it;
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

    before.clear();
    after.clear(); after.insert(j[2]); after.insert(j[3]);
    SpawnCopyJob(j[0], d[2], d[6], before, after, par);

    before.clear();
    after.clear(); after.insert(j[2]); after.insert(j[3]);
    SpawnCopyJob(j[1], d[3], d[7], before, after, par);

    read.clear(); read.insert(d[0]); read.insert(d[1]); read.insert(d[2]); read.insert(d[7]);
    write.clear(); write.insert(d[1]); write.insert(d[2]);
    before.clear(); before.insert(j[0]); before.insert(j[1]);
    after.clear(); after.insert(j[4]);
    SpawnComputeJob("applyLeft", j[2], read, write, before, after, par);

    read.clear(); read.insert(d[3]); read.insert(d[4]); read.insert(d[5]); read.insert(d[6]);
    write.clear(); write.insert(d[3]); write.insert(d[4]);
    before.clear(); before.insert(j[0]); before.insert(j[1]);
    after.clear(); after.insert(j[4]);
    SpawnComputeJob("applyRight", j[3], read, write, before, after, par);

    read.clear();
    write.clear();
    before.clear(); before.insert(j[2]); before.insert(j[3]);
    after.clear();
    IDSet<data_id_t> temp_set;
    temp_set.insert(d[0]); temp_set.insert(d[1]); temp_set.insert(d[2]); temp_set.insert(d[3]);
    temp_set.insert(d[4]); temp_set.insert(d[5]); temp_set.insert(d[6]); temp_set.insert(d[7]);
    par = temp_set.toString();
    par += ("-" + ID<data_id_t>(counter - 1).toString());
    SpawnComputeJob("forLoop", j[4], read, write, before, after, par);
  } else {
    GetNewJobID(&j, 2);

    read.clear(); read.insert(d[1]); read.insert(d[2]);
    write.clear();
    before.clear();
    after.clear();
    SpawnComputeJob("print", j[0], read, write, before, after, par);

    read.clear(); read.insert(d[3]); read.insert(d[4]);
    write.clear();
    before.clear();
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
  uint32_t base_val;
  nimbus::ParseID(params, base_val);
  std::cout << "Executing the init job\n";
  Vec *d = reinterpret_cast<Vec*>(da[0]);
  for (int i = 0; i < d->size() ; i++)
    d->arr()[i] = base_val + i;
};


Print::Print() {
};

Job * Print::Clone() {
  std::cout << "Cloning print job!\n";
  return new Print();
};

void Print::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the print job\n";
  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);
  for (int i = 0; i < d1->size(); i++)
    std::cout << d1->arr()[i] << ", ";
  for (int i = 0; i < d2->size(); i++)
    std::cout << d2->arr()[i] << ", ";
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
  int sten[] = {-1, +2, -1};

  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);
  Vec *d3 = reinterpret_cast<Vec*>(da[2]);
  Vec *d4 = reinterpret_cast<Vec*>(da[3]);

  int len = d2->size();
  int temp_middle[len]; // NOLINT

  temp_middle[0] = sten[0] * d1->arr()[0] + sten[1] * d2->arr()[0] + sten[2] * d2->arr()[1];
  temp_middle[len-1] = sten[0] * d2->arr()[len-2] + sten[1] * d2->arr()[len-1] + sten[2] * d3->arr()[0]; // NOLINT
  for (int i = 1; i < (len - 1); i++)
    temp_middle[i] = sten[0] * d2->arr()[i-1] + sten[1] * d2->arr()[i] + sten[2] * d2->arr()[i+1];

  d3->arr()[0] =  sten[0] * d2->arr()[len-1] + sten[1] * d3->arr()[0] + sten[2] * d4->arr()[0];

  for (int i = 0; i < len; i++)
    d2->arr()[i] = temp_middle[i];
};

ApplyRight::ApplyRight() {
};

Job * ApplyRight::Clone() {
  std::cout << "Cloning applyRight job!\n";
  return new ApplyRight();
};

void ApplyRight::Execute(std::string params, const DataArray& da) {
  std::cout << "Executing the applyRight job\n";
  int sten[] = {-1, +2, -1};

  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);
  Vec *d3 = reinterpret_cast<Vec*>(da[2]);
  Vec *d4 = reinterpret_cast<Vec*>(da[3]);

  int len = d2->size();
  int temp_middle[len]; // NOLINT

  temp_middle[0] = sten[0] * d1->arr()[0] + sten[1] * d2->arr()[0] + sten[2] * d2->arr()[1];
  temp_middle[len-1] = sten[0] * d2->arr()[len-2] + sten[1] * d2->arr()[len-1] + sten[2] * d3->arr()[0]; // NOLINT
  for (int i = 1; i < (len - 1); i++)
    temp_middle[i] = sten[0] * d2->arr()[i-1] + sten[1] * d2->arr()[i] + sten[2] * d2->arr()[i+1];

  d1->arr()[0] =  sten[0] * d4->arr()[0] + sten[1] * d1->arr()[0] + sten[2] * d2->arr()[0];

  for (int i = 0; i < len; i++)
    d2->arr()[i] = temp_middle[i];
};




