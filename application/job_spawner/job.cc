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
 * Job classes for job spawner application.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include "./job.h"
#include "./data.h"

Main::Main(Application* app) {
  set_application(app);
};

Job * Main::Clone() {
  std::cout << "Cloning main job!\n";
  return new Main(application());
};

void Main::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the main job\n";

  std::vector<job_id_t> j;
  std::vector<logical_data_id_t> d;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  IDSet<partition_id_t> neighbor_partitions;
  GeometricRegion r_0(0, 0, 0, 1, 1, 0);
  GeometricRegion r_1(1, 0, 0, 3, 1, 0);
  GeometricRegion r_2(4, 0, 0, 1, 1, 0);
  GeometricRegion r_3(5, 0, 0, 1, 1, 0);
  GeometricRegion r_4(6, 0, 0, 3, 1, 0);
  GeometricRegion r_5(9, 0, 0, 1, 1, 0);
  ID<partition_id_t> p_0(0);
  ID<partition_id_t> p_1(1);
  ID<partition_id_t> p_2(2);
  ID<partition_id_t> p_3(3);
  ID<partition_id_t> p_4(4);
  ID<partition_id_t> p_5(5);
  Parameter par;
  IDSet<param_id_t> param_idset;

  GetNewJobID(&j, 7);
  GetNewLogicalDataID(&d, 8);

  DefinePartition(p_0, r_0, par);
  DefinePartition(p_1, r_1, par);
  DefinePartition(p_2, r_2, par);
  DefinePartition(p_3, r_3, par);
  DefinePartition(p_4, r_4, par);
  DefinePartition(p_5, r_5, par);

  DefineData("side", d[0], p_0.elem(), neighbor_partitions, par);
  DefineData("middle", d[1], p_1.elem(), neighbor_partitions, par);
  DefineData("side", d[2], p_2.elem(), neighbor_partitions, par);
  DefineData("side", d[3], p_3.elem(), neighbor_partitions, par);
  DefineData("middle", d[4], p_4.elem(), neighbor_partitions, par);
  DefineData("side", d[5], p_5.elem(), neighbor_partitions, par);

  read.clear(); read.insert(d[0]);
  write.clear(); write.insert(d[0]);
  before.clear();
  after.clear(); after.insert(j[6]);
  param_idset.clear(); param_idset.insert(0);
  par.set_idset(param_idset);
  SpawnComputeJob("init", j[0], read, write, before, after, par);

  read.clear(); read.insert(d[1]);
  write.clear(); write.insert(d[1]);
  before.clear();
  after.clear(); after.insert(j[6]);
  param_idset.clear(); param_idset.insert(0);
  par.set_idset(param_idset);
  SpawnComputeJob("init", j[1], read, write, before, after, par);

  read.clear(); read.insert(d[2]);
  write.clear(); write.insert(d[2]);
  before.clear();
  after.clear(); after.insert(j[6]);
  param_idset.clear(); param_idset.insert(ML - 1);
  par.set_idset(param_idset);
  SpawnComputeJob("init", j[2], read, write, before, after, par);

  read.clear(); read.insert(d[3]);
  write.clear(); write.insert(d[3]);
  before.clear();
  after.clear(); after.insert(j[6]);
  param_idset.clear(); param_idset.insert(0);
  par.set_idset(param_idset);
  SpawnComputeJob("init", j[3], read, write, before, after, par);

  read.clear(); read.insert(d[4]);
  write.clear(); write.insert(d[4]);
  before.clear();
  after.clear(); after.insert(j[6]);
  param_idset.clear(); param_idset.insert(1);
  par.set_idset(param_idset);
  SpawnComputeJob("init", j[4], read, write, before, after, par);

  read.clear(); read.insert(d[5]);
  write.clear(); write.insert(d[5]);
  before.clear();
  after.clear(); after.insert(j[6]);
  param_idset.clear(); param_idset.insert(0);
  par.set_idset(param_idset);
  SpawnComputeJob("init", j[5], read, write, before, after, par);

  read.clear();
  write.clear();
  before.clear(); before.insert(j[0]); before.insert(j[1]); before.insert(j[2]);
  before.insert(j[3]); before.insert(j[4]); before.insert(j[5]);
  after.clear();
  param_idset.clear();
  param_idset.insert(d[0]); param_idset.insert(d[1]);
  param_idset.insert(d[2]); param_idset.insert(d[3]);
  param_idset.insert(d[4]); param_idset.insert(d[5]);
  param_idset.insert(LOOP_COUNTER);
  par.set_idset(param_idset);
  SpawnComputeJob("forLoop", j[6], read, write, before, after, par);

  // *******************************************************************
  // Cheching to query  ldo_map with in the same job that defines the data
  // *******************************************************************
  if (GetLogicalObject(d[0]) == NULL) {
    dbg(DBG_TEMP, "Application: FAIL - did not find the ldo in ldo_map.\n");
  } else {
    dbg(DBG_TEMP, "Application: Found the ldo in ldo_map.\n");
  }
  // ****
};

ForLoop::ForLoop(Application* app) {
  set_application(app);
};

Job * ForLoop::Clone() {
  std::cout << "Cloning forLoop job!\n";
  return new ForLoop(application());
};

void ForLoop::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the forLoop job\n";

  // Check if partitions are synched or not.
  GeometricRegion r_temp;
  if (GetPartition(3, &r_temp)) {
    dbg(DBG_TEMP, "Application: Found the partition in ldo_map.\n");
  } else {
    dbg(DBG_TEMP, "Application: FAIL - did not find the partition in ldo_map.\n");
  }
  // *******************************************************************
  // Cheching the ldo_map and creation of PhysicalDataInstance from Data
  // *******************************************************************
  // logical_data_id_t l_id = da[0]->logical_id();
  // std::string l_name = da[0]->name();
  // Data* l_data = da[0];
  logical_data_id_t l_id = 100001;
  std::string l_name = "some_name";
  Data* l_data = new Data();
  if (GetLogicalObject(l_id) == NULL) {
    dbg(DBG_TEMP, "Application: FAIL - did not find the ldo in ldo_map.\n");
  } else {
    dbg(DBG_TEMP, "Application: Found the ldo in ldo_map.\n");
    GeometricRegion* r_t = new GeometricRegion(*(GetLogicalObject(l_id)->region()));
    LogicalDataObject ldo_t(l_id, l_name, r_t);
    PhysicalDataInstance pdi_t(l_id, &ldo_t, l_data, data_version_t(0));
  }
  // ****

  std::vector<job_id_t> j;
  std::vector<logical_data_id_t> d;
  IDSet<logical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter par;
  IDSet<param_id_t> param_idset;

  IDSet<param_id_t>::IDSetContainer::iterator it;
  IDSet<param_id_t> temp_set = params.idset();
  for (it = temp_set.begin(); it != temp_set.end(); it++) {
    d.push_back(*it);
  }

  if (d[6] > LOOP_CONDITION) {
    GetNewJobID(&j, 1);

    read.clear();
    write.clear();
    before.clear();
    after.clear();
    param_idset.clear();
    param_idset.insert(d[0]); param_idset.insert(d[1]);
    param_idset.insert(d[2]); param_idset.insert(d[3]);
    param_idset.insert(d[4]); param_idset.insert(d[5]);
    param_idset.insert(d[6] - 1);
    par.set_idset(param_idset);
    SpawnComputeJob("forLoop", j[0], read, write, before, after, par);
  } else {
    GetNewJobID(&j, 2);

    read.clear(); read.insert(d[1]); read.insert(d[2]);
    write.clear();
    before.clear(); before.insert(id().elem());
    after.clear();
    SpawnComputeJob("print", j[0], read, write, before, after, par);

    read.clear(); read.insert(d[3]); read.insert(d[4]);
    write.clear();
    before.clear(); before.insert(id().elem());
    after.clear();
    SpawnComputeJob("print", j[1], read, write, before, after, par);

    TerminateApplication();
  }
};

Init::Init() {
};

Job * Init::Clone() {
  std::cout << "Cloning init job!\n";
  return new Init();
};

void Init::Execute(Parameter params, const DataArray& da) {
  uint32_t base_val;
  base_val = *(params.idset().begin());
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

void Print::Execute(Parameter params, const DataArray& da) {
  std::cout << "Executing the print job\n";
  Vec *d1 = reinterpret_cast<Vec*>(da[0]);
  Vec *d2 = reinterpret_cast<Vec*>(da[1]);
  std::cout << "OUTPUT: ";
  for (int i = 0; i < d1->size(); i++)
    std::cout << d1->arr()[i] << ", ";
  for (int i = 0; i < d2->size(); i++)
    std::cout << d2->arr()[i] << ", ";
  std::cout << std::endl;
};




