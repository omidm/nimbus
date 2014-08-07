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
 * An Example application that spawns a lot of jobs in the 3d space.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include "./utils.h"
#include "protobuf_compiled/parameter_msg.pb.h"

using parameter_msg::ParameterMsg;

void LoadLogicalIdsInSet(Job* job,
    IDSet<logical_data_id_t>* set,
    const GeometricRegion& region, ...) {
  CLdoVector result;
  va_list vl;
  va_start(vl, region);
  char* arg = va_arg(vl, char*);
  while (arg != NULL) {
    job->GetCoveredLogicalObjects(&result,
        arg, &region);
    for (size_t i = 0; i < result.size(); ++i) {
      set->insert(result[i]->id());
      dbg(DBG_WORKER, "Loaded logical id %d of variable %s to the set.\n",
          result[i]->id(), arg);
    }
    arg = va_arg(vl,
        char*);
  }
  va_end(vl);
}

bool CompareData(Data* i, Data* j) {
  return i->region().x() < j->region().x();
}


void LoadDataFromNimbus(Job* job,
    const DataArray& da, std::vector<int>* result) {
  DataArray da_t;
  std::set<physical_data_id_t> seen;
  for (size_t i = 0; i < da.size(); ++i) {
    physical_data_id_t id = da[i]->physical_id();
    if ((job->read_set().contains(id)) && (seen.count(id) == 0)) {
      seen.insert(id);
      da_t.push_back(da[i]);
    }
  }
  std::sort(da_t.begin(), da_t.end(), CompareData);
  result->clear();
  for (size_t i = 0; i < da_t.size(); ++i) {
    if (i < (da_t.size() - 1)) {
      assert((da_t[i]->region().x() + da_t[i]->region().dx()) == (da_t[i + 1]->region().x()));
    }
    Vec *d = reinterpret_cast<Vec*>(da_t[i]);
    for (int j = 0; j < d->size(); ++j) {
      result->push_back(d->arr()[j]);
    }
  }
}

void SaveDataToNimbus(Job* job,
    const DataArray& da, std::vector<int>* vec) {
  DataArray da_t;
  std::set<physical_data_id_t> seen;
  for (size_t i = 0; i < da.size(); ++i) {
    physical_data_id_t id = da[i]->physical_id();
    if ((job->write_set().contains(id)) && (seen.count(id) == 0)) {
      seen.insert(id);
      da_t.push_back(da[i]);
    }
  }
  std::sort(da_t.begin(), da_t.end(), CompareData);
  size_t write_length = 0;
  size_t cursor = 0;
  for (size_t i = 0; i < da_t.size(); ++i) {
    if (i < (da_t.size() - 1)) {
      assert((da_t[i]->region().x() + da_t[i]->region().dx()) == (da_t[i + 1]->region().x()));
    }
    Vec *d = reinterpret_cast<Vec*>(da_t[i]);
    for (int j = 0; j < d->size(); ++j) {
      d->arr()[j] = vec->operator[](cursor++);
    }
    write_length += d->size();
  }

  assert(write_length == vec->size());
}

bool LoadParameter(Parameter *parameter, size_t *value) {
  std::string params_str(parameter->ser_data().data_ptr_raw(),
                         parameter->ser_data().size());
  ParameterMsg vec_msg;
  vec_msg.ParseFromString(params_str);
  assert(vec_msg.elem_size() == 1);
  *value = vec_msg.elem(0);
  return true;
}

bool SerializeParameter(Parameter *parameter, size_t value) {
  std::string params_str;
  ParameterMsg vec_msg;
  vec_msg.add_elem(value);
  vec_msg.SerializeToString(&params_str);
  parameter->set_ser_data(SerializedData(params_str));
  return true;
}




