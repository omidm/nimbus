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
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <string>

#include "application/water_multiple/physbam_include.h"
#include "application/water_multiple/physbam_tools.h"
#include "data/app_data/app_var.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

#include "application/water_multiple/projection/data_raw_vector_nd.h"
#include "application/water_multiple/projection/app_data_vector.h"

namespace application {

AppDataVector::AppDataVector() {
  data_ = NULL;
}

AppDataVector::AppDataVector(const nimbus::GeometricRegion &global_reg,
                         bool make_proto,
                         const std::string& name)
    : global_region_(global_reg) {
  set_name(name);
  data_ = NULL;
  if (make_proto)
    MakePrototype();
}

AppDataVector::AppDataVector(const nimbus::GeometricRegion &global_reg,
                         const nimbus::GeometricRegion &ob_reg)
    : AppVar(ob_reg), global_region_(global_reg), local_region_(ob_reg) {
  data_ = new DATA_TYPE;
}

AppDataVector::~AppDataVector() {
  Destroy();
}

void AppDataVector::Destroy() {
  if (data_)
    delete data_;
}

nimbus::AppVar* AppDataVector::CreateNew(
    const nimbus::GeometricRegion &ob_reg) const {
  nimbus::AppVar* temp = new AppDataVector(global_region_, ob_reg);
  temp->set_name(name());
  return temp;
}

void AppDataVector::ReadAppData(const nimbus::DataArray &read_set,
                              const nimbus::GeometricRegion &read_reg) {
  if (read_set.size() == 0) {
    return;
  }
  assert(read_set.size() == 1);
  assert(read_set[0]->region().IsEqual(&read_reg));
  assert(object_region().IsEqual(&read_reg));
  assert(dynamic_cast<DataRawVectorNd*>(read_set[0]));
  dynamic_cast<DataRawVectorNd*>(read_set[0])->LoadFromNimbus(data_);
}

void AppDataVector::WriteAppData(
    const nimbus::DataArray &write_set,
    const nimbus::GeometricRegion &write_reg) const {
  if (write_set.size() == 0) {
    return;
  }
  if (write_reg.dx() <= 0 || write_reg.dy() <= 0 || write_reg.dz() <= 0)
    return;
  assert(write_set.size() == 1);
  assert(write_set[0]->region().IsEqual(&write_reg));
  assert(object_region().IsEqual(&write_reg));
  assert(dynamic_cast<DataRawVectorNd*>(write_set[0]));
  dynamic_cast<DataRawVectorNd*>(write_set[0])->SaveToNimbus(*data_);
}

} // namespace application
