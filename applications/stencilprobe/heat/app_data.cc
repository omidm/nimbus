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
 * Data classes are implemented here.
 *
 */

#include <stdlib.h>
#include "./app_data.h"

AppDataVec::AppDataVec() {
  data_ = NULL;
}

AppDataVec::AppDataVec(const nimbus::GeometricRegion &global_reg,
                       const int ghost_width,
                       bool make_proto,
                       const std::string& name)
    : global_region_(global_reg),
      ghost_width_(ghost_width) {
    set_name(name);
    data_ = NULL;
    if (make_proto)
        MakePrototype();
}

AppDataVec::AppDataVec(const nimbus::GeometricRegion &global_reg,
                       const nimbus::GeometricRegion &ob_reg,
                       const int ghost_width,
                       const std::string& name)
    : AppVar(ob_reg),
      global_region_(global_reg),
      // local_region_(ob_reg.NewEnlarged(-ghost_width)),
      local_region_(ob_reg),
      ghost_width_(ghost_width) {
    set_name(name);
    shift_.x = local_region_.x() - global_reg.x();
    shift_.y = local_region_.y() - global_reg.y();
    shift_.z = local_region_.z() - global_reg.z();
    size_t size = local_region_.dx() * local_region_.dy() * local_region_.dz();
    if (size > 0) {
      data_ = new double[size];
    } else {
      data_ = NULL;
    }
}

AppDataVec::~AppDataVec() {
    Destroy();
}

void AppDataVec::Destroy() {
    if (data_) {
        delete data_;
        data_ = NULL;
    }
}

nimbus::AppVar* AppDataVec::CreateNew(const nimbus::GeometricRegion &ob_reg) const {
    nimbus::AppVar* temp = new AppDataVec(global_region_,
                                          ob_reg,
                                          ghost_width_,
                                          name());
    return temp;
}

void AppDataVec::ReadAppData(const nimbus::DataArray &read_set,
                             const nimbus::GeometricRegion &read_reg) {
    dbg(DBG_APP, "--- Reading %i elements into app data for region %s\n", read_set.size(), read_reg.ToNetworkData().c_str());
    nimbus::GeometricRegion ob_reg = object_region();
    nimbus::GeometricRegion final_read_reg =
        nimbus::GeometricRegion::GetIntersection(read_reg, ob_reg);
    assert(final_read_reg.dx() > 0 && final_read_reg.dy() > 0 && final_read_reg.dz() > 0);
    // Translator::template
    //     ReadFaceArray<T>(final_read_reg, local_region_, shift_, read_set, data_);
}

void AppDataVec::WriteAppData(const nimbus::DataArray &write_set,
               const nimbus::GeometricRegion &write_reg) const {
    dbg(DBG_APP, "--- Writing %i elements from app data for region %s\n", write_set.size(), write_reg.ToNetworkData().c_str());
    if (write_reg.dx() <= 0 || write_reg.dy() <= 0 || write_reg.dz() <= 0)
        return;
    nimbus::GeometricRegion ob_reg = object_region();
    nimbus::GeometricRegion final_write_reg =
        nimbus::GeometricRegion::GetIntersection(write_reg, ob_reg);
    assert(final_write_reg.dx() > 0 && final_write_reg.dy() > 0 && final_write_reg.dz() > 0);
    // Translator::template
    //     WriteFaceArray<T>(final_write_reg, shift_, write_set, data_);
}




