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
#include "./data.h"
#include "./app_data.h"

#define Index3DShifted(_nx, _ny, _i, _j, _k, _s) ((_i+_s.x-1)+_nx*((_j+_s.y-1)+_ny*(_k+_s.z-1)))

AppDataVec::AppDataVec() {
  data_ = NULL;
}

AppDataVec::AppDataVec(const nimbus::GeometricRegion &global_reg,
                       const std::string& name,
                       bool make_proto)
    : global_region_(global_reg) {
    set_name(name);
    data_ = NULL;
    size_ = 0;
    if (make_proto)
        MakePrototype();
}

AppDataVec::AppDataVec(const nimbus::GeometricRegion &global_reg,
                       const nimbus::GeometricRegion &local_reg,
                       const std::string& name)
    : AppVar(local_reg),
      global_region_(global_reg),
      local_region_(local_reg) {
    set_name(name);
    shift_.x = global_region_.x() - local_region_.x();
    shift_.y = global_region_.y() - local_region_.y();
    shift_.z = global_region_.z() - local_region_.z();
    size_ = local_region_.dx() * local_region_.dy() * local_region_.dz();
    if (size_ > 0) {
      data_ = new double[size_];
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
    size_ = 0;
}

nimbus::AppVar* AppDataVec::CreateNew(const nimbus::GeometricRegion &local_reg) const {
    nimbus::AppVar* temp = new AppDataVec(global_region_,
                                          local_reg,
                                          name());
    return temp;
}

void AppDataVec::ReadAppData(const nimbus::DataArray &read_set,
                             const nimbus::GeometricRegion &read_reg) {
  size_t counter = 0;
  assert(local_region_ == object_region());
  nimbus::GeometricRegion final_read_reg =
    nimbus::GeometricRegion::GetIntersection(read_reg, local_region_);
  if (final_read_reg.dx() <= 0 || final_read_reg.dy() <= 0 || final_read_reg.dz() <= 0)
    return;

  DataArray::const_iterator iter = read_set.begin();
  for (; iter != read_set.end(); ++iter) {
    Vec *d = static_cast<Vec*>(*iter);
    GeometricRegion data_reg = d->region();
    nimbus::Coord data_shift;
    data_shift.x = global_region_.x() - data_reg.x();
    data_shift.y = global_region_.y() - data_reg.y();
    data_shift.z = global_region_.z() - data_reg.z();
    nimbus::GeometricRegion r =
      nimbus::GeometricRegion::GetIntersection(final_read_reg, data_reg);
    int_dimension_t xs = r.x(); int_dimension_t xe = r.x() + r.dx();
    int_dimension_t ys = r.y(); int_dimension_t ye = r.y() + r.dy();
    int_dimension_t zs = r.z(); int_dimension_t ze = r.z() + r.dz();
    for (int_dimension_t k = zs; k < ze; ++k) {
      int_dimension_t s_dz = r.dy() * (k + data_shift.z - 1);
      int_dimension_t d_dz = local_region_.dy() * (k + shift_.z - 1);
      for (int_dimension_t j = ys; j < ye; ++j) {
        int_dimension_t s_dyz = r.dx() * ((j + data_shift.y - 1) + s_dz);
        int_dimension_t d_dyz = local_region_.dx() * ((j + shift_.y - 1) + d_dz);
        for (int_dimension_t i = xs; i < xe; ++i) {
          data_[i + shift_.x - 1 + d_dyz] = d->data()[i + data_shift.x - 1 + s_dyz];
          // assert(size_ > Index3DShifted(local_region_.dx(), local_region_.dy(), i, j, k, shift_));
          // assert(d->size() > Index3DShifted(r.dx(), r.dy(), i, j, k, data_shift));
          // data_[Index3DShifted(local_region_.dx(), local_region_.dy(), i, j, k, shift_)] =
          //   d->data()[Index3DShifted(r.dx(), r.dy(), i, j, k, data_shift)];
          ++counter;
        }
      }
    }
  }
  dbg(DBG_APP, ">>> Reading %4i elements into app data for region %s\n", counter, read_reg.ToNetworkData().c_str());
}

void AppDataVec::WriteAppData(const nimbus::DataArray &write_set,
               const nimbus::GeometricRegion &write_reg) const {
  size_t counter = 0;
  assert(local_region_ == object_region());
  nimbus::GeometricRegion final_write_reg =
    nimbus::GeometricRegion::GetIntersection(write_reg, local_region_);
  if (final_write_reg.dx() <= 0 || final_write_reg.dy() <= 0 || final_write_reg.dz() <= 0)
    return;

  DataArray::const_iterator iter = write_set.begin();
  for (; iter != write_set.end(); ++iter) {
    Vec *d = static_cast<Vec*>(*iter);
    GeometricRegion data_reg = d->region();
    nimbus::Coord data_shift;
    data_shift.x = global_region_.x() - data_reg.x();
    data_shift.y = global_region_.y() - data_reg.y();
    data_shift.z = global_region_.z() - data_reg.z();
    nimbus::GeometricRegion r =
      nimbus::GeometricRegion::GetIntersection(final_write_reg, data_reg);
    int_dimension_t xs = r.x(); int_dimension_t xe = r.x() + r.dx();
    int_dimension_t ys = r.y(); int_dimension_t ye = r.y() + r.dy();
    int_dimension_t zs = r.z(); int_dimension_t ze = r.z() + r.dz();
    for (int_dimension_t k = zs; k < ze; ++k) {
      int_dimension_t d_dz = r.dy() * (k + data_shift.z - 1);
      int_dimension_t s_dz = local_region_.dy() * (k + shift_.z - 1);
      for (int_dimension_t j = ys; j < ye; ++j) {
        int_dimension_t d_dyz = r.dx() * ((j + data_shift.y - 1) + d_dz);
        int_dimension_t s_dyz = local_region_.dx() * ((j + shift_.y - 1) + s_dz);
        for (int_dimension_t i = xs; i < xe; ++i) {
          d->data()[i + data_shift.x - 1 + d_dyz] = data_[i + shift_.x - 1 + s_dyz];
          // assert(size_ > Index3DShifted(local_region_.dx(), local_region_.dy(), i, j, k, shift_));
          // assert(d->size() > Index3DShifted(r.dx(), r.dy(), i, j, k, data_shift));
          // d->data()[Index3DShifted(r.dx(), r.dy(), i, j, k, data_shift)] =
          //   data_[Index3DShifted(local_region_.dx(), local_region_.dy(), i, j, k, shift_)];
          counter++;
        }
      }
    }
  }
  dbg(DBG_APP, "<<< Writing %4i elements from app data for region %s\n", counter, write_reg.ToNetworkData().c_str());
}


