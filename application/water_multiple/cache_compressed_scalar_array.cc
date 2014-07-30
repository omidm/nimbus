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
#include "data/cache/cache_var.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

#include "application/water_multiple/cache_compressed_scalar_array.h"

namespace application {

template<class T> CacheCompressedScalarArray<T>::
CacheCompressedScalarArray(const nimbus::GeometricRegion &global_reg,
                           const int ghost_width,
                           bool make_proto)
    : global_region_(global_reg),
      ghost_width_(ghost_width) {
  data_ = new DataType();
  index_data_ = NULL;
  data_length_ = -1;
  if (make_proto)
    MakePrototype();
}

template<class T> CacheCompressedScalarArray<T>::
CacheCompressedScalarArray(const nimbus::GeometricRegion &global_reg,
                           const nimbus::GeometricRegion &ob_reg,
                           const int ghost_width)
    : CacheVar(ob_reg),
      global_region_(global_reg),
      local_region_(ob_reg.NewEnlarged(-ghost_width)),
      ghost_width_(ghost_width) {
  data_ = new DataType();
  index_data_ = NULL;
  data_length_ = -1;
  shift_.x = local_region_.x() - global_reg.x();
  shift_.y = local_region_.y() - global_reg.y();
  shift_.z = local_region_.z() - global_reg.z();
}

template<class T> nimbus::CacheVar *CacheCompressedScalarArray<T>::
CreateNew(const nimbus::GeometricRegion &ob_reg) const {
  return new CacheCompressedScalarArray(global_region_,
                                        ob_reg,
                                        ghost_width_);
}

template<class T> void CacheCompressedScalarArray<T>::
ReadToCache(const nimbus::DataArray &read_set,
            const nimbus::GeometricRegion &read_reg) {
  nimbus::GeometricRegion ob_reg = object_region();
  nimbus::GeometricRegion final_read_reg =
      nimbus::GeometricRegion::GetIntersection(read_reg, ob_reg);
  assert(final_read_reg.dx() > 0 && final_read_reg.dy() > 0 && final_read_reg.dz() > 0);
  // Loop through each element in read set, and fetch it to the cache object.
  Translator::template
      ReadCompressedScalarArray<T>(final_read_reg, shift_, read_set, data_,
                                   data_length_, *index_data_);
}

template<class T> void CacheCompressedScalarArray<T>::
WriteFromCache(const nimbus::DataArray &write_set,
               const nimbus::GeometricRegion &write_reg) const {
  if (write_reg.dx() <= 0 || write_reg.dy() <= 0 || write_reg.dz() <= 0)
    return;
  nimbus::GeometricRegion ob_reg = object_region();
  nimbus::GeometricRegion final_write_reg =
      nimbus::GeometricRegion::GetIntersection(write_reg, ob_reg);
  assert(final_write_reg.dx() > 0 && final_write_reg.dy() > 0 && final_write_reg.dz() > 0);
  // Loop through each element in write_set, look up the region using index, and
  // then write.
  Translator::template
      WriteCompressedScalarArray<T>(final_write_reg, shift_, write_set, *data_,
                                    data_length_, *index_data_);
}

template class CacheCompressedScalarArray<float>;

} // namespace application
