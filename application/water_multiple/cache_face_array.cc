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
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <string>

#include "application/water_multiple/cache_face_array.h"
#include "application/water_multiple/physbam_include.h"
#include "application/water_multiple/physbam_tools.h"
#include "data/cache/cache_object.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace application {

template<class T, class TS> CacheFaceArray<T, TS>::
CacheFaceArray(std::string type,
               const nimbus::GeometricRegion &local_region,
               const nimbus::GeometricRegion &global_region,
               int ghost_width)
    : CacheObject(type, local_region, global_region),
      ghost_width_(ghost_width),
      data_region_(local_region) {
      data_region_.Enlarge(ghost_width_);
      shift_.x = local_region.x() - global_region.x();
      shift_.y = local_region.y() - global_region.y();
      shift_.z = local_region.z() - global_region.z();
      Range domain = RangeFromRegions<TV>(local_region, global_region);
      TV_INT count = CountFromRegion(local_region);
      mac_grid.Initialize(count, domain, true);
      data_ = new PhysBAMFaceArray(mac_grid, ghost_width, false);
}

template<class T, class TS> void CacheFaceArray<T, TS>::
ReadToCache(const nimbus::DataSet &read_set) {
    Translator::template ReadFaceArray<T>(data_region_, shift_, read_set, data_);
}

template<class T, class TS> void CacheFaceArray<T, TS>::
WriteFromCache(const nimbus::DataSet &write_set) const {
    Translator::template WriteFaceArray<T>(data_region_, shift_, write_set, data_);
}

template<class T, class TS> nimbus::CacheObject *CacheFaceArray<T, TS>::
CreateNew() const {
    return new CacheFaceArray(type(),
                              local_region(),
                              global_region(),
                              ghost_width_);
}

template class CacheFaceArray<float, float>;
template class CacheFaceArray<bool, float>;

} // namespace application
