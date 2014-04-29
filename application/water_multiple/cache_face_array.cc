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
               const nimbus::GeometricRegion &global_region,
               const int ghost_width,
               const nimbus::GeometricRegion &app_region)
    : CacheObject(type, app_region),
      ghost_width_(ghost_width),
      global_region_(global_region),
      local_region_(app_region.NewEnlarged(-ghost_width_)) {
      shift_.x = local_region_.x() - global_region.x();
      shift_.y = local_region_.y() - global_region.y();
      shift_.z = local_region_.z() - global_region.z();
      if (local_region_.dx() > 0 && local_region_.dy() > 0 && local_region_.dz() > 0) {
        Range domain = RangeFromRegions<TV>(global_region, local_region_);
        TV_INT count = CountFromRegion(local_region_);
        mac_grid_.Initialize(count, domain, true);
        data_ = new PhysBAMFaceArray(mac_grid_, ghost_width, false);
      }
}

template<class T, class TS> void CacheFaceArray<T, TS>::
ReadToCache(const nimbus::DataArray &read_set,
            const nimbus::GeometricRegion &reg) {
    //dbg(DBG_WARN, "\n--- Reading %i elements into face array for region %s\n", read_set.size(), reg.toString().c_str());
    Translator::template ReadFaceArray<T>(reg, local_region_, shift_, read_set, data_);
}

template<class T, class TS> void CacheFaceArray<T, TS>::
WriteFromCache(const nimbus::DataArray &write_set,
               const nimbus::GeometricRegion &reg) const {
    //dbg(DBG_WARN, "\n Writing %i elements into face array for region %s\n", write_set.size(), reg.toString().c_str());
    Translator::template WriteFaceArray<T>(reg, shift_, write_set, data_);
}

template<class T, class TS> nimbus::CacheObject *CacheFaceArray<T, TS>::
CreateNew(const nimbus::GeometricRegion &ar) const {
    return new CacheFaceArray(type(),
                              global_region_,
                              ghost_width_,
                              ar);
}

template class CacheFaceArray<float, float>;
template class CacheFaceArray<bool, float>;

} // namespace application
