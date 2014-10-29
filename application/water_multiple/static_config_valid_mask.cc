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

#include "application/water_multiple/static_config_valid_mask.h"

namespace application {

StaticConfigValidMask::StaticConfigValidMask(
    const GeometricRegion& global_region) {
  global_region_ = global_region;
  physbam_data_structure_ = NULL;
}

StaticConfigValidMask::~StaticConfigValidMask() {
  if (physbam_data_structure_) {
    Destroy();
  }
}

StaticConfigVariable* StaticConfigValidMask::CreateNew(
    const GeometricRegion& local_region) const {
  StaticConfigValidMask* valid_mask = new StaticConfigValidMask(global_region_);
  valid_mask->local_region_ = local_region;
  // TODO
  // valid_mask->physbam_data_structure_->?//
  return valid_mask;
}

void StaticConfigValidMask::Destroy() {
  if (physbam_data_structure_) {
    delete physbam_data_structure_;
  }
}

StaticConfigValidMask::DataType* StaticConfigValidMask::GetData() const {
  return physbam_data_structure_;
}

}  // namespace application
