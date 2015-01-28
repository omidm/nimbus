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

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_VECTOR_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_VECTOR_H_

#include <string>

#include "application/water_multiple/parameters.h"
#include "application/water_multiple/physbam_include.h"
#include "data/app_data/app_var.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace application {

class AppDataVector : public nimbus::AppVar {
 public:
  typedef PhysBAM::VECTOR_ND<float> DATA_TYPE;
  explicit AppDataVector(const nimbus::GeometricRegion& global_reg,
                       bool make_proto,
                       const std::string& name);

  DATA_TYPE* data() {
    return data_;
  }
  void set_data(DATA_TYPE* d) {
    data_ = d;
  }

  virtual size_t memory_size() {
    return data_ ? sizeof(*this) + data_->memory_size() : sizeof(*this);
  }

 protected:
  explicit AppDataVector(const nimbus::GeometricRegion& global_reg,
                       const nimbus::GeometricRegion& ob_reg);

  virtual nimbus::AppVar* CreateNew(const nimbus::GeometricRegion &ob_reg) const;

  virtual void ReadAppData(const nimbus::DataArray& read_set,
                           const nimbus::GeometricRegion& read_reg);
  virtual void WriteAppData(const nimbus::DataArray& write_set,
                              const nimbus::GeometricRegion& write_reg) const;

 private:
  nimbus::GeometricRegion global_region_;
  nimbus::GeometricRegion local_region_;
  DATA_TYPE* data_;
};  // class AppDataVector

}  // namespace application

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_VECTOR_H_
