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

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_COMPRESSED_SCALAR_ARRAY_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_COMPRESSED_SCALAR_ARRAY_H_

#include <string>

#include "application/water_multiple/physbam_include.h"
#include "data/app_data/app_var.h"
#include "data/physbam/translator_physbam.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace application {

template<class T>
class AppDataCompressedScalarArray : public nimbus::AppVar {
  typedef PhysBAM::VECTOR_ND<float> DataType;
  // Cell index to matrix index.
  typedef typename PhysBAM::VECTOR<int, 3> TV_INT;
  typedef PhysBAM::ARRAY<int, TV_INT> IndexType;
  typedef typename nimbus::TranslatorPhysBAM<float> Translator;

 public:
  explicit AppDataCompressedScalarArray(const nimbus::GeometricRegion &global_reg,
                                      const int ghost_width,
                                      bool make_proto,
                                      const std::string& name);

  DataType* data() { return data_; }
  void set_data(DataType* d) { data_ = d; }
  IndexType* index_data() { return index_data_; }
  void set_index_data(IndexType* d);
  nimbus::int_dimension_t data_length() { return data_length_; }
  void set_data_length(nimbus::int_dimension_t l) { data_length_ = l;  }
  static long CalculateHashCode(IndexType& index);
  virtual size_t memory_size() {
    size_t temp = sizeof(*this);
    if (data_) {
      temp += data_->memory_size();
    }
    if (index_data_) {
      temp += index_data_->memory_size();
    }
    return temp;
  }

 protected:
  explicit AppDataCompressedScalarArray(
      const nimbus::GeometricRegion &global_reg,
      const nimbus::GeometricRegion &ob_reg,
      const int ghost_width);

  virtual nimbus::AppVar *CreateNew(const nimbus::GeometricRegion &ob_reg) const;

  // The data should be DataCompressedScalarArray (corresponding nimbus type).
  virtual void ReadAppData(const nimbus::DataArray &read_set,
                           const nimbus::GeometricRegion &read_reg);
  virtual void WriteAppData(const nimbus::DataArray &write_set,
                              const nimbus::GeometricRegion &write_reg) const;

 private:
  nimbus::GeometricRegion global_region_;
  nimbus::GeometricRegion local_region_;
  int ghost_width_;
  nimbus::Coord shift_;
  DataType* data_;
  // Index data should be external. This object does not have ownership.
  IndexType* index_data_;
  nimbus::int_dimension_t data_length_;
}; // class AppDataCompressedScalarArray

} // namespace application

#endif // NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_COMPRESSED_SCALAR_ARRAY_H_
