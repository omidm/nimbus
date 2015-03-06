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
 * Nimbus data type for "cell_index_to_matrix_index", which is a mapping from
 * cell index to matrix index in projection calculation.
 * Translating functionality is implemented in the class.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_DATA_RAW_GRID_ARRAY_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_DATA_RAW_GRID_ARRAY_H_

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include "data/physbam/physbam_data.h"
#include "shared/nimbus.h"

namespace application {

class DataRawGridArray : public nimbus::PhysBAMData {
 private:
  struct Header {
    int64_t n;
  };
 public:
  typedef int T;
  typedef PhysBAM::VECTOR<int, 3> TV_INT;
  explicit DataRawGridArray(std::string name);
  virtual nimbus::Data* Clone();
  float FloatingHash();

  // Saves the PhysBAM data structure to this nimbus data instance.
  bool SaveToNimbus(const PhysBAM::ARRAY<T, TV_INT>& array_input);
  // Loads the PhysBAM data structure from this nimbus data instance.
  bool LoadFromNimbus(PhysBAM::ARRAY<T, TV_INT>* array);
};

}  // namespace application

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_DATA_RAW_GRID_ARRAY_H_
