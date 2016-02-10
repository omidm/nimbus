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

#include "src/data/physbam/physbam_data.h"
#include "src/shared/nimbus.h"
#include "string.h"
#include "applications/physbam/water//data_compressed_scalar_array.h"

namespace application {

template<typename T> DataCompressedScalarArray<T>::DataCompressedScalarArray(
    std::string name) {
  set_name(name);
}

template<typename T> nimbus::Data* DataCompressedScalarArray<T>::Clone() {
  return (new DataCompressedScalarArray<T>(name()));
}

template<typename T> void DataCompressedScalarArray<T>::Create() {
  nimbus::PhysBAMDataWithMeta::Create();
}

template<typename T> float DataCompressedScalarArray<T>::FloatingHash() {
  const T *b = (T *)buffer();
  size_t s = size();
  T sum = 0;
  for (size_t i = 0; i < (s/sizeof(T)); ++i) {
    sum += b[i];
  }
  return(float(sum));
}

template class DataCompressedScalarArray<float>;
template class DataCompressedScalarArray<int>;

} // namespace application
