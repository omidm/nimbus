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

#include "applications/physbam/water//projection/translator_util.h"
#include "src/data/physbam/physbam_data.h"
#include "src/shared/nimbus.h"

#include "applications/physbam/water//projection/data_raw_array_m2c.h"

namespace application {

DataRawArrayM2C::DataRawArrayM2C(std::string name) {
  set_name(name);
}

nimbus::Data* DataRawArrayM2C::Clone() {
  return (new DataRawArrayM2C(name()));
}

bool DataRawArrayM2C::SaveToNimbus(const PhysBAM::ARRAY<TV_INT>& array_input) {
  Header header;
  header.n = PhysBAM::Value(array_input.m);
  ClearTempBuffer();
  AddToTempBuffer(reinterpret_cast<char*>(&header), sizeof(header));
  Buffer buffer;
  SerializePhysBAMArray(array_input, &buffer);
  AddToTempBuffer(reinterpret_cast<char*>(buffer.pointer), buffer.size);
  assert(buffer.size == header.n * sizeof(TV_INT));
  buffer.Clean();
  CommitTempBuffer();
  return true;
}

bool DataRawArrayM2C::LoadFromNimbus(PhysBAM::ARRAY<TV_INT>* array) {
  assert(array != NULL);
  char* pointer = buffer();
  assert(pointer != NULL);
  const Header &header = *(reinterpret_cast<const Header*>(pointer));
  array->m = header.n;
  array->buffer_size = array->m;
  pointer += sizeof(Header);
  Buffer buffer;
  buffer.pointer = reinterpret_cast<void*>(pointer);
  buffer.size = header.n * sizeof(TV_INT);
  DeserializePhysBAMArray(buffer, array);
  buffer.Reset();
  return true;
}

float DataRawArrayM2C::FloatingHash() {
  const int *b = (int *)buffer();
  size_t s = size();
  int sum = 0;
  for (size_t i = 0; i < (s/sizeof(int)); ++i) {
    sum += b[i];
  }
  return(float(sum));
}

}  // namespace application
