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

#include "application/smoke/projection/translator_util.h"
#include "data/physbam/physbam_data.h"
#include "shared/nimbus.h"

#include "application/smoke/projection/data_sparse_matrix.h"

namespace application {

DataSparseMatrix::DataSparseMatrix(std::string name) {
  set_name(name);
}

nimbus::Data* DataSparseMatrix::Clone() {
  return (new DataSparseMatrix(name()));
}


bool DataSparseMatrix::SaveToNimbus(
    const PhysBAM::SPARSE_MATRIX_FLAT_NXN<float>& matrix) {
  Header header;
  header.n = matrix.n;
  header.length_offset = matrix.offsets.m;
  header.length_element = matrix.A.m;
  header.length_diagonal = matrix.diagonal_index.m;
  ClearTempBuffer();
  AddToTempBuffer(reinterpret_cast<char*>(&header), sizeof(header));
  Buffer buffer;
  SerializePhysBAMArray(matrix.offsets, &buffer);
  AddToTempBuffer(reinterpret_cast<char*>(buffer.pointer), buffer.size);
  buffer.Clean();
  SerializePhysBAMArray(matrix.A, &buffer);
  AddToTempBuffer(reinterpret_cast<char*>(buffer.pointer), buffer.size);
  buffer.Clean();
  SerializePhysBAMArray(matrix.diagonal_index, &buffer);
  AddToTempBuffer(reinterpret_cast<char*>(buffer.pointer), buffer.size);
  buffer.Clean();
  CommitTempBuffer();
  return true;
}

bool DataSparseMatrix::LoadFromNimbus(
    PhysBAM::SPARSE_MATRIX_FLAT_NXN<float>* matrix) {
  assert(matrix != NULL);
  char* pointer = buffer();
  assert(pointer != NULL);
  const Header &header = *(reinterpret_cast<const Header*>(pointer));
  matrix->n = header.n;
  assert((size_t)size() == sizeof(Header)
         + header.length_offset *
           sizeof(typename PhysBAM::ARRAY<int>::ELEMENT)
         + header.length_element *
           sizeof(typename PhysBAM::ARRAY<PhysBAM::SPARSE_MATRIX_ENTRY<float> >
                  ::ELEMENT)
	 + header.length_diagonal *
	   sizeof(typename PhysBAM::ARRAY<int>::ELEMENT));
  pointer += sizeof(Header);
  Buffer buffer;
  buffer.pointer = reinterpret_cast<void*>(pointer);
  buffer.size = header.length_offset *
      sizeof(typename PhysBAM::ARRAY<int>::ELEMENT);
  pointer += buffer.size;
  DeserializePhysBAMArray<int>(buffer, &matrix->offsets);
  buffer.Reset();

  buffer.pointer = reinterpret_cast<void*>(pointer);
  buffer.size = header.length_element *
      sizeof(typename PhysBAM::ARRAY<PhysBAM::SPARSE_MATRIX_ENTRY<float> >
             ::ELEMENT);
  pointer += buffer.size;
  DeserializePhysBAMArray<PhysBAM::SPARSE_MATRIX_ENTRY<float> >(
      buffer, &matrix->A);
  buffer.Reset();

  buffer.pointer = reinterpret_cast<void*>(pointer);
  buffer.size = header.length_diagonal * 
      sizeof(typename PhysBAM::ARRAY<int>::ELEMENT);
  pointer += buffer.size;
  DeserializePhysBAMArray<int>(buffer, &matrix->diagonal_index);
  buffer.Reset();

  // matrix->Initialize_Diagonal_Index();
  
  return true;
}

}  // namespace application
