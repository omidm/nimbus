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
  * Utility functions to serialze and deserialze PhysBAM array and vector. Most
  * PhysBAM data structures are ultimately serialized to a linear array. By
  * providing functionality to serialze the linear array, Nimbus can pass
  * PhysBAM data around with little translation cost compared to translators in
  * data directory.
  *
  * The implementation here operates at a low level, and thus it cannot
  * understand geometry compared to the implentation in data directory. This is
  * the case when Nimbus performs projection calculation.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_APPLICATION_SMOKE_PROJECTION_TRANSLATOR_UTIL_H_
#define NIMBUS_APPLICATION_SMOKE_PROJECTION_TRANSLATOR_UTIL_H_

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>

#include "shared/nimbus.h"

// An internal data structure used to represent a buffer.
struct Buffer {
  void* pointer;
  size_t size;
  Buffer() : pointer(NULL), size(0) {}
  void Reset() {
    pointer = NULL;
    size = 0;
  }
  void Clean() {
    if (pointer) {
      free(pointer);
    }
    Reset();
  }
};

// Serializes the PhysBAM array into a buffer. Allocates memory in the buffer.
template<class T> bool SerializePhysBAMArray(
    const PhysBAM::ARRAY<T>& array_input,
    Buffer* buffer) {
  assert(buffer != NULL);
  int64_t length = array_input.m;
  assert(buffer->pointer == NULL);
  buffer->pointer = malloc(length * sizeof(T));
  buffer->size = length * sizeof(T);
  T* pointer = reinterpret_cast<T*>(buffer->pointer);
  memcpy(pointer, array_input.base_pointer, length * sizeof(T));
  return true;
}

// Deserialzes the data in the buffer into the PhysBAM array. Allocates memory
// for the PhysBAM array.
template<class T> bool DeserializePhysBAMArray(
    const Buffer& buffer_input,
    PhysBAM::ARRAY<T>* array) {
  assert(array != NULL);
  array->Clean_Memory();
  int64_t length = buffer_input.size / sizeof(T);
  array->m = length;
  array->base_pointer = new T[length];
  const T* pointer = reinterpret_cast<T*>(buffer_input.pointer);
  memcpy(array->base_pointer, pointer, length * sizeof(T));
  return true;
}

// Serializes the PhysBAM vector into a buffer. Allocates memory in the buffer.
template<class T> bool SerializePhysBAMVector(
    const PhysBAM::VECTOR_ND<T>& vector_input,
    Buffer* buffer) {
  assert(buffer != NULL);
  int64_t length = vector_input.n;
  assert(buffer->pointer == NULL);
  buffer->pointer = malloc(length * sizeof(T));
  buffer->size = length * sizeof(T);
  T* pointer = reinterpret_cast<T*>(buffer->pointer);
  memcpy(pointer, vector_input.x, length * sizeof(T));
  return true;
}

// Deserialzes the data in the buffer into the PhysBAM vector. Allocates memory
// for the PhysBAM vector.
template<class T> bool DeserializePhysBAMVector(
    const Buffer& buffer_input,
    PhysBAM::VECTOR_ND<T>* vector) {
  assert(vector != NULL);
  assert(vector->Owns_Data());
  if (vector->x != NULL && vector->Owns_Data()) {
    delete[] vector->x;
  }
  int64_t length = buffer_input.size / sizeof(T);
  vector->n = length;
  vector->x = new T[length];
  const T* pointer = reinterpret_cast<T*>(buffer_input.pointer);
  memcpy(vector->x, pointer, length * sizeof(T));
  return true;
}

#endif  // NIMBUS_APPLICATION_SMOKE_PROJECTION_TRANSLATOR_UTIL_H_
