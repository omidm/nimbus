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
  * Object representation of an identifires.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/id.h"

using namespace nimbus; // NOLINT


template<typename T>
ID<T>::ID() {
  elem_ = 0;
}

template<typename T>
ID<T>::ID(const T& elem) {
  elem_ = elem;
}

template<typename T>
ID<T>::ID(const ID<T>& other)
: elem_(other.elem_) {
}

template<typename T>
ID<T>::~ID() {}

template<typename T>
bool ID<T>::Parse(const std::string& input) {
  std::stringstream ss(input);
  T num;
  ss >> num;
  if (ss.fail()) {
    std::cout << "ERROR: wrong element as ID." << std::endl;
    return false;
  }
  elem_ = num;
  return true;
}

template<typename T>
std::string ID<T>::toString() {
  std::ostringstream ss;
  ss << elem_;
  std::string rval;
  return ss.str();
}

template<typename T>
void ID<T>::set_elem(T elem) {
  elem_ = elem;
}

template<typename T>
T ID<T>::elem() {
  return elem_;
}

template<typename T>
ID<T>& ID<T>::operator= (const ID<T>& right) {
  elem_ = right.elem_;
  return *this;
}

template class ID<uint64_t>;
template class ID<uint32_t>;
template class ID<int32_t>;



