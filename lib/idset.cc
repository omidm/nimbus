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
  * Object representation of a set of identifires.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "lib/idset.h"

using namespace nimbus; // NOLINT


IDSet::IDSet() {}

IDSet::IDSet(std::string s) {
  parseIDSetFromString(s, identifiers_);
}

IDSet::~IDSet() {}

std::string IDSet::toString() {
  bool empty = true;
  std::string rval = "{";
  IDSetIter iter =  identifiers_.begin();
  for (; iter !=  identifiers_.end(); ++iter) {
    empty = false;
    std::ostringstream ss;
    ss << *iter;
    rval += ss.str();
    rval += ",";
  }
  if (empty)
    rval += "}";
  else
    rval[rval.length() - 1] = '}';
  return rval;
}

void IDSet::insert(int n) {
  identifiers_.insert(n);
}

void IDSet::clear() {
  identifiers_.clear();
}

int IDSet::size() {
  return identifiers_.size();
}


IDSet::IDSetIter IDSet::begin() {
  return identifiers_.begin();
}

IDSet::IDSetIter IDSet::end() {
  return identifiers_.end();
}


// template<typename T>
// IDSet<T>::IDSet() {}
//
// template<typename T>
// IDSet<T>::IDSet(std::string s) {
//   parseIDSetFromString(s, identifiers_);
// }
//
// template<typename T>
// IDSet<T>::~IDSet() {}
//
// template<typename T>
// std::string IDSet<T>::toString() {
//   bool empty = true;
//   std::string rval = "{";
//   IDSetIter iter =  identifiers_.begin();
//   for (; iter !=  identifiers_.end(); ++iter) {
//     empty = false;
//     std::ostringstream ss;
//     ss << *iter;
//     rval += ss.str();
//     rval += ",";
//   }
//   if (empty)
//     rval += "}";
//   else
//     rval[rval.length() - 1] = '}';
//   return rval;
// }
//
// template<typename T>
// void IDSet<T>::insert(int n) {
//   identifiers_.insert(n);
// }
//
// template<typename T>
// void IDSet<T>::clear() {
//   identifiers_.clear();
// }
//
// template<typename T>
// int IDSet<T>::size() {
//   return identifiers_.size();
// }
//
//
// template<typename T>
// IDSet<T>::IDSetIter IDSet<T>::begin() {
//   return identifiers_.begin();
// }
//
// template<typename T>
// IDSet<T>::IDSetIter IDSet<>::end() {
//   return identifiers_.end();
// }


