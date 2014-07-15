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

#include <boost/tokenizer.hpp>
#include "shared/idset.h"

using namespace nimbus; // NOLINT
using boost::tokenizer;
using boost::char_separator;

// Following two functions are temporarily added, since we are in transition of
// changing the parsers and constructors. Until then, these two would be used
// in the constructor of the IDSet. Ideally the constructor should not call
// parser, which is the case in the coming implementation. I could not include
// these two from parser.h due to cyclic #include problem.

using boost::tokenizer;
using boost::char_separator;


// TODO(omidm): make sure that these are not used anymore and then remove.

// void tempParseIDSetFromString(const std::string& input,
//     IDSet<uint64_t>::IDSetContainer& list) {
//   int num;
//   std::string str = input.substr(1, input.length() - 2);
//   char_separator<char> separator(",");
//   tokenizer<char_separator<char> > tokens(str, separator);
//   tokenizer<char_separator<char> >::iterator iter = tokens.begin();
//   for (; iter != tokens.end(); ++iter) {
//     std::stringstream ss(*iter);
//     ss >> num;
//     list.push_back(num);
//   }
// }
//
// void tempParseIDSetFromString(const std::string& input,
//     IDSet<uint32_t>::IDSetContainer& list) {
//   int num;
//   std::string str = input.substr(1, input.length() - 2);
//   char_separator<char> separator(",");
//   tokenizer<char_separator<char> > tokens(str, separator);
//   tokenizer<char_separator<char> >::iterator iter = tokens.begin();
//   for (; iter != tokens.end(); ++iter) {
//     std::stringstream ss(*iter);
//     ss >> num;
//     list.push_back(num);
//   }
// }
//
// template<typename T>
// IDSet<T>::IDSet(std::string s) {
//   tempParseIDSetFromString(s, identifiers_);
// }


template<typename T>
IDSet<T>::IDSet() {}

template<typename T>
IDSet<T>::IDSet(const IDSetContainer& ids) {
  identifiers_ = ids;
}

template<typename T>
IDSet<T>::IDSet(const IDSet<T>& other)
: identifiers_(other.identifiers_) {
}

template<typename T>
IDSet<T>::~IDSet() {}

template<typename T>
bool IDSet<T>::Parse(const std::string& input) {
  T num;
  identifiers_.clear();
  if (input[0] != '{' || input[input.length() - 1] != '}') {
    std::cout << "ERROR: wrong format for IDSet." << std::endl;
    return false;
  }
  std::string str = input.substr(1, input.length() - 2);
  char_separator<char> separator(",");
  tokenizer<char_separator<char> > tokens(str, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (; iter != tokens.end(); ++iter) {
    std::stringstream ss(*iter);
    ss >> num;
    if (ss.fail()) {
      std::cout << "ERROR: wrong element in IDSet." << std::endl;
      identifiers_.clear();
      return false;
    }
    // identifiers_.push_back(num);
    identifiers_.insert(num);
  }
  return true;
}

template<typename T>
std::string IDSet<T>::toString() {
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

template<typename T>
void IDSet<T>::insert(T n) {
  // IDSetIter iter =  identifiers_.begin();
  // for (; iter !=  identifiers_.end(); ++iter) {
  //   if (*iter == n) {
  //     return;
  //   }
  // }
  // identifiers_.push_back(n);
  identifiers_.insert(n);
}

template<typename T>
void IDSet<T>::insert(const IDSet<T>& add_set) {
  IDSetContainer ids =  add_set.identifiers_;
  IDSetIter iter =  ids.begin();
  for (; iter !=  ids.end(); ++iter) {
    this->insert(*iter);
  }
}

template<typename T>
void IDSet<T>::remove(T n) {
  // IDSetIter iter =  identifiers_.begin();
  // for (; iter !=  identifiers_.end(); ++iter) {
  //   if (*iter == n) {
  //     identifiers_.erase(iter);
  //     break;
  //   }
  // }
  identifiers_.erase(n);
}

template<typename T>
void IDSet<T>::remove(IDSetIter it) {
  identifiers_.erase(it);
}

template<typename T>
void IDSet<T>::remove(const IDSet<T>& remove_set) {
  IDSetContainer ids =  remove_set.identifiers_;
  IDSetIter iter =  ids.begin();
  for (; iter !=  ids.end(); ++iter) {
    this->remove(*iter);
  }
}

template<typename T>
void IDSet<T>::clear() {
  identifiers_.clear();
}

template<typename T>
bool IDSet<T>::contains(T n) const {
  // ConstIter iter =  identifiers_.begin();
  // for (; iter !=  identifiers_.end(); ++iter) {
  //   if (*iter == n)
  //     return true;
  // }
  // return false;
  return (identifiers_.count(n) != 0);
}

template<typename T>
int IDSet<T>::size() const {
  return identifiers_.size();
}

template<typename T> void IDSet<T>::swap(IDSet<T>& idset) {
  identifiers_.swap(idset.identifiers_);
}

template<typename T>
typename IDSet<T>::IDSetIter IDSet<T>::begin() {
  return identifiers_.begin();
}

template<typename T>
typename IDSet<T>::IDSetIter IDSet<T>::end() {
  return identifiers_.end();
}

template<typename T>
typename IDSet<T>::ConstIter IDSet<T>::begin() const {
  return identifiers_.begin();
}

template<typename T>
typename IDSet<T>::ConstIter IDSet<T>::end() const {
  return identifiers_.end();
}

template<typename T>
IDSet<T>& IDSet<T>::operator= (const IDSet<T>& right) {
  identifiers_ = right.identifiers_;
  return *this;
}

template class IDSet<uint64_t>;
template class IDSet<uint32_t>;



