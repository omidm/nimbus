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

#ifndef NIMBUS_SHARED_IDSET_H_
#define NIMBUS_SHARED_IDSET_H_

#include <boost/tokenizer.hpp>
#include <boost/unordered_set.hpp>
#include <google/protobuf/repeated_field.h>
#include <algorithm>
#include <sstream> // NOLINT
#include <iostream> // NOLINT
#include <string>
#include <list>
#include <set>
#include "shared/nimbus_types.h"

namespace nimbus {

template<typename T>
class IDSet {
 public:
  // typedef typename std::list<T> IDSetContainer;
  // typedef typename std::list<T>::iterator IDSetIter;
  // typedef typename std::list<T>::const_iterator ConstIter;
  typedef typename boost::unordered_set<T> IDSetContainer;
  typedef typename boost::unordered_set<T>::iterator IDSetIter;
  typedef typename boost::unordered_set<T>::const_iterator ConstIter;

  IDSet();
  explicit IDSet(const IDSetContainer& ids);
  IDSet(const IDSet<T>& other);
  virtual ~IDSet();

  // TODO(omidm): remove this obsolete constructor.
  // explicit IDSet(std::string s);

  bool Parse(const std::string& input);
  virtual std::string ToNetworkData();
  virtual void insert(T entry);
  virtual void insert(const IDSet<T>& add_set);
  virtual void remove(T entry);
  virtual void remove(IDSetIter it);
  virtual void remove(const IDSet<T>& remove_set);
  virtual void clear();
  virtual bool contains(T entry) const;
  virtual int size() const;
  virtual void swap(IDSet<T>& idset);

  IDSetIter begin();
  IDSetIter end();

  ConstIter begin() const;
  ConstIter end() const;

  IDSet<T>& operator= (const IDSet<T>& right);

  virtual void ConvertToRepeatedField(google::protobuf::RepeatedField<T>* b);
  virtual void ConvertFromRepeatedField(const google::protobuf::RepeatedField<T>& b);

 private:
  IDSetContainer identifiers_;
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_IDSET_H_
