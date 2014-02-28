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
  * The graph edge abstraction class, build as template class so that one can
  * choose the data structure to keep the data at each vertex dynamically.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_EDGE_H_
#define NIMBUS_SHARED_EDGE_H_

#include <sstream> // NOLINT
#include <iostream> // NOLINT
#include <string>
#include <vector>
#include <map>
#include <set>
#include "shared/nimbus_types.h"



namespace nimbus {

template<typename T, typename key_t>
class Vertex;

template<typename T, typename key_t>
class Edge {
  public:
    typedef typename std::map<key_t, Edge<T, key_t>*> Map;
    typedef typename std::map<key_t, Edge<T, key_t>*>::iterator Iter;
    typedef typename std::map<key_t, Edge<T, key_t>*>::const_iterator ConstIter;

    Edge(Vertex<T, key_t>* start_vertex, Vertex<T, key_t>* end_vertex);
    Edge(const Edge<T, key_t>& other);
    virtual ~Edge();

    virtual Vertex<T, key_t>* start_vertex();
    virtual Vertex<T, key_t>* end_vertex();

    Edge<T, key_t>& operator=(const Edge<T, key_t>& other);

  private:
    Vertex<T, key_t>* start_vertex_;
    Vertex<T, key_t>* end_vertex_;
};


template<typename T, typename key_t>
Edge<T, key_t>::Edge(Vertex<T, key_t>* start_vertex, Vertex<T, key_t>* end_vertex)
  :start_vertex_(start_vertex), end_vertex_(end_vertex) {
}

template<typename T, typename key_t>
Edge<T, key_t>::Edge(const Edge<T, key_t>& other) {
  start_vertex_ = other.start_vertex_;
  end_vertex_ = other.end_vertex_;
}

template<typename T, typename key_t>
Edge<T, key_t>::~Edge() {
}

template<typename T, typename key_t>
Vertex<T, key_t>* Edge<T, key_t>::start_vertex() {
  return start_vertex_;
}

template<typename T, typename key_t>
Vertex<T, key_t>* Edge<T, key_t>::end_vertex() {
  return end_vertex_;
}

template<typename T, typename key_t>
Edge<T, key_t>& Edge<T, key_t>::operator=(const Edge<T, key_t>& other) {
  start_vertex_ = other.start_vertex_;
  end_vertex_ = other.end_vertex_;
  return *this;
}


}  // namespace nimbus

#endif  // NIMBUS_SHARED_EDGE_H_
