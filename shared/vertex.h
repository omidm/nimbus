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
  * The graph vertex abstraction class, build as template class so that one can
  * choose the data structure to keep the data at each vertex dynamically.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_VERTEX_H_
#define NIMBUS_SHARED_VERTEX_H_

#include <sstream> // NOLINT
#include <iostream> // NOLINT
#include <string>
#include <vector>
#include <map>
#include <set>
#include "shared/nimbus_types.h"



namespace nimbus {

template<typename T, typename key_t>
class Graph;

template<typename T, typename key_t>
class Edge;

template<typename T, typename key_t>
class Vertex {
  friend class Graph<T, key_t>;
  public:
    explicit Vertex(T* entry) {}
    Vertex(const Vertex<T, key_t>& other) {}
    virtual ~Vertex() {}

    virtual std::map<key_t, Edge<T, key_t>*>* outgoing_edges() {}
    virtual std::map<key_t, Edge<T, key_t>*>* incoming_edges() {}
    virtual key_t id() {}
    virtual T* entry() {}

    Vertex<T, key_t>& operator=(const Vertex<T, key_t>& other) {}

  protected:
    virtual void AddOutgoingEdge(Edge<T, key_t>* e) {}

    virtual void AddIncomingEdge(Edge<T, key_t>* e) {}

    virtual void RemoveOutgoingEdge(Edge<T, key_t>* e) {}

    virtual void RemoveIncomingEdge(Edge<T, key_t>* e) {}

  private:
    std::map<key_t, Edge<T, key_t>*> outgoing_edges_;
    std::map<key_t, Edge<T, key_t>*> incoming_edges_;
    key_t id_;
    T* entry_;
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_VERTEX_H_
