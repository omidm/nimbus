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
  * The graph abstraction which defines the vertex, edge and graph class. They
  * are all build as template class so that one can choose the data structure
  * to keep the data at each vertex dynamically.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_GRAPH_H_
#define NIMBUS_SHARED_GRAPH_H_

#include <sstream> // NOLINT
#include <iostream> // NOLINT
#include <string>
#include <vector>
#include <map>
#include <set>
#include "shared/nimbus_types.h"



namespace nimbus {

template<typename T>
class Edge;

template<typename T>
class Vertex;

template<typename T>
class Vertex {
  friend class Graph;
  public:
    explicit Vertex(T* entry) {}
    Vertex(const Vertex<T>& other) {}
    virtual ~Vertex() {}

    virtual std::map<edge_id_t, Edge<T>*>* outgoing_edges() {}
    virtual std::map<edge_id_t, Edge<T>*>* incoming_edges() {}
    virtual T* entry() {}

    Vertex<T>& operator=(const Vertex<T>& other) {}

  protected:
    virtual void AddOutgoingEdge(Edge<T>* e) {}

    virtual void AddIncomingEdge(Edge<T>* e) {}

    virtual void RemoveOutgoingEdge(Edge<T>* e) {}

    virtual void RemoveIncomingEdge(Edge<T>* e) {}

  private:
    std::map<edge_id_t, Edge<T>*> outgoing_edges_;
    std::map<edge_id_t, Edge<T>*> incoming_edges_;
    T* entry_;
};

template<typename T>
class Edge {
  public:
    Edge(edge_id_t id, Vertex<T>* startVertex, Vertex<T>* endVertex) {}
    Edge(const Edge<T>& other) {}
    virtual ~Edge() {}

    virtual Vertex<T>* start_vertex() {}
    virtual Vertex<T>* end_vertex() {}
    virtual edge_id_t id();

    Edge<T>& operator=(const Edge<T>& other) {}

  private:
    Vertex<T>* start_vertex_;
    Vertex<T>* end_vertex_;
    edge_id_t id_;
};

template<typename T, typename key_t>
class Graph {
  public:
    Graph() {}
    virtual ~Graph() {}

    virtual std::map<key_t, Vertex<T>*>* vertices() {}
    virtual std::map<edge_id_t, Edge<T>*>* edges() {}

    virtual bool AddVertex(key_t key, T* entry) {}

    virtual bool AddVertex(key_t key, T* entry, Vertex<T>** vertex) {}

    virtual bool GetVertex(key_t key, Vertex<T>** vertex) {}

    virtual bool RemoveVertex(key_t key) {}

    virtual bool AddEdge(Vertex<T>* start, Vertex<T>* end) {}

    virtual bool AddEdge(key_t start_key, key_t end_key) {}

    virtual bool RemoveEdge(Vertex<T>* start, Vertex<T>* end) {}

    virtual bool RemoveEdge(key_t start_key, key_t end_key) {}

  private:
    std::map<key_t, Vertex<T>*> vertices_;
    std::map<edge_id_t, Edge<T>*> edges_;
};


}  // namespace nimbus

#endif  // NIMBUS_SHARED_GRAPH_H_
