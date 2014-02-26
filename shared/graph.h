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
  * The graph graph abstraction class, build as template class so that one can
  * choose the data structure to keep the data at each vertex dynamically.
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
#include "shared/dbg.h"
#include "shared/nimbus_types.h"
#include "shared/edge.h"
#include "shared/vertex.h"

namespace nimbus {

template<typename T, typename key_t>
class Graph {
  public:
    Graph();
    virtual ~Graph();

    virtual std::map<key_t, Vertex<T, key_t>*>* vertices();

    virtual bool AddVertex(key_t key, T* entry);

    virtual bool AddVertex(key_t key, T* entry, Vertex<T, key_t>** vertex);

    virtual bool HasVertex(key_t key);

    virtual bool GetVertex(key_t key, Vertex<T, key_t>** vertex);

    virtual bool RemoveVertex(key_t key);

    virtual bool AddEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end);

    virtual bool AddEdge(key_t start_key, key_t end_key);

    virtual bool RemoveEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end);

    virtual bool RemoveEdge(key_t start_key, key_t end_key);

  private:
    std::map<key_t, Vertex<T, key_t>*> vertices_;
};


template<typename T, typename key_t>
Graph<T, key_t>::Graph() {
}

template<typename T, typename key_t>
Graph<T, key_t>::~Graph() {
  typename std::map<key_t, Vertex<T, key_t>*>::iterator iter;

  for (iter = vertices_.begin(); iter != vertices_.end(); ++iter) {
    typename std::map<key_t, Edge<T, key_t>*>::iterator it;
    std::map<key_t, Edge<T, key_t>*>* edges = iter->second->outgoing_edges();
    for (it = edges->begin(); it != edges->end(); ++it) {
      delete (it->second);
    }
    delete (iter->second);
  }
}

template<typename T, typename key_t>
std::map<key_t, Vertex<T, key_t>*>* Graph<T, key_t>::vertices() {
  return &vertices_;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::AddVertex(key_t key, T* entry) {
  if (HasVertex(key)) {
    dbg(DBG_ERROR, "ERROR: vertex with id %lu already exist.\n", key);
    return false;
  }

  Vertex<T, key_t>* new_vertex = new Vertex<T, key_t>(key, entry);
  vertices_[key] = new_vertex;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::AddVertex(key_t key, T* entry, Vertex<T, key_t>** vertex) {
  if (HasVertex(key)) {
    dbg(DBG_ERROR, "ERROR: vertex with id %lu already exist.\n", key);
    *vertex = NULL;
    return false;
  }

  Vertex<T, key_t>* new_vertex = new Vertex<T, key_t>(key, entry);
  vertices_[key] = new_vertex;
  *vertex = new_vertex;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::HasVertex(key_t key) {
  return (vertices_.find(key) != vertices_.end());
}

template<typename T, typename key_t>
bool Graph<T, key_t>::GetVertex(key_t key, Vertex<T, key_t>** vertex) {
  if (!HasVertex(key)) {
    dbg(DBG_ERROR, "ERROR: vertex with id %lu does not exist.\n", key);
    *vertex = NULL;
    return false;
  }

  *vertex = vertices_[key];
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::RemoveVertex(key_t key) {
  if (!HasVertex(key)) {
    dbg(DBG_WARN, "WARNING: vertex with id %lu does not exist.\n", key);
    return true;
  }

  Vertex<T, key_t>* vertex = vertices_[key];

  // first remove all the edges
  typename std::map<key_t, Edge<T, key_t>*>::iterator it;

  std::map<key_t, Edge<T, key_t>*>* outgoing_edges = vertex->outgoing_edges();
  for (it = outgoing_edges->begin(); it != outgoing_edges->end(); ++it) {
    RemoveEdge(it->second->start_vertex(), it->second->end_vertex());
  }

  std::map<key_t, Edge<T, key_t>*>* incoming_edges = vertex->incoming_edges();
  for (it = incoming_edges->begin(); it != incoming_edges->end(); ++it) {
    RemoveEdge(it->second->start_vertex(), it->second->end_vertex());
  }

  vertices_.erase(key);
  delete vertex;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::AddEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end) {
  if (start->HasOutgoingEdgeTo(end)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu already exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* new_edge = new Edge<T, key_t>(start, end);
  start->AddOutgoingEdge(new_edge);
  end->AddIncomingEdge(new_edge);
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::AddEdge(key_t start_key, key_t end_key) {
  if (!HasVertex(start_key)) {
    dbg(DBG_ERROR, "ERROR: start vertex with id %lu does not exist.\n", start_key);
    return false;
  }
  if (!HasVertex(end_key)) {
    dbg(DBG_ERROR, "ERROR: end vertex with id %lu does not exist.\n", end_key);
    return false;
  }
  Vertex<T, key_t>* start = vertices_[start_key];
  Vertex<T, key_t>* end = vertices_[end_key];

  if (start->HasOutgoingEdgeTo(end)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu already exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* new_edge = new Edge<T, key_t>(start, end);
  start->AddOutgoingEdge(new_edge);
  end->AddIncomingEdge(new_edge);
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::RemoveEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end) {
  if (start->HasOutgoingEdgeTo(end)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu does not exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* to_delete_edge =
    start->outgoing_edges()->operator[](end->id());
  start->RemoveOutgoingEdge(to_delete_edge);
  end->RemoveIncomingEdge(to_delete_edge);
  delete to_delete_edge;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::RemoveEdge(key_t start_key, key_t end_key) {
  if (!HasVertex(start_key)) {
    dbg(DBG_ERROR, "ERROR: start vertex with id %lu does not exist.\n", start_key);
    return false;
  }
  if (!HasVertex(end_key)) {
    dbg(DBG_ERROR, "ERROR: end vertex with id %lu does not exist.\n", end_key);
    return false;
  }
  Vertex<T, key_t>* start = vertices_[start_key];
  Vertex<T, key_t>* end = vertices_[end_key];

  if (start->HasOutgoingEdgeTo(end)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu does not exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* to_delete_edge =
    start->outgoing_edges()->operator[](end->id());
  start->RemoveOutgoingEdge(to_delete_edge);
  end->RemoveIncomingEdge(to_delete_edge);
  delete to_delete_edge;
  return true;
}








}  // namespace nimbus

#endif  // NIMBUS_SHARED_GRAPH_H_
