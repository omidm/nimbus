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

    virtual typename Vertex<T, key_t>::Map* vertices();

    virtual typename Vertex<T, key_t>::Iter begin();
    virtual typename Vertex<T, key_t>::ConstIter begin() const;
    virtual typename Vertex<T, key_t>::Iter end();
    virtual typename Vertex<T, key_t>::ConstIter end() const;

    virtual bool AddVertex(key_t key, T* entry);

    virtual bool AddVertex(key_t key, T* entry, Vertex<T, key_t>** vertex);

    virtual bool HasVertex(key_t key);

    virtual bool HasVertex(key_t key, typename Vertex<T, key_t>::Iter *iter);

    virtual bool GetVertex(key_t key, Vertex<T, key_t>** vertex);

    virtual bool RemoveVertex(key_t key);

    virtual bool AddEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end);

    virtual bool AddEdge(key_t start_key, key_t end_key);

    virtual bool AddEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end, Edge<T, key_t>** edge);

    virtual bool AddEdge(key_t start_key, key_t end_key, Edge<T, key_t>** edge);

    virtual bool RemoveEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end);

    virtual bool RemoveEdge(key_t start_key, key_t end_key);

  private:
    typename Vertex<T, key_t>::Map vertices_;
};


template<typename T, typename key_t>
Graph<T, key_t>::Graph() {
}

template<typename T, typename key_t>
Graph<T, key_t>::~Graph() {
  typename Vertex<T, key_t>::Iter iter;

  for (iter = vertices_.begin(); iter != vertices_.end(); ++iter) {
    typename Edge<T, key_t>::Iter it;
    typename Edge<T, key_t>::Map* edges = iter->second->outgoing_edges();
    for (it = edges->begin(); it != edges->end(); ++it) {
      delete (it->second);
    }
    delete (iter->second);
  }
}

template<typename T, typename key_t>
typename Vertex<T, key_t>::Map* Graph<T, key_t>::vertices() {
  return &vertices_;
}

template<typename T, typename key_t>
typename Vertex<T, key_t>::Iter Graph<T, key_t>::begin() {
  return vertices_.begin();
}

template<typename T, typename key_t>
typename Vertex<T, key_t>::ConstIter Graph<T, key_t>::begin() const {
  return vertices_.begin();
}

template<typename T, typename key_t>
typename Vertex<T, key_t>::Iter Graph<T, key_t>::end() {
  return vertices_.end();
}

template<typename T, typename key_t>
typename Vertex<T, key_t>::ConstIter Graph<T, key_t>::end() const {
  return vertices_.end();
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
bool Graph<T, key_t>::HasVertex(key_t key, typename Vertex<T, key_t>::Iter *iter) {
  *iter =  vertices_.find(key);
  return ((*iter) != vertices_.end());
}


template<typename T, typename key_t>
bool Graph<T, key_t>::GetVertex(key_t key, Vertex<T, key_t>** vertex) {
  typename Vertex<T, key_t>::Iter iter;
  if (!HasVertex(key, &iter)) {
    *vertex = NULL;
    return false;
  }

  *vertex = iter->second;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::RemoveVertex(key_t key) {
  typename Vertex<T, key_t>::Iter iter;
  if (!HasVertex(key, &iter)) {
    dbg(DBG_WARN, "WARNING: vertex with id %lu does not exist.\n", key);
    dbg(DBG_WARN, "WARNING: Nothing removed from graph.\n");
    return false;
  }

  Vertex<T, key_t>* vertex = iter->second;

  // first remove all the edges
  typename Edge<T, key_t>::Iter it;
  typename Edge<T, key_t>::Iter temp_it;

  typename Edge<T, key_t>::Map* outgoing_edges = vertex->outgoing_edges();
  for (it = outgoing_edges->begin(); it != outgoing_edges->end();) {
    temp_it = it++;
    RemoveEdge(temp_it->second->start_vertex(), temp_it->second->end_vertex());
  }

  typename Edge<T, key_t>::Map* incoming_edges = vertex->incoming_edges();
  for (it = incoming_edges->begin(); it != incoming_edges->end();) {
    temp_it = it++;
    RemoveEdge(temp_it->second->start_vertex(), temp_it->second->end_vertex());
  }

  vertices_.erase(iter);
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
  typename Vertex<T, key_t>::Iter start_iter;
  if (!HasVertex(start_key, &start_iter)) {
    dbg(DBG_ERROR, "ERROR: start vertex with id %lu does not exist.\n", start_key);
    return false;
  }
  typename Vertex<T, key_t>::Iter end_iter;
  if (!HasVertex(end_key, &end_iter)) {
    dbg(DBG_ERROR, "ERROR: end vertex with id %lu does not exist.\n", end_key);
    return false;
  }
  Vertex<T, key_t>* start = start_iter->second;
  Vertex<T, key_t>* end = end_iter->second;

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
bool Graph<T, key_t>::AddEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end,
    Edge<T, key_t>** edge) {
  if (start->HasOutgoingEdgeTo(end)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu already exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* new_edge = new Edge<T, key_t>(start, end);
  start->AddOutgoingEdge(new_edge);
  end->AddIncomingEdge(new_edge);
  *edge = new_edge;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::AddEdge(key_t start_key, key_t end_key, Edge<T, key_t>** edge) {
  typename Vertex<T, key_t>::Iter start_iter;
  if (!HasVertex(start_key, &start_iter)) {
    dbg(DBG_ERROR, "ERROR: start vertex with id %lu does not exist.\n", start_key);
    return false;
  }
  typename Vertex<T, key_t>::Iter end_iter;
  if (!HasVertex(end_key, &end_iter)) {
    dbg(DBG_ERROR, "ERROR: end vertex with id %lu does not exist.\n", end_key);
    return false;
  }
  Vertex<T, key_t>* start = start_iter->second;
  Vertex<T, key_t>* end = end_iter->second;

  if (start->HasOutgoingEdgeTo(end)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu already exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* new_edge = new Edge<T, key_t>(start, end);
  start->AddOutgoingEdge(new_edge);
  end->AddIncomingEdge(new_edge);
  *edge = new_edge;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::RemoveEdge(Vertex<T, key_t>* start, Vertex<T, key_t>* end) {
  typename Edge<T, key_t>::Iter iter;
  if (!start->HasOutgoingEdgeTo(end, &iter)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu does not exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* to_delete_edge = iter->second;
  start->RemoveOutgoingEdge(to_delete_edge);
  end->RemoveIncomingEdge(to_delete_edge);
  delete to_delete_edge;
  return true;
}

template<typename T, typename key_t>
bool Graph<T, key_t>::RemoveEdge(key_t start_key, key_t end_key) {
  typename Vertex<T, key_t>::Iter start_iter;
  if (!HasVertex(start_key, &start_iter)) {
    dbg(DBG_ERROR, "ERROR: start vertex with id %lu does not exist.\n", start_key);
    return false;
  }
  typename Vertex<T, key_t>::Iter end_iter;
  if (!HasVertex(end_key, &end_iter)) {
    dbg(DBG_ERROR, "ERROR: end vertex with id %lu does not exist.\n", end_key);
    return false;
  }
  Vertex<T, key_t>* start = start_iter->second;
  Vertex<T, key_t>* end = end_iter->second;

  typename Edge<T, key_t>::Iter iter;
  if (!start->HasOutgoingEdgeTo(end, &iter)) {
    dbg(DBG_WARN, "WARNING: edge from %lu to %lu does not exist.\n", start->id(), end->id());
    return true;
  }

  Edge<T, key_t>* to_delete_edge = iter->second;
  start->RemoveOutgoingEdge(to_delete_edge);
  end->RemoveIncomingEdge(to_delete_edge);
  delete to_delete_edge;
  return true;
}

}  // namespace nimbus

#endif  // NIMBUS_SHARED_GRAPH_H_
